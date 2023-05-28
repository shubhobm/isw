# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("data.table","harmonicmeanp","Matrix","RSpectra","sf","spdep")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

#parameters for approximate Hilbert embedding
evec_cnt = 20 #L in paper
invCDFnsamples = 5 #Dprime in paper

#load subset of Chicago crime data 
load("../data/chicago_theft2018.RData")
#data contents:
#beat_map --- map of beats (geographic area subdivision used by police)
#distros --- contains three matrices (Tuesday, Thursday, and Saturday) 
#each matrix row corresponds to a day in 2018, each column corresponds to a beat
#the values give the row-normalized counts of Theft crime for day & beat
#thus, each matrix row is a probability distribution on the beat graph

#adjacency matrix for the underlying graph
#graph nodes are beats, edges capture whether the beats share a boundary
A = nb2mat(poly2nb(as(beat_map, "Spatial")), style="B", zero.policy = TRUE)
#Laplacian matrix
D = Diagonal(x = rowSums(A))
Lap = D - A

#eigendecomposition
ev = eigs(Lap, k=evec_cnt+1, which = "SM")
#drop the constant eigenvector with zero eigenvalue & reorder
evals = ev$values[evec_cnt:1]
evecs = ev$vectors[,evec_cnt:1]

#approximate Hilbert embedding for a probability distribution on the graph
#a distribution is represented as a vector of weights with sum=1
#we use rolling joins to efficiently compute the quantiles
approx_hilbert_emb <- function(wts, Dprime, lambda, phi){
    L = ncol(phi) #number of eigenvectors
    alpha = 1/lambda**2 #biharmonic scaling
    
    #quantiles to be extracted, for each eigenvector
    q = CJ(ell=1:L, s_k=(1:Dprime)/(Dprime+1))
    setkey(q, "ell", "s_k")
    
    #compute pushforward by each eigenvector
    df = data.table(x=as.numeric(phi), wt=rep(wts, L), 
                    ell = rep(1:L, each=length(wts)))
    df = df[,.(wt=sum(wt)), keyby=.(ell, x)] #aggregate & sort
    
    #quantile definition has strict inequality, so subtract eps
    eps = 1e-8
    df = df[, s_k:=cumsum(wt)/sum(wt)-eps, by=ell]
    setkey(df, "ell", "s_k")
    
    #extract quantiles via rolling join for all eigenvectors at once
    df[q,,roll=-Inf][order(ell,s_k)][, x:=sqrt(alpha[ell]/Dprime)*x]$x
}

#approximate Hilbert embeddings for all distributions (=rows of the matrix)
embed_all_rows <- function(R, ...){
    do.call('rbind', lapply(seq_len(nrow(R)), 
            function(i) approx_hilbert_emb(R[i,], ...)))
}

combo_test <- function(E1, E2){
    #t-test on each dimension of the Hilbert embedding
    p_ttests = sapply(seq_len(ncol(E1)), function(i)
        tryCatch({t.test(E1[, i], E2[, i])$p.value}, error=function(e){NaN}))
    if(all(is.nan(p_ttests))){
        stop("data are essentially constant")
    } else {
        #NaN values provide no information towards acceptance/rejection
        p_ttests = p_ttests[!is.nan(p_ttests)]
        #overall p-value using harmonic mean combo
        unname(p.hmp(p_ttests, L=length(p_ttests), multilevel=FALSE))
    }
}

ETue = embed_all_rows(distros$Tuesday, invCDFnsamples, evals, evecs)
EThu = embed_all_rows(distros$Thursday, invCDFnsamples, evals, evecs)
ESat = embed_all_rows(distros$Saturday, invCDFnsamples, evals, evecs)

cat("Tuesday vs Thursday:", combo_test(ETue, EThu), "\n") #0.4518086 
cat("Tuesday vs Saturday:", combo_test(ETue, ESat), "\n") #4.728311e-06 