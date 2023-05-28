# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("data.table","parallel","harmonicmeanp")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

# ----------------------------------------------------------------------
# auxiliary functions
# ======================================================================

mat2dt <- function (x) {
  x = as.matrix(x)
  data.table(
    i = rep(1:nrow(x), ncol(x)), 
    j =  rep(1:ncol(x), each = nrow(x)), 
    x = as.vector(x)
  )
}

# calculate inverse CDF using histogram bins
# histogram is represented as $idx --- case/person id, $x -- the bin center locations, $wt the bin counts/weights
invCDF_data_wtd <- function(df, invCDFnsamples){
    levs = seq(0, 1 - 1e-16, length.out=invCDFnsamples) #want to capture 100th percentile too
    
    #just in case: aggregate weights by location
    df = df[,.(wt=sum(wt)), keyby=.(idx, x)][, lev:=cumsum(wt)/sum(wt), by=idx]
    q = CJ(idx=unique(df$idx), lev=levs)
    
    setkey(df, "idx", "lev")
    setkey(q, "idx", "lev")
    
    qq = df[q,,roll=-Inf][order(idx, lev)]
    
    M = matrix(qq$x, ncol=length(levs), byrow=TRUE)
    M
}

#input is in the dataframe form
pushForward = function(df, T=NULL, type="circle", ell=1, ...){
    # assign T and shift if null
    if(is.null(T)){
        T = max(df$x)
    }
    
    # push forward using eigenfunction
    if(type=="circle"){
        lambda = (2*pi*ell/T)^2
        
        Phi1 = copy(df)[, x:=sqrt(2/T)*cos(2*pi*ell*x/T)]
        Phi2 = copy(df)[, x:=sqrt(2/T)*sin(2*pi*ell*x/T)]
        
        # inverse map
        E = cbind(invCDF_data_wtd(Phi1, ...), invCDF_data_wtd(Phi2, ...))
    } else{
        lambda = (pi*ell/T)^2
        Phi = copy(df)[, x:= sqrt(2/T)*cos(pi*ell*x/T)]
        
        # inverse map
        E = invCDF_data_wtd(Phi, ...)
    }
    E
}

# generate truncated normal
rnormt = function(n, t=3){
    out = rnorm(n)
    at = abs(t)
    out[which(out > at)] = at
    out[which(out < -at)] = -at
    out
}

# ----------------------------------------------------------------------
# Simulation parameters and wrapper code
# ======================================================================
set.seed(10022020)
m = 30
locs = t = (1:m)/(m+1)
a0 = 1.2
a1 = 2.3
a2 = 4.2
mu = a0 + a1*cos(2*pi*t) + a2*cos(2*pi*t)
n1 = 60
n2 = 40
shift = 3*(1+sqrt(2)+sqrt(3))

get_pvals = function(seed, delta){
    set.seed(1e3*seed)
    X1 = matrix(0, n1, m)
    for(i in 1:n1){
        alpha = rnormt(m) + sqrt(2)*rnormt(m)*cos(2*pi*t) + sqrt(3)*rnormt(m)*sin(2*pi*t)
        X1[i,] = mu + alpha
    }
    X2 = matrix(0, n2, m)
    for(i in 1:n2){
        alpha = rnormt(m) + sqrt(2)*rnormt(m)*cos(2*pi*t) + sqrt(3)*rnormt(m)*sin(2*pi*t)
        X2[i,] = mu + alpha + delta
    }
    s1 = matrix(unlist(lapply(1:n1, function(i) rmultinom(1, size=1e3, prob=X1[i,]+shift))), nrow=n1, byrow=T)
    s2 = matrix(unlist(lapply(1:n2, function(i) rmultinom(1, size=1e3, prob=X2[i,]+shift))), nrow=n2, byrow=T)
    
    #data frame form with weights and locations
    df1 = mat2dt(s1)[order(i)][,idx:=i][, wt:=x][,x:=locs[j]][,c("i","j"):=NULL]
    df2 = mat2dt(s2)[order(i)][,idx:=i][, wt:=x][,x:=locs[j]][,c("i","j"):=NULL]  
    
    ############## SLICED VERSION ##############
    #notice how just keeping the first coordinate in sliced version gives almost same result as no slicing!!
    evaluate = function(L,D){
        E1 = Reduce(cbind, lapply(1:L, function(i) pushForward(df1, ell=i, type="line", invCDFnsamples=D)))
        E2 = Reduce(cbind, lapply(1:L, function(i) pushForward(df2, ell=i, type="line", invCDFnsamples=D)))
        
        #apply t-test to each dimension and combine
        p_ttests = sapply(seq_len(ncol(E1)), function(i)
            tryCatch({t.test(E1[, i], E2[, i])$p.value}, error=function(e){1.0}))
        
        #overall p-value using harmonic mean combo
        p_harmonic_sliced = p.hmp(p_ttests, L=length(p_ttests), multilevel=FALSE)
        p_harmonic_sliced
    }
    
    evaluate(3,10) # evaluate for L = 3, D' = 10
}

# ----------------------------------------------------------------------
# Run the experiments
# ======================================================================
# consider using a lower nrep for faster runs
nrep = 1e3
delta_list = seq(0,5,by=.1)
out_list = vector("list",length(delta_list))
for(r in 2:length(delta_list)){
    cat("Doing delta =",delta_list[[r]],"-> ")
    out_list[[r]] = mclapply(1:nrep, function(s) get_pvals(s, delta=delta_list[[r]]), mc.cores=16)
    out_list[[r]] = matrix(unlist(out_list[[r]]), ncol=1, byrow=T)
    cat("done\n")
}
out_list[[1]] = out_list[[2]]

# calculate power curve
z1 = lapply(out_list, function(y) apply(y,2,function(x) mean(x<=.05)))
z1 = matrix(unlist(z1), ncol=1, byrow=T)
z1
