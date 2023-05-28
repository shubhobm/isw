# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("data.table","parallel","harmonicmeanp","circular")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

# ----------------------------------------------------------------------
# auxiliary functions
# ======================================================================

# calculate inverse CDF: use the actual data points instead of histogram
# histogram is represented as $idx --- case/person id, $x -- the locations
invCDF_data <- function(df, invCDFnsamples = 20){
    levs = seq(0,1,length.out=invCDFnsamples+1)
    levs = levs[levs<1]
    qq = df[, .(q = quantile(x, levs, type=1)), by = idx][order(idx,q)]
    M = matrix(qq$q, ncol=length(levs), byrow=TRUE)
    M
}

# basic data manipulation
circ2mat = function(circ) matrix(unlist(lapply(circ, as.numeric)), nrow=length(circ), byrow=T)
vec2hist = function(vec) hist(vec, breaks=break_vec)$counts
mat2dt <- function (x) {x = as.matrix(x); data.table(i = rep(1:nrow(x), ncol(x)), j =  rep(1:ncol(x), each = nrow(x)), x = as.vector(x))[, .(idx = i, x)]}

# alpha for heat kernel
alphafun = function(lambda, type=1, t=NULL){
    if(type==1){ # heat kernel
        ret = exp(-t*lambda)
    } else{ # diffusion kernel
        ret = ifelse(lambda==0, 0, 1/lambda^2)
    }
    ret
}

# function for pushforward
pushForward = function(X, T, type="circle", ell=1, ...){
    # push forward using eigenfunction
    lambda = (2*pi*ell/T)^2
    Phi1 = sqrt(2/T)*cos(2*pi*ell*X/T)
    Phi2 = sqrt(2/T)*sin(2*pi*ell*X/T)
    
    # inverse map
    cbind(invCDF_data(mat2dt(Phi1), ...), invCDF_data(mat2dt(Phi2), ...)) # don't scale
}

# ----------------------------------------------------------------------
# Simulation parameters and wrapper code
# ======================================================================
set.seed(10022020)
m = 100
t = (1:m)/(m+1)
n1 = 60
n2 = 40
p = 100
qseq = seq(0.01,0.99,(0.99-0.01)/50)

get_pvals = function(seed, delta){
    set.seed(1e3*seed)

    # generate data
    delta1 = delta*pi/180 # convert degree to radians
    mu1 = rnorm(n1, mean=0, sd=0.1)
    mu2 = rnorm(n2, mean=delta1, sd=0.1)
    X1_list = lapply(mu1, function(mu) rvonmises(p, circular(mu), 2))
    X2_list = lapply(mu2, function(mu) rvonmises(p, circular(mu), 2))

    # pushforward maps
    X1 = circ2mat(X1_list)
    X2 = circ2mat(X2_list)
    E1 = Reduce(cbind, lapply(1:10, function(i) pushForward(X1, T=2*pi, ell=i)))
    E2 = Reduce(cbind, lapply(1:10, function(i) pushForward(X2, T=2*pi, ell=i)))

    #apply t-test to each dimension and combine
    p_ttests = sapply(seq_len(ncol(E1)), function(i) try(t.test(E1[, i], E2[, i])$p.value))
    tol = 1e-5
    p_ttests[which(is.na(p_ttests))] = 1-tol # set NA p-values to close to 1
    p_ttests[which(p_ttests >= 1-tol)] = 1-tol # set large to close to 1
    p_ttests[which(p_ttests <= tol)] = tol # set very small p-values tp small positive value

    # overall p-value using harmonic mean combo
    p_harmonic = p.hmp(p_ttests, L=length(p_ttests), multilevel=FALSE)

    # #overall p-value using Cauchy combo
    p_ttests1 = p_ttests
    p_cauchy = 0.5 - atan(mean(tan((0.5-p_ttests1)*pi)))/pi
    c(p_harmonic, p_cauchy)
}

# ----------------------------------------------------------------------
# Run the experiments
# ======================================================================
delta_list = seq(0,30,by=.25) # 0-8 degrees, will convert to radians inside
OUTPUT_FILE = '/path/to/output/file.rds'

out_list = vector("list",length(delta_list))
for(r in 1:length(delta_list)){
    cat("Doing delta =",delta_list[[r]],"-> ")
    out_list[[r]] = matrix(unlist(
        mclapply(
            1:1e3, 
            function(s) get_pvals(s, delta=delta_list[[r]]), 
            mc.cores=80)
        ), ncol=5, byrow=T)
    cat("done\n")
}
saveRDS(out_list, OUTPUT_FILE)


# calculate power curve
z1 = lapply(out_list, function(y) apply(y,2,function(x) mean(x<=.05)))
z1 = matrix(unlist(z1), ncol=5, byrow=T)
z1
