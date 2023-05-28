# ----------------------------------------------------------------------
# load libraries
# ======================================================================
Required.Packages <- c("data.table","parallel","harmonicmeanp","circular","squash")
sapply(Required.Packages, FUN = function(x) {suppressMessages(require(x, character.only = TRUE))})

# ----------------------------------------------------------------------
# auxiliary functions
# ======================================================================
circ2mat = function(circ) matrix(unlist(lapply(circ, as.numeric)), nrow=length(circ), byrow=T)
getMid = function(vec) vec[-length(vec)] + diff(vec) / 2
vec2hist2dt = function(vec, ...){
    n = length(vec)/2
    h = hist2(vec[1:n], vec[n+(1:n)], plot=F, ...)
    dt = data.table(wt=as.numeric(h$z), expand.grid(getMid(h$x),getMid(h$y)))
    names(dt)[-1] = paste0("x",1:2)
    dt
} 

# calculate inverse CDF using bivariate histogram bins
#histogram is represented as $idx --- case/person id,
#($x1,$x2) -- bivariate bin center locations, $wt the bin counts/weights
invCDF_data_wtd2 <- function(df, invCDFnsamples=20){
    #levs = seq(0.05, 0.95, length.out=invCDFnsamples) #want to capture 100th percentile too
    levs = seq(0, 1, length.out=invCDFnsamples+2)
    levs = levs[levs>0 & levs<1]    
    #just in case: aggregate weights by location
    df = df[,.(wt=sum(wt)), keyby=.(idx, x1, x2)]
    sumfun = function(w,xx1,xx2){ # bivariate cumsum
        sapply(1:length(w),function(i) sum(w*(xx1<=xx1[i] & xx2<=xx2[i])))
    }
    df = df[,lev:=sumfun(wt,x1,x2)/sum(wt), by=idx][order(idx, lev)]
    
    q = CJ(idx=unique(df$idx), lev=levs)
    setkey(df, "idx", "lev")
    setkey(q, "idx", "lev")
    qq = df[q,,roll=-Inf][order(idx, lev)]
    # one entry per (idx,lev): last one
    qq = qq[,.(
        x1 = rev(x1)[1],
        x2 = rev(x2)[1]), by=.(idx,lev)]
    
    Mx1 = matrix(qq$x1, ncol=length(levs), byrow=TRUE)
    Mx2 = matrix(qq$x2, ncol=length(levs), byrow=TRUE)
    cbind(Mx1,Mx2)
}

#two-dimensional pushforward: input is in the dataframe form
pushForward2 = function(df, invCDFnsamples=20, type="circle", T1, T2=NULL, ell1=1, ell2=1){
    
    # push forward using eigenfunction
    if(type=="circle"){
        T = T1
        ell = ell1
        lambda = (2*pi*ell/T)^2
        
        Phi1 = copy(df)[, x:=sqrt(2/T)*cos(2*pi*ell*x/T)]
        Phi2 = copy(df)[, x:=sqrt(2/T)*sin(2*pi*ell*x/T)]
        
        # inverse map
        E = cbind(invCDF_data_wtd(Phi1, invCDFnsamples), invCDF_data_wtd(Phi2, invCDFnsamples))
    } else if(type=="cylinder"){
        lambda = (2*pi*ell1/T1)^2 + (pi*ell2/T2)^2
        x = as.matrix(df[,c("x1","x2")])
        Phi1 = sqrt(2/T1)*cbind(cos(2*pi*ell1*x/T1),sin(2*pi*ell1*x/T1))
        Phi2 = sqrt(2/T2)*cos(pi*ell2*x/T2)
        Phi_cos = data.table(idx=df[,idx],wt=df[,wt],Phi1[,1:2]*Phi2)
        Phi_sin = data.table(idx=df[,idx],wt=df[,wt],Phi1[,3:4]*Phi2)
        
        E = cbind(invCDF_data_wtd2(Phi_cos, invCDFnsamples),
                  invCDF_data_wtd2(Phi_sin, invCDFnsamples))
    } else{
        T = T1
        ell = ell1
        lambda = (pi*ell/T)^2
        Phi = copy(df)[, x:= sqrt(2/T)*cos(pi*ell*x/T)]
        
        # inverse map
        E = invCDF_data_wtd(Phi, invCDFnsamples)
    }
    E
}

# get histogram wt for cylindrical distribution-- mardia-sutton 1978 distribution on cylinder
# notation as per "Extensions of probability distributions on torus, cylinder and disc", pg35
rnormt = function(n, mu=0, sigma=1, lb=0, ub=3){
    out = rnorm(n,mu,sigma)
    out[which(out > ub)] = ub
    out[which(out < lb)] = lb
    out
}

pcyl = function(theta, x, mu, mu0, kappa, rho1, rho2, sigma){
    muc = mu + sqrt(kappa)*sigma*(rho1*(cos(theta) - cos(mu0)) + rho2*(sin(theta) - sin(mu0)))
    rho = sqrt(rho1^2+rho2^2)
    sigmac = sigma^2*(1-rho^2)
    pvonmises(circular(theta), circular(mu0), kappa)*pnorm(x, muc, sigmac)
}

rcyl = function(n, mu, mu0, kappa, rho1, rho2, sigma){
    theta = rvonmises(n, circular(mu0), kappa)
    muc = mu + sqrt(kappa)*sigma*(rho1*(cos(theta) - cos(mu0)) + rho2*(sin(theta) - sin(mu0)))
    rho = sqrt(rho1^2+rho2^2)
    sigmac = sigma^2*(1-rho^2)
    x = rnormt(n, muc, sigmac, lb=0, ub=2*pi)
    c(as.numeric(theta), x)
}

# ----------------------------------------------------------------------
# Simulation parameters and wrapper code
# ======================================================================
set.seed(10022020)
# 20 bins per circle, 5 such circles
m1 = 20
m2 = 5
n = 500
xmid = 2*pi*(1:m1)/(m1+1)
ymid = pi*(1:m2)/(m2+1)
locs = as.matrix(expand.grid(x=xmid,y=ymid))
n1 = 60
n2 = 40
rho1 = .5
rho2 = .5
sigma = 1
kappa = 2

get_pvals = function(seed, delta, L, D){
    set.seed(1e3*seed)
    
    # generate data
    mu = runif(n1+n2,0,.5)
    mu0 = runif(n1+n2,0,.5)
    delta1 = delta*pi/180 # convert degree to radians
    X1 = do.call(rbind, lapply(1:n1, function(i)
        rcyl(n,mu[i],mu0[i], kappa, rho1, rho2, sigma)))
    X2 = do.call(rbind, lapply(1:n2, function(i)
        rcyl(n,mu[n1+i]+delta1,mu0[n1+i]+delta1, kappa, rho1, rho2, sigma)))
    
    # generate bivariate histograms, then data table of index, weights and locations
    df1 = rbindlist(lapply(1:n1, function(i) data.table(idx=i,vec2hist2dt(X1[i,], nx=20, ny=10))))
    df1[is.na(wt),wt := 0]
    df2 = rbindlist(lapply(1:n2, function(i) data.table(idx=i,vec2hist2dt(X2[i,], nx=20, ny=10))))
    df2[is.na(wt),wt := 0]
    
    # pushforward maps
    E1 = Reduce(cbind, lapply(1:L, function(i)
        Reduce(cbind, lapply(1:L, function(j)
            pushForward2(df1, type="cylinder", T1=2*pi, T2=2*pi, ell1=i, ell2=j, invCDFnsamples=D)))))
    E2 = Reduce(cbind, lapply(1:L, function(i)
        Reduce(cbind, lapply(1:L, function(j) 
            pushForward2(df2, type="cylinder", T1=2*pi, T2=2*pi, ell1=i, ell2=j, invCDFnsamples=D)))))
    
    #apply t-test to each dimension and combine
    p_ttests = sapply(seq_len(ncol(E1)), function(i)
        tryCatch({t.test(E1[, i], E2[, i])$p.value}, error=function(e){1.0}))
    #overall p-value using harmonic mean combo
    p_harmonic = p.hmp(p_ttests, L=length(p_ttests), multilevel=FALSE)
    p_harmonic
}

# ----------------------------------------------------------------------
# Run the experiments
# ======================================================================
delta_list = seq(0,20,by=.1) # 0-8 degrees, will convert to radians inside
LD_grid = expand.grid(L=4, D=24)
OUTPUT_FILE = '/path/to/output/file.rds'

# this runs for the full grid. consider taking out a single iteration and run for test purposes
all_out = vector("list",nrow(LD_grid))
for(i in 1:nrow(LD_grid)){
    L = LD_grid[i,1]
    D = LD_grid[i,2]
    cat("=====\nDoing L =",L,", D=",D,"\n=====\n")
    all_out[[i]] = lapply(delta_list, function(delta){
        cat("Doing delta =",delta,"-> ")
        out = mclapply(1:1e3, function(s) get_pvals(s, delta=delta, L=L, D=D), mc.cores=16)
        cat("done\n")
        out
    })
    saveRDS(list(LD_grid,all_out), file=OUTPUT_FILE) # save per iteration for any runtime errors
}

saveRDS(list(LD_grid,all_out), file=OUTPUT_FILE)

# calculate power curve for one value of (L,D')
z1 = lapply(all_out[[1]], function(y) apply(y,2,function(x) mean(x<=.05)))
z1 = matrix(unlist(z1), ncol=1, byrow=T)
z1