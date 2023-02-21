## Here are the functions to call within R to run K-Mix (kernel mixture GMMs)
##   or the Bayesian implementation with diffuse prior "bmix"
# Before restarting efforts on 1/7/18, this file was the latest version from 6/11/17
# v34 and v36 had earlier dates (last updated April 2017)
#

## datagen_v3.fun <- function(Xdim, caseno, seedno){
## sq.err <- function(x,y=test[,1],printout=TRUE){
## phi2p.fun <- function(phi){
## indivpred.fun <- function(par, mat, d, xind, yind){
## iterpred.fun <- function(par2, mat, d, xind, yind){
## getpred.fun <- function(parmat, mat, d=NULL, n=n, yind=NULL){
## get.sse <- function(par, p0){
## find.alpha <- function(prob){
## trymclust <- function(Xmat,kclust){
## run.mclust <- function(X, sampsize=NULL, seedno=NULL, cmax=NULL, printlevel=2){
## KMix <- function(X, S, wl, n=NULL, Burn=1000, J=100, Thin=10, rho=NULL, a=1, b=1, seed="patrick",
##                  iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind, out.type=4){
## BMix <- function(X, S, wl, n=NULL, Burn=1000, J=100, Thin=10, rho=NULL, a=1, b=1, seed="patrick",
##                  iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind=1, out.type=4){
## KMix.mod and BMix.mod are not more efficient, restarts are unnecessary for the Gibbs sampler.

## 
## FUNCTIONS FROM V1:
##
sq.err <- function(x,y=test[,1],printout=TRUE){
    if (printout) {
        cat("Error is: ",sum((x-y)^2)/sum((mean(y)-y)^2),"\n")
    } else {
        return(sum((x-y)^2)/sum((mean(y)-y)^2))
    }
}
phi2p.fun <- function(phi){
    ## converts phi variables to mixture proportions p
    phi <- cbind(phi,1)
    k <- ncol(phi)
    phi2 <- 1 - phi
    for (j in 2:k){
        phi2[,j] <- phi2[,j-1] * phi2[,j]
        phi[,j] <- phi2[,j-1] * phi[,j]
    }	
    return(phi)
}
indivlik.fun <- function(par, mat, d, xind, yind){
    ## likelihood calculation
    return(par[1] * mvtnorm::dmvnorm(mat, mean=par[1+xind], 
                            sigma=as.matrix(matrix(par[-(1:(d+1))],d,d)[-yind,-yind])))
}
indivpred.fun <- function(par, mat, d, xind, yind){
    ## predicted value calculation
    sig <- matrix(par[-(1:(d+1))],d,d)
    return(sig[yind,-yind] %*% solve(sig[-yind,-yind]) %*% (t(mat) - par[1+xind]) + par[1+yind])
}
iterpred.fun <- function(par2, mat, d, xind, yind){
    ## get prediction per iteration
    lik <- apply(par2, 2, indivlik.fun, mat=as.matrix(mat), d=d, xind, yind)
    lik[,1] <- lik[,1] + 1e-300
    pred <- apply(par2, 2, indivpred.fun, mat=mat, d=d, xind, yind)
    rowSums( pred * lik * (1 / rowSums(lik)) )
}
getpred.fun <- function(parmat, mat, d=NULL, n=n, yind=NULL){
    ## get overall prediction: the mean of all iteration predicted values
    if (is.null(d)) {d <- ncol(mat)}
    if (is.null(yind)) {yind <- d}
    xind <- setdiff(1:d,yind)
    apply(parmat[,1:n,], 3, iterpred.fun, mat=mat[,xind], d=d, xind=xind, yind=yind)
}
get.sse <- function(par, p0){
    p0 <- p0[order(p0, decreasing=TRUE)]
    if (length(p0) < 100) {
        p0 <- c(p0, rep(0,100-length(p0)))
    } else p0 <- p0[1:100]
    sum((p0 - par^(0:99)/(par+1)^(1:100))^2)
}
find.alpha <- function(probmat, ranklist, G=3){
    ranklist <- setdiff(ranklist, 1)
    alphaval <- rep(NA, G)
    for (i in 1:G){
        k <- ranklist[i]
        alphaval[i] <- optimize(get.sse, c(0,5), p0=probmat[k,1:k])$min
        if (alphaval[i] > 4)
            alphaval[i] <- optimize(get.sse, c(0,50), p0=probmat[k,1:k])$min
    }
    ## return(weighted.mean(alphaval,w=G:1))
    return(mean(alphaval))
}
