## Here are the functions to call within R to run K-Mix (kernel mixture GMMs)
##   or the Bayesian implementation with diffuse prior "bmix"
# Before restarting efforts on 1/7/18, this file was the latest version from 6/11/17
# v34 and v36 had earlier dates (last updated April 2017)
#
library(mclust)
#library(MixSim)
library(MASS)
library(mvtnorm)
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

datagen_v3.fun <- function(Xdim, caseno, seedno){
    ## generate test cases for Xdim = 8, 16 predictors
    ## Version 2 of Test Case generation, based on Maitra's comments
    ## 1-2 (r), 3-4 (u), 5-6 (s), 7-8 (h), 9-10 (f)
    ## Scenarios:
    ##  r - MVN generation of 5-density dsn
    ##  u - mixture of increasing variance to simulate M.V. t-dsn
    ##  s - presence of spurious (noisy) predictors
    ##  h - high component overlap case (avg .05/max .10 instead of avg .02/max .04)
    ##  f - few component dsn (G=4 instead G=8)
    ## seedno 1-25
    
    G <- 8  # an 8-component mixture dsn
    ## my thinking is that if this works on a 5-component mixture, the MVN assumption is more acceptable.
    unif <- spur <- FALSE
    ## MixSim parameters:
    ## BarOmega <- 0.05  # this is average component pairwise overlap
    ## MaxOmega <- 0.15  # this is maximum allowed comp. pairwise overlap
    BarOmega <- 0.02  # this is average component pairwise overlap
    MaxOmega <- 0.06  # this is maximum allowed comp. pairwise overlap

    if (caseno %in% 3:4) unif <- TRUE
    if (caseno %in% 5:6) spur <- TRUE
    if (caseno %in% 7:8){
        BarOmega <- 0.05  # this is average component separation
        MaxOmega <- 0.15  # this is maximum allowed comp. separation
    }
    if (caseno %in% 9:10) G <- 4

    noisevar <- 0
    p <- Xdim + 1
    if (spur){
        noisevar <- round(Xdim * .25)
        p <- p - noisevar
    }

    ## sample size from Banks (2003). k=4 is his small set, to
    ##   show how well the method performs under sparse conditions
    ## k=10 is the medium-sized case, for a plentiful data comparison
    ## (it's Bayesian, so should focus on smallish data conditions)
    ## And test points always 10000 in Banks (2003) no matter the dimension.
    #klist <- c(.04, .02, .01)
    klist <- c(0.5,2)*(1+0.5*(Xdim+2)*(Xdim+1)+(Xdim+1))*8  # Even when G=4, use 8
    k <- klist[(caseno-1) %% 2 + 1]
    ## Ntr <- k*ceiling((40*Xdim)^1.1) # p / (n^.9) = .01,.02,.04
          ## changed from 2^Xdim * k
    Ntr <- k
    Nte <- 5000  # reduced from 10000 for testing speed,
                                        # and considering size of training sets

    ## have larger data sets be an extension of the small sets, same MixSim arrangement
    set.seed(25*floor((caseno-1)/4)+seedno)   # seeds 1 to 25 for scenario 1:4, ...,
                                        #  26 to 50 for scenario 5:6 and 7:8 (fund. different)
                                        #  and 51 to 75 for scenario 9:10 
    a <- MixSim(BarOmega=BarOmega, MaxOmega=MaxOmega, K=G, p=p,
                PiLow=1/(2*G), int=c(-1,1), resN=5000)

    set.seed(75+seedno) # reuse seeds 76 to 100 to have similar data per seed across scenarios
    b <- simdataset(Ntr,a$Pi,a$Mu,a$S, n.noise=noisevar, n.out=0)
    ## for poor results of rgamma(n,2,1), see run_mpdpm_v34.Rout, scenarios 3 and 4.
    delta <- 1/rgamma(nrow(b$X),4,4)
    if (unif){
        ## draw from a multivariate t-8 dsn (see Gamerman, pp.24-25) 
        delvals <- tabulate(ceiling(delta*20))
        for (i in 1:length(delvals)){
            if (length(delvals[i]) == 0){
                next
            } else {
                nsim = max(delvals[i], 40*G)
                b1 <- simdataset(nsim,a$Pi,a$Mu,((i-.5)/20)*a$S, n.noise=noisevar, n.out=0)
                drawrows <- sample(1:nsim, delvals[i])
                b$X <- rbind(b$X, b1$X[drawrows,])
                b$id <- c(b$id, b1$id[drawrows])
                ## cat(i/10,delvals[i],':')
            }
        }
        b$X <- b$X[-(1:Ntr),]
        b$id <- b$id[-(1:Ntr)]
    }
    ## see bottom of mvt_testing.R for confirmation of functionality:
    ## 1e6 POINTS USING PACKAGE MVTNORM FUNCTION rmvt
    ## > sapply(1:10*50, function(i) sum(s2 > i)/length(s2))
    ## [1] 0.053411 0.006331 0.001553 0.000582 0.000265 0.000134 0.000075 0.000050 0.000034 0.000023
    ## 1e5 POINTS USING ABOVE METHOD
    ## > sapply(1:10*50, function(i) sum(sX > i)/length(sX))
    ## [1] 0.05335 0.00617 0.00162 0.00066 0.00033 0.00016 0.00004 0.00002 0.00000 0.00000

    set.seed(100+seedno) # scenarios 1-4,7-8 will have same test data.
                                        # (5-6 different due to noisevar, 9-10 due to G=10)
    b2 <- simdataset(Nte,a$Pi,a$Mu,a$S, n.noise=noisevar, n.out=0)

    ## shuffle of order or rows and columns of training data:
    set.seed(125+seedno)  # this won't matter, but is a good idea
    cid <- c(1,sample(2:ncol(b$X)))
    xid <- sample(1:nrow(b$X))
    tr <- b$X[xid,cid]
    mem <- list()
    mem$tr <- b$id[xid]

    ## first calculate E(Y) given membership for test set:
    ## remove the spurious predictors not in gmdist 'a' (placed at last columns)
    cr <- 2:p
    trueval <- rep(NA, Nte)
    for (groupno in unique(b2$id)){
        trueval[b2$id == groupno] <- a$Mu[groupno,1] +
            t(t(b2$X[b2$id == groupno,cr]) - a$Mu[groupno,cr]) %*%
               (solve(a$S[cr,cr,groupno]) %*% a$S[1,cr,groupno])
    }
    ## Adding inv(SigmaX) is correct above as we have MVN densities.
    ## for (groupno in unique(b2$id)){
    ##     trueval[b2$id == groupno] <- a$Mu[groupno,1] +
    ##         t(t(b2$X[b2$id == groupno,cr]) - a$Mu[groupno,cr]) %*% a$S[1,cr,groupno]
    ## } 

    ## then transform columns as done for training set:
    te <- b2$X[,cid]
    mem$te <- b2$id

##     ## can't explain this:
## plot(te[,1],col=2)
## points(b2$id/4, pch=16, cex=.2)
## points(b2$X[,1], col=3, cex=.5)
## points(a$Mu[b2$id,1],col=4,cex=.3)
## points(trueval, col=5, cex=.25)

    return(list(tr=tr, te=te, mem=mem, orac=trueval))
}
##
## NOTE: datagen_v37.fun is NOT USED for Paper 1 test cases
##
datagen_v37.fun <- function(Xdim, caseno, seedno){
    ## generate test cases for Xdim = 8, 16 predictors
    ## Version 2 of Test Case generation, based on Maitra's comments
    ## 1-2 (r), 3-4 (u), 5-6 (s), 7-8 (h), 9-10 (f)
    ## Scenarios:
    ##  r - MVN generation of 5-density dsn
    ##  u - mixture of increasing variance to simulate M.V. t-dsn
    ##  s - presence of spurious (noisy) predictors
    ##  h - high component overlap case (avg .05/max .10 instead of avg .02/max .04)
    ##  f - few component dsn (G=4 instead G=8)
    ## seedno 1-25
    
    G <- 8  # an 8-component mixture dsn
    ## my thinking is that if this works on a 5-component mixture, the MVN assumption is more acceptable.
    unif <- spur <- FALSE
    ## MixSim parameters:
    ## BarOmega <- 0.05  # this is average component pairwise overlap
    ## MaxOmega <- 0.15  # this is maximum allowed comp. pairwise overlap
    BarOmega <- 0.075  # this is average component pairwise overlap
    MaxOmega <- 0.15 # this is maximum allowed comp. pairwise overlap

    if (caseno %in% 3:4) unif <- TRUE
    if (caseno %in% 5:6) spur <- TRUE
    if (caseno %in% 7:8){
        BarOmega <- 0.05  # this is average component separation
        MaxOmega <- 0.15  # this is maximum allowed comp. separation
    }
    if (caseno %in% 9:10) G <- 4

    noisevar <- 0
    p <- Xdim + 1
    if (spur){
        noisevar <- round(Xdim * .25)
        p <- p - noisevar
    }

    ## sample size from Banks (2003). k=4 is his small set, to
    ##   show how well the method performs under sparse conditions
    ## k=10 is the medium-sized case, for a plentiful data comparison
    ## (it's Bayesian, so should focus on smallish data conditions)
    ## And test points always 10000 in Banks (2003) no matter the dimension.
    #klist <- c(.04, .02, .01)
    klist <- c(0.5,2)*(1+0.5*(Xdim+2)*(Xdim+1)+(Xdim+1))*8  # Even when G=4, use 8
    k <- klist[(caseno-1) %% 2 + 1]
    ## Ntr <- k*ceiling((40*Xdim)^1.1) # p / (n^.9) = .01,.02,.04
          ## changed from 2^Xdim * k
    Ntr <- k
    Nte <- 5000  # reduced from 10000 for testing speed,
                                        # and considering size of training sets

    ## have larger data sets be an extension of the small sets, same MixSim arrangement
    set.seed(25*floor((caseno-1)/4)+seedno)   # seeds 1 to 25 for scenario 1:4, ...,
                                        #  26 to 50 for scenario 5:6 and 7:8 (fund. different)
                                        #  and 51 to 75 for scenario 9:10 
    a <- MixSim(BarOmega=BarOmega, MaxOmega=MaxOmega, K=G, p=p,
                PiLow=1/(2*G), int=c(-1,1), resN=5000)

    set.seed(75+seedno) # reuse seeds 76 to 100 to have similar data per seed across scenarios
    b <- simdataset(Ntr,a$Pi,a$Mu,a$S, n.noise=noisevar, n.out=0)
    delta <- 1/rgamma(nrow(b$X),2,1)
    if (unif){
        ## draw from a multivariate t-2 dsn 
        delvals <- tabulate(ceiling(delta*20))
        for (i in 1:length(delvals)){
            if (length(delvals[i]) == 0){
                next
            } else {
                nsim = max(delvals[i], 40*G)
                b1 <- simdataset(nsim,a$Pi,a$Mu,((i-.5)/20)*a$S, n.noise=noisevar, n.out=0)
                drawrows <- sample(1:nsim, delvals[i])
                b$X <- rbind(b$X, b1$X[drawrows,])
                b$id <- c(b$id, b1$id[drawrows])
                ## cat(i/10,delvals[i],':')
            }
        }
        b$X <- b$X[-(1:Ntr),]
        b$id <- b$id[-(1:Ntr)]
    }
    set.seed(100+seedno) # scenarios 1-4,7-8 will have same test data.
                                        # (5-6 different due to noisevar, 9-10 due to G=10)
    b2 <- simdataset(Nte,a$Pi,a$Mu,a$S, n.noise=noisevar, n.out=0)

    ## shuffle of order or rows and columns of training data:
    set.seed(125+seedno)  # this won't matter, but is a good idea
    cid <- c(1,sample(2:ncol(b$X)))
    xid <- sample(1:nrow(b$X))
    tr <- b$X[xid,cid]
    mem <- list()
    mem$tr <- b$id[xid]

    ## first calculate E(Y) given membership for test set:
    ## remove the spurious predictors not in gmdist 'a' (placed at last columns)
    cr <- 2:p
    trueval <- rep(NA, Nte)
    for (groupno in unique(b2$id)){
        trueval[b2$id == groupno] <- a$Mu[groupno,1] +
            t(t(b2$X[b2$id == groupno,cr]) - a$Mu[groupno,cr]) %*%
               (solve(a$S[cr,cr,groupno]) %*% a$S[1,cr,groupno])
    }
    ## Adding inv(SigmaX) is correct above as we have MVN densities.
    ## for (groupno in unique(b2$id)){
    ##     trueval[b2$id == groupno] <- a$Mu[groupno,1] +
    ##         t(t(b2$X[b2$id == groupno,cr]) - a$Mu[groupno,cr]) %*% a$S[1,cr,groupno]
    ## } 

    ## then transform columns as done for training set:
    te <- b2$X[,cid]
    mem$te <- b2$id

##     ## can't explain this:
## plot(te[,1],col=2)
## points(b2$id/4, pch=16, cex=.2)
## points(b2$X[,1], col=3, cex=.5)
## points(a$Mu[b2$id,1],col=4,cex=.3)
## points(trueval, col=5, cex=.25)

    return(list(tr=tr, te=te, mem=mem, orac=trueval))
}