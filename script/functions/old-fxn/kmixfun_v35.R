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
## a <- MixSim(BarOmega=BarOmega, MaxOmega=MaxOmega, K=3, p=2, PiLow=1/6, int=c(-1,1), resN=250)
## b <- simdataset(100,a$Pi,a$Mu,a$S, n.noise=0, n.out=20, alpha=.001, max.out=1e5)
## plot(b$X[,1], b$X[,2], col=b$id+1)

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
trymclust <- function(Xmat,kclust,d){
    tryCatch({
        if (d==1)
            out <- Mclust(Xmat, G=kclust)
        else
            out <- Mclust(Xmat, G=kclust, modelNames=c("VVV", "EEV", "VEV", "VVI", "EVI", "VEI"))
        return(out)
    }, error=function(cond) {
        message(cond,"\n")
        return(0)
    })
}
trymclust_noI <- function(Xmat,kclust,d){
  tryCatch({
    if (d==1)
      out <- Mclust(Xmat, G=kclust)
    else
      out <- Mclust(Xmat, G=kclust, modelNames=c("VVV", "EEV", "VEV"))
    return(out)
  }, error=function(cond) {
    message(cond,"\n")
    return(0)
  })
}
trymclust_vvv <- function(Xmat,kclust,d){
  tryCatch({
    if (d==1)
      out <- Mclust(Xmat, G=kclust)
    else
      out <- Mclust(Xmat, G=kclust, modelNames="VVV")
    return(out)
  }, error=function(cond) {
    message(cond,"\n")
    return(0)
  })
}
run.mclust_old <- function(X, sampsize=NULL, seedno=NULL, cmax=NULL, printlevel=2){
    ## NOTE: I DON'T THINK THIS FUNCTION IS USED ANYMORE. SEE run.mclust_v2() THAT FOLLOWS:
    ## output 'prob' = mixture proportions for best BIC, used to find optimal a/b ratio
    d <- ncol(X)
    if (is.null(cmax)) cmax <- 3*d
    if (!is.null(seedno)) set.seed(seedno)
    if (is.null(sampsize)) sampsize <- nrow(X)
    if (sampsize > nrow(X)) sampsize <- nrow(X)
    idx <- sample(1:nrow(X), sampsize)
    idx <- idx[order(idx)]
    Xsub <- X[idx,]
    priorwt <- NULL
    priormat <- array(NA, dim=c(d,d,sum(1:cmax)))
    biclist <- rep(NA,cmax)
    biclimit <- bestbic <- -Inf
    bicstop <- FALSE
    prokeep <- matrix(NA, cmax, cmax)
    k <- 0
    if (printlevel > 1) cat(cmax,":: ",sep="")
    while (k < cmax & !bicstop){
        k <- k + 1
        out <- trymclust(Xmat=Xsub, kclust=k, d=d)
        if (is.list(out)){
            if (k==1) biclimit <- as.double(out$bic)
            if (out$bic < biclimit & k >= 12) bicstop <- TRUE
            biclist[k] <- out$bic
            prokeep[k,1:k] <- out$param$pro
            if (out$bic >= bestbic){
                bestbic <- out$bic
            }
            if (printlevel > 1) cat(k,out$modelName,round(out$bic - biclimit,1),"/ ")
            priorwt <- c(priorwt, out$param$pro)
            priormat[,,1:k+sum(1:k)-k] <- out$param$variance$sigma[,,]
        } else {
            bicstop <- TRUE
            k <- k - 1
        }
    }
    if (printlevel > 1) cat("\n")
    priormat <- priormat[,,1:sum(1:k)]
    biclist <- biclist[1:k]
    bicrank <- rank(biclist)
    for (i in 1:k){
      	isrange <- 1:i+sum(1:i)-i
	      priorwt[isrange] <- priorwt[isrange] * bicrank[i] 
    }

    ranklist <- (1:length(biclist))[order(biclist,decreasing=T)]
    rankn <- min(length(ranklist),8)
    if (printlevel > 0){
        cat("For seedno", seedno, "mclust finished (", length(biclist),
            ")  Best bic groups are: ", ranklist[1:rankn], "\n")
    }
    return(list(mat = priormat, wt = priorwt, bic = biclist, prob=prokeep, rl=ranklist))
}
run.mclust_v2 <- function(X, sampsize=NULL, seedno=NULL, cmax=NULL, cmin=12, printlevel=2){
    ## output 'prob' = mixture proportions for best BIC, used to find optimal a/b ratio
    d <- ncol(X)
    if (is.null(cmax)) cmax <- 3*d
    if (!is.null(seedno)) set.seed(seedno)
    if (is.null(sampsize)) sampsize <- nrow(X)
    if (sampsize > nrow(X)) sampsize <- nrow(X)
    idx <- sample(1:nrow(X), sampsize)
    idx <- idx[order(idx)]
    Xsub <- X[idx,]
    priorwt <- NULL
    priormat <- array(NA, dim=c(d,d,sum(1:cmax)))
    biclist <- rep(NA,cmax)
	  logliklist <- rep(NA,cmax)
    biclimit <- bestbic <- -Inf
    bicstop <- FALSE
    prokeep <- matrix(NA, cmax, cmax)
    k <- 0
    if (printlevel > 1) cat(cmax,":: ",sep="")
    while (k < cmax & !bicstop){
        k <- k + 1
        idx <- sample(1:nrow(X), sampsize)
        idx <- idx[order(idx)]
        Xsub <- X[idx,]

        out <- trymclust(Xmat=Xsub, kclust=k, d=d)
        if (is.list(out)){
            if (k==1) biclimit <- as.double(out$bic)
            if (out$bic < biclimit & k >= cmin) bicstop <- TRUE
            biclist[k] <- out$bic
		      	logliklist[k] <- out$loglik
            prokeep[k,1:k] <- out$param$pro
            if (out$bic >= bestbic){
                bestbic <- out$bic
            }
            if (printlevel > 1) cat(k,out$modelName,round(out$bic - biclimit,1),"/ ")
            priorwt <- c(priorwt, rep(k^(-0.5),k))
            ## priorwt <- c(priorwt, rep(1,k))
            priormat[,,1:k+sum(1:k)-k] <- out$param$variance$sigma[,,]
        } else {
            bicstop <- TRUE
            k <- k - 1
        }
    }
    if (printlevel > 1) cat("\n")
    priormat <- priormat[,,1:sum(1:k)]
    biclist <- biclist[1:k]
    bicrank <- rank(biclist)
    ## for (i in 1:k){
    ##     isrange <- 1:i+sum(1:i)-i
    ##     priorwt[isrange] <- priorwt[isrange] * bicrank[i] 
    ## }

    ranklist <- (1:length(biclist))[order(biclist,decreasing=T)]
    rankn <- min(length(ranklist),8)
    if (printlevel > 0){
        cat("For seedno", seedno, "mclust finished (", length(biclist),
            ")  Best bic groups are: ", ranklist[1:rankn], "\n")
    }
    return(list(mat = priormat, wt = priorwt, bic = biclist, loglik = logliklist, n = sampsize, prob=prokeep, rl=ranklist))
}
run.mclust_v2B <- function(X, sampsize=NULL, seedno=NULL, cmax=NULL, cmin=12, printlevel=2){
    ## output 'prob' = mixture proportions for best BIC, used to find optimal a/b ratio
    d <- ncol(X)
    if (is.null(cmax)) cmax <- 3*d
    if (!is.null(seedno)) set.seed(seedno)
    if (is.null(sampsize)) sampsize <- nrow(X)
    if (sampsize > nrow(X)) sampsize <- nrow(X)
    idx <- sample(1:nrow(X), sampsize)
    idx <- idx[order(idx)]
    Xsub <- X[idx,]
    priorwt <- NULL
    priormat <- array(NA, dim=c(d,d,sum(1:cmax)))
    biclist <- rep(NA,cmax)
	  logliklist <- rep(NA,cmax)
    biclimit <- bestbic <- -Inf
    bicstop <- FALSE
    prokeep <- matrix(NA, cmax, cmax)
    k <- 0
    if (printlevel > 1) cat(cmax,":: ",sep="")
    while (k < cmax & !bicstop){
        k <- k + 1
        idx <- sample(1:nrow(X), sampsize)
        idx <- idx[order(idx)]
        Xsub <- X[idx,]

        out <- trymclust(Xmat=Xsub, kclust=k, d=d)
        if (is.list(out)){
            if (k==1) biclimit <- as.double(out$bic)
            if (out$bic < biclimit & k >= cmin) bicstop <- TRUE
            biclist[k] <- out$bic
	      		logliklist[k] <- out$loglik
            prokeep[k,1:k] <- out$param$pro
            if (out$bic >= bestbic){
                bestbic <- out$bic
            }
            if (printlevel > 1) cat(k,out$modelName,round(out$bic - biclimit,1),"/ ")
            ##priorwt <- c(priorwt, rep(k^(-0.5),k))
            ## priorwt <- c(priorwt, rep(1,k))
            priorwt <- c(priorwt, out$param$pro)
            priormat[,,1:k+sum(1:k)-k] <- out$param$variance$sigma[,,]
        } else {
            bicstop <- TRUE
            k <- k - 1
        }
    }
    if (printlevel > 1) cat("\n")
    priormat <- priormat[,,1:sum(1:k)]
    biclist <- biclist[1:k]
    bicrank <- rank(biclist)
    ## for (i in 1:k){
    ##     isrange <- 1:i+sum(1:i)-i
    ##     priorwt[isrange] <- priorwt[isrange] * bicrank[i] 
    ## }

    ranklist <- (1:length(biclist))[order(biclist,decreasing=T)]
    rankn <- min(length(ranklist),8)
    if (printlevel > 0){
        cat("For seedno", seedno, "mclust finished (", length(biclist),
            ")  Best bic groups are: ", ranklist[1:rankn], "\n")
    }
    return(list(mat = priormat, wt = priorwt, bic = biclist, loglik = logliklist, n = sampsize, prob=prokeep, rl=ranklist))
}
run.mclust_v3 <- function(X, sampsize=NULL, seedno=NULL, cmax=NULL, cmin=12, printlevel=2, noIflag=TRUE, VVVonly=FALSE){
    ## output 'prob' = mixture proportions for best BIC, used to find optimal a/b ratio
    ## obs = number of observations used for that clustering
    d <- ncol(X)
    if (is.null(cmax)) cmax <- 3*d
    if (!is.null(seedno)) set.seed(seedno)
    if (is.null(sampsize)) sampsize <- round(nrow(X)/4)
    if (sampsize > nrow(X)) sampsize <- nrow(X)
    # idx <- sample(1:nrow(X), sampsize)
    # idx <- idx[order(idx)]
    # Xsub <- X[idx,]
    priorwt <- NULL
    priormat <- array(NA, dim=c(d,d,sum(1:cmax)))
    biclist <- rep(NA,cmax)
    biclimit <- bestbic <- -Inf
    bicstop <- FALSE
    obs <- NULL
    prokeep <- matrix(NA, cmax, cmax)
    k <- 0
    if (printlevel > 1) cat(cmax,":: ",sep="")
    while (k < cmax & !bicstop){
        k <- k + 1
        idx <- sample(1:nrow(X), sampsize)
        idx <- idx[order(idx)]
        Xsub <- X[idx,]
		if (VVVonly){
			out <- trymclust_vvv(Xmat=Xsub, kclust=k, d=d)
		} else if (noIflag) {
			out <- trymclust_noI(Xmat=Xsub, kclust=k, d=d)
		} else {
			out <- trymclust(Xmat=Xsub, kclust=k, d=d)
		}
        if (is.list(out)){
            if (k==1) biclimit <- as.double(out$bic)
            if (out$bic < biclimit & k >= cmin) bicstop <- TRUE
            biclist[k] <- out$bic
            prokeep[k,1:k] <- out$param$pro
            if (out$bic >= bestbic){
                bestbic <- out$bic
            }
            if (printlevel > 1) cat(k,out$modelName,round(out$bic - biclimit,1),"/ ")
		      	pmult <- as.double(exp(0.5*(out$bic - biclimit)/nrow(X)))
            priorwt <- c(priorwt, rep(pmult*k^(-1),k))
            ## priorwt <- c(priorwt, rep(1,k))
            priormat[,,1:k+sum(1:k)-k] <- out$param$variance$sigma[,,]
			obsout <- nrow(Xsub) * out$parameters$pro
            obs <- c(obs, obsout)
        } else {
            bicstop <- TRUE
            k <- k - 1
        }
    }
    if (printlevel > 1) cat("\n")
    priormat <- priormat[,,1:sum(1:k)]
    biclist <- biclist[1:k]
    bicrank <- rank(biclist)
    ## for (i in 1:k){
    ##     isrange <- 1:i+sum(1:i)-i
    ##     priorwt[isrange] <- priorwt[isrange] * bicrank[i] 
    ## }

    ranklist <- (1:length(biclist))[order(biclist,decreasing=T)]
    rankn <- min(length(ranklist),8)
    if (printlevel > 0){
        cat("For seedno", seedno, "mclust finished (", length(biclist),
            ")  Best bic groups are: ", ranklist[1:rankn], "\n")
    }
    return(list(mat = priormat, wt = priorwt, bic = biclist, obs = obs, prob=prokeep, n = sampsize, rl=ranklist))
}
run.mclust_v3noI <- function(X, sampsize=NULL, seedno=NULL, cmax=NULL, cmin=12, printlevel=2, noIflag=NULL){
    error('Run run.mclust_v3 with flag noI=TRUE') 
}
MPDPM <- function(X, S, obs, wl, n=NULL, Burn=500, J=250, Thin=2, rho=NULL, a=1, b=1, seed="patrick",
                  iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind, out.type=4){
    ## Runs the K-Mix algorithm (in C code)
    ## X = N-by-d data matrix
    ## S = d-by-d-by-L prior covariance matrices
	## obs = observations used to create S
    ## wl = weights on prior covariance matrices
    ## Burn = burn-in iterations
    ## J = samples
    ## Thin = thinning iterations per sample
    ## rho = Wishart prior d.f. parameter (default=d+2) 
    ## a,b = beta prior parameters 
    ## seed = starting seed phrase
    ## iterprint = print progress statements if printlevel > 0
    ## printlevel = {0 NONE, 1, 2 DETAILED} for level of printout information 
    ## so.dir = character string for path of kmix.so file
    ## so.file = the kmix.so file
    ## Xtest = data matrix for the test set
    ## yind = index of response (to be predicted)
    ## out.type = output level:
    ##   1-K-Mix output parameters as list
    ##   2-K-Mix output parameters as array
    ##   3-prediction matrix for Y
    ##   4-prediction vector for Y (rowMeans of prediction matrix)
    
    ## error checking
    if (!(out.type %in% 1:4)) stop("Error: out.type must be in 1:4")
    if (is.null(n)) n <- ncol(X)*2
    if (is.null(rho)) rho <- ncol(X)+2
    if (is.null(so.dir)) so.dir <- paste(getwd(),"/code",sep="")
    if (so.dir == "pwd") so.dir <- getwd()
    if (substr(so.dir,nchar(so.dir),nchar(so.dir)) != "/")
        so.dir <-paste(so.dir, "/", sep="")
    if (is.null(so.file)) so.file <- "mpdpm.so"
    ## the difference between mpdpm.c and mpdpm2.c is multiplying the input scale
    ##   matrices by rho instead of rho+k+1, to put into form of Wishart
    ##   (rho*S is the mean and mode of the Wishart distribution, not (rho+k+1)*S...)
    if (is.null(dim(Xtest))) Xtest <- matrix(Xtest, ncol=1)
    
    ##-set up parameters for C program
    N <- nrow(X)
    k <- ncol(X)
    L <- dim(S)[3]
	
	## S smoothing
	Xm <- as.matrix(X)
	G <- (t(Xm) %*% Xm) * (nrow(X)^(-9/7))
	for (l in 1:L){
		if (obs[l] > 1){
			S[,,l] <- (1/obs[l]) * ( (obs[l]-1) * S[,,l] + G )
		} else {
		    S[,,l] <- G	
		}				
	}

    ## format covariance matrix and weight input for C program
    Smat <- matrix(NA, k*(k+1)/2+1, L)
    for (i in 1:k)
	for (j in 1:i)
            Smat[sum(1:i)-i+j,] <- S[i,j,]
    Smat[k*(k+1)/2+1,] <- wl

    ## load C program
    dyn.load(paste(so.dir, so.file, sep=""))

    ## call C program
    out <- .C("mpdpm",
              as.integer(c(N,k,n,rho,L,Burn,Thin,J,iterprint,printlevel)),
              as.double(c(a,b)),
              as.double(stack(as.data.frame(t(X)))$values),
              as.double(stack(as.data.frame(Smat))$values),
              ind = integer(N*J),
              phi = double((n-1)*J),
              mu = double(n*k*J),
              sigma = double(n*k*(k+1)/2*J),
              as.character(seed)
              )

    ## unload C program
    dyn.unload(paste(so.dir, so.file, sep=""))

    ## reformat C output for easier use (this will eventually be done in C itself)
    ## ind is fine; split par (parameter values) into phi, mu, sigma
    phimat <- matrix(out$phi,J,n-1,byrow=T)
    muvect <- matrix(out$mu,J,n*k,byrow=T)
    sigmawidth <- k*(k+1)/2 
    sigmavect <- matrix(out$sigma,J,n*sigmawidth,byrow=T)
    mumat <- array(NA, dim=c(J,k,n))
    sigmamat <- array(NA, dim=c(J,sigmawidth,n))
    for (i in 1:n){
        mumat[,,i] <- muvect[,1:k+(i-1)*k]
        sigmamat[,,i] <- sigmavect[,1:sigmawidth+(i-1)*sigmawidth]
    }

    ## return C output to R program
    if (out.type == 1){
        indmat <- matrix(out$ind,J,N,byrow=T)
        return(list(index=indmat, phi=phimat, mu=mumat, sigma=sigmamat))
    }
    if (out.type > 1){
        J <- dim(mumat)[1]
        d <- dim(mumat)[2]
        n <- dim(mumat)[3]

        parmat <- array(NA, dim=c(1+d+d^2,n,J))
        pmat <- phi2p.fun(phimat)
        parmat[1,,] <- t(pmat)
        parmat[1+1:d,,] <- sapply(1:J, function(i){mumat[i,,]})
        sm <- diag(d)
        sm[upper.tri(sm, diag=TRUE)] <- 1:(d*(d+1)/2)
        sm[sm == 0] <- t(sm)[sm == 0]
        parmat[-(1:(d+1)),,] <- sapply(1:J, function(i){sapply(1:n, function(k){(sigmamat[i,,k])})[as.double(sm),]})
    }
    if (out.type == 2) return(parmat)
    if (out.type > 2){
        pred <- matrix(NA, nrow(Xtest), J)
        if (printlevel > 1) cat("Making predictions ")

        SN <- ceiling(nrow(Xtest) / 5000)
        for (segno in 1:SN){ 
            seg <- seq(ceiling((segno-1)*nrow(Xtest)/SN+1e-9), segno*nrow(Xtest)/SN, by=1)
            pred[seg,] <- getpred.fun(parmat, mat=Xtest[seg,], d=d, n=n, yind=yind)
            if (printlevel > 1 & segno != SN) cat("* ")
        }
        if (printlevel > 1) cat("\n")
    }
    if (out.type == 3) return(pred)
	if (out.type == 4) return(rowMeans(pred))
}
MPDPM_MOD <- function(X, S, obs, wl, n=NULL, Burn=500, J=250, Thin=2, rho=NULL, a=1, b=1, seed="patrick",
                      iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind, out.type=4){
    ## Runs the K-Mix algorithm (in C code)
    ## X = N-by-d data matrix
    ## S = d-by-d-by-L prior covariance matrices
	## obs = observations used to create S
    ## wl = weights on prior covariance matrices
    ## Burn = burn-in iterations
    ## J = samples
    ## Thin = thinning iterations per sample
    ## rho = Wishart prior d.f. parameter (default=d+2) 
    ## a,b = beta prior parameters 
    ## seed = starting seed phrase
    ## iterprint = print progress statements if printlevel > 0
    ## printlevel = {0 NONE, 1, 2 DETAILED} for level of printout information 
    ## so.dir = character string for path of kmix.so file
    ## so.file = the kmix.so file
    ## Xtest = data matrix for the test set
    ## yind = index of response (to be predicted)
    ## out.type = output level:
    ##   1-K-Mix output parameters as list
    ##   2-K-Mix output parameters as array
    ##   3-prediction matrix for Y
    ##   4-prediction vector for Y (rowMeans of prediction matrix)
    
    ## error checking
    if (!(out.type %in% 1:4)) stop("Error: out.type must be in 1:4")
    if (is.null(n)) n <- ncol(X)*2
    if (is.null(rho)) rho <- ncol(X)+2
    if (is.null(so.dir)) so.dir <- paste(getwd(),"/code",sep="")
    if (so.dir == "pwd") so.dir <- getwd()
    if (substr(so.dir,nchar(so.dir),nchar(so.dir)) != "/")
        so.dir <-paste(so.dir, "/", sep="")
    if (is.null(so.file)) so.file <- "mpdpm.so"
    ## the difference between mpdpm.c and mpdpm2.c is multiplying the input scale
    ##   matrices by rho instead of rho+k+1, to put into form of Wishart
    ##   (rho*S is the mean and mode of the Wishart distribution, not (rho+k+1)*S...)
    if (is.null(dim(Xtest))) Xtest <- matrix(Xtest, ncol=1)
    
    ##-set up parameters for C program
    N <- nrow(X)
    k <- ncol(X)
    L <- dim(S)[3]
	
	## S smoothing
	Xm <- as.matrix(X)
	G <- (t(Xm) %*% Xm) * (nrow(X)^(-9/7))
	for (l in 1:L){
		if (obs[l] >= 1){
			S[,,l] <- (1/obs[l]) * ( (obs[l] - 0.5) * S[,,l] + 0.5 * G )
		} else {
		    S[,,l] <- 0.5 * (S[,,l] + G)	
		}				
	}

    ## format covariance matrix and weight input for C program
    Smat <- matrix(NA, k*(k+1)/2+1, L)
    for (i in 1:k)
	for (j in 1:i)
            Smat[sum(1:i)-i+j,] <- S[i,j,]
    Smat[k*(k+1)/2+1,] <- wl

    ## load C program
    dyn.load(paste(so.dir, so.file, sep=""))

    ## call C program
    out <- .C("mpdpm",
              as.integer(c(N,k,n,rho,L,Burn,Thin,J,iterprint,printlevel)),
              as.double(c(a,b)),
              as.double(stack(as.data.frame(t(X)))$values),
              as.double(stack(as.data.frame(Smat))$values),
              ind = integer(N*J),
              phi = double((n-1)*J),
              mu = double(n*k*J),
              sigma = double(n*k*(k+1)/2*J),
              as.character(seed)
              )

    ## unload C program
    dyn.unload(paste(so.dir, so.file, sep=""))

    ## reformat C output for easier use (this will eventually be done in C itself)
    ## ind is fine; split par (parameter values) into phi, mu, sigma
    phimat <- matrix(out$phi,J,n-1,byrow=T)
    muvect <- matrix(out$mu,J,n*k,byrow=T)
    sigmawidth <- k*(k+1)/2 
    sigmavect <- matrix(out$sigma,J,n*sigmawidth,byrow=T)
    mumat <- array(NA, dim=c(J,k,n))
    sigmamat <- array(NA, dim=c(J,sigmawidth,n))
    for (i in 1:n){
        mumat[,,i] <- muvect[,1:k+(i-1)*k]
        sigmamat[,,i] <- sigmavect[,1:sigmawidth+(i-1)*sigmawidth]
    }

    ## return C output to R program
    if (out.type == 1){
        indmat <- matrix(out$ind,J,N,byrow=T)
        return(list(index=indmat, phi=phimat, mu=mumat, sigma=sigmamat))
    }
    if (out.type > 1){
        J <- dim(mumat)[1]
        d <- dim(mumat)[2]
        n <- dim(mumat)[3]

        parmat <- array(NA, dim=c(1+d+d^2,n,J))
        pmat <- phi2p.fun(phimat)
        parmat[1,,] <- t(pmat)
        parmat[1+1:d,,] <- sapply(1:J, function(i){mumat[i,,]})
        sm <- diag(d)
        sm[upper.tri(sm, diag=TRUE)] <- 1:(d*(d+1)/2)
        sm[sm == 0] <- t(sm)[sm == 0]
        parmat[-(1:(d+1)),,] <- sapply(1:J, function(i){sapply(1:n, function(k){(sigmamat[i,,k])})[as.double(sm),]})
    }
    if (out.type == 2) return(parmat)
    if (out.type > 2){
        pred <- matrix(NA, nrow(Xtest), J)
        if (printlevel > 1) cat("Making predictions ")

        SN <- ceiling(nrow(Xtest) / 5000)
        for (segno in 1:SN){ 
            seg <- seq(ceiling((segno-1)*nrow(Xtest)/SN+1e-9), segno*nrow(Xtest)/SN, by=1)
            pred[seg,] <- getpred.fun(parmat, mat=Xtest[seg,], d=d, n=n, yind=yind)
            if (printlevel > 1 & segno != SN) cat("* ")
        }
        if (printlevel > 1) cat("\n")
    }
    if (out.type == 3) return(pred)
	if (out.type == 4) return(rowMeans(pred))
}
MPDPM_OLD <- function(X, S, obs=NULL, wl, n=NULL, Burn=500, J=250, Thin=2, rho=NULL, a=1, b=1, seed="patrick",
                      iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind, out.type=4){
    ## Runs the K-Mix algorithm (in C code)
    ## X = N-by-d data matrix
    ## S = d-by-d-by-L prior covariance matrices
    ## wl = weights on prior covariance matrices
    ## Burn = burn-in iterations
    ## J = samples
    ## Thin = thinning iterations per sample
    ## rho = Wishart prior d.f. parameter (default=d+2) 
    ## a,b = beta prior parameters 
    ## seed = starting seed phrase
    ## iterprint = print progress statements if printlevel > 0
    ## printlevel = {0 NONE, 1, 2 DETAILED} for level of printout information 
    ## so.dir = character string for path of kmix.so file
    ## so.file = the kmix.so file
    ## Xtest = data matrix for the test set
    ## yind = index of response (to be predicted)
    ## out.type = output level:
    ##   1-K-Mix output parameters as list
    ##   2-K-Mix output parameters as array
    ##   3-prediction matrix for Y
    ##   4-prediction vector for Y (rowMeans of prediction matrix)
    
    ## error checking
    if (!(out.type %in% 1:4)) stop("Error: out.type must be in 1:4")
    if (is.null(n)) n <- ncol(X)*2
    if (is.null(rho)) rho <- ncol(X)+2
    if (is.null(so.dir)) so.dir <- paste(getwd(),"/code",sep="")
    if (so.dir == "pwd") so.dir <- getwd()
    if (substr(so.dir,nchar(so.dir),nchar(so.dir)) != "/")
        so.dir <-paste(so.dir, "/", sep="")
    if (is.null(so.file)) so.file <- "mpdpm.so"
    ## the difference between mpdpm.c and mpdpm2.c is multiplying the input scale
    ##   matrices by rho instead of rho+k+1, to put into form of Wishart
    ##   (rho*S is the mean and mode of the Wishart distribution, not (rho+k+1)*S...)
    if (is.null(dim(Xtest))) Xtest <- matrix(Xtest, ncol=1)
    
    ##-set up parameters for C program
    N <- nrow(X)
    k <- ncol(X)
    L <- dim(S)[3]

    ## format covariance matrix and weight input for C program
    Smat <- matrix(NA, k*(k+1)/2+1, L)
    for (i in 1:k)
	for (j in 1:i)
            Smat[sum(1:i)-i+j,] <- S[i,j,]
    Smat[k*(k+1)/2+1,] <- wl

    ## load C program
    dyn.load(paste(so.dir, so.file, sep=""))

    ## call C program
    out <- .C("mpdpm",
              as.integer(c(N,k,n,rho,L,Burn,Thin,J,iterprint,printlevel)),
              as.double(c(a,b)),
              as.double(stack(as.data.frame(t(X)))$values),
              as.double(stack(as.data.frame(Smat))$values),
              ind = integer(N*J),
              phi = double((n-1)*J),
              mu = double(n*k*J),
              sigma = double(n*k*(k+1)/2*J),
              as.character(seed)
              )

    ## unload C program
    dyn.unload(paste(so.dir, so.file, sep=""))

    ## reformat C output for easier use (this will eventually be done in C itself)
    ## ind is fine; split par (parameter values) into phi, mu, sigma
    phimat <- matrix(out$phi,J,n-1,byrow=T)
    muvect <- matrix(out$mu,J,n*k,byrow=T)
    sigmawidth <- k*(k+1)/2 
    sigmavect <- matrix(out$sigma,J,n*sigmawidth,byrow=T)
    mumat <- array(NA, dim=c(J,k,n))
    sigmamat <- array(NA, dim=c(J,sigmawidth,n))
    for (i in 1:n){
        mumat[,,i] <- muvect[,1:k+(i-1)*k]
        sigmamat[,,i] <- sigmavect[,1:sigmawidth+(i-1)*sigmawidth]
    }

    ## return C output to R program
    if (out.type == 1){
        indmat <- matrix(out$ind,J,N,byrow=T)
        return(list(index=indmat, phi=phimat, mu=mumat, sigma=sigmamat))
    }
    if (out.type > 1){
        J <- dim(mumat)[1]
        d <- dim(mumat)[2]
        n <- dim(mumat)[3]

        parmat <- array(NA, dim=c(1+d+d^2,n,J))
        pmat <- phi2p.fun(phimat)
        parmat[1,,] <- t(pmat)
        parmat[1+1:d,,] <- sapply(1:J, function(i){mumat[i,,]})
        sm <- diag(d)
        sm[upper.tri(sm, diag=TRUE)] <- 1:(d*(d+1)/2)
        sm[sm == 0] <- t(sm)[sm == 0]
        parmat[-(1:(d+1)),,] <- sapply(1:J, function(i){sapply(1:n, function(k){(sigmamat[i,,k])})[as.double(sm),]})
    }
    if (out.type == 2) return(parmat)
    if (out.type > 2){
        pred <- matrix(NA, nrow(Xtest), J)
        if (printlevel > 1) cat("Making predictions ")

        SN <- ceiling(nrow(Xtest) / 5000)
        for (segno in 1:SN){ 
            seg <- seq(ceiling((segno-1)*nrow(Xtest)/SN+1e-9), segno*nrow(Xtest)/SN, by=1)
            pred[seg,] <- getpred.fun(parmat, mat=Xtest[seg,], d=d, n=n, yind=yind)
            if (printlevel > 1) cat("* ")
        }
        if (printlevel > 1) cat("\n")
    }
    if (out.type == 3) return(pred)
    if (out.type == 4) return(rowMeans(pred))
}
DPMG <- function(X, S, wl, n=NULL, kappa=1, Burn=500, J=250, Thin=2, rho=NULL, a=1, b=1, seed="patrick",
                 iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind, out.type=4){
    ## FINAL FOR S RUNS.
    ## Runs the K-Mix algorithm (in C code)
    ## X = N-by-d data matrix
    ## S = d-by-d prior covariance matrix
    ## wl = weights on prior covariance matrices
    ## Burn = burn-in iterations
    ## J = samples
    ## Thin = thinning iterations per sample
    ## rho = Wishart prior d.f. parameter (default=d+2) 
    ## a,b = beta prior parameters 
    ## seed = starting seed phrase
    ## iterprint = print progress statements if printlevel > 0
    ## printlevel = {0 NONE, 1, 2 DETAILED} for level of printout information 
    ## so.dir = character string for path of kmix.so file
    ## so.file = the kmix.so file
    ## Xtest = data matrix for the test set
    ## yind = index of response (to be predicted)
    ## out.type = output level:
    ##   1-K-Mix output parameters as list
    ##   2-K-Mix output parameters as array
    ##   3-prediction matrix for Y
    ##   4-prediction vector for Y (rowMeans of prediction matrix)
    
    ## error checking
    if (!(out.type %in% 1:4)) stop("Error: out.type must be in 1:4")
    if (is.null(n)) n <- ncol(X)*2
    if (is.null(rho)) rho <- ncol(X)+2
    if (is.null(so.dir)) so.dir <- paste(getwd(),"/code",sep="")
    if (so.dir == "pwd") so.dir <- getwd()
    if (substr(so.dir,nchar(so.dir),nchar(so.dir)) != "/")
        so.dir <-paste(so.dir, "/", sep="")
    if (is.null(so.file)) so.file <- "dpmg2.so"
    if (is.null(dim(Xtest))) Xtest <- matrix(Xtest, ncol=1)
    
    ##-set up parameters for C program
    N <- nrow(X)
    k <- ncol(X)
    L <- 1

    ## format covariance matrix and weight input for C program
    Smat <- matrix(NA, k*(k+1)/2+1, L)
    for (i in 1:k)
	for (j in 1:i)
            Smat[sum(1:i)-i+j,1] <- S[i,j]
    Smat[k*(k+1)/2+1,L] <- 1
    
    ## load C program
    dyn.load(paste(so.dir, so.file, sep=""))

    ## call C program
    out <- .C("dpmg",
              as.integer(c(N,k,n,rho,L,Burn,Thin,J,iterprint,printlevel)),
              as.double(c(a,b,kappa)),
              as.double(stack(as.data.frame(t(X)))$values),
              as.double(stack(as.data.frame(Smat))$values),
              as.double(colMeans(X)),
              ind = integer(N*J),
              phi = double((n-1)*J),
              mu = double(n*k*J),
              sigma = double(n*k*(k+1)/2*J),
              as.character(seed)
              )

    ## unload C program
    dyn.unload(paste(so.dir, so.file, sep=""))

    ## reformat C output for easier use (this will eventually be done in C itself)
    ## ind is fine; split par (parameter values) into phi, mu, sigma
    phimat <- matrix(out$phi,J,n-1,byrow=T)
    muvect <- matrix(out$mu,J,n*k,byrow=T)
    sigmawidth <- k*(k+1)/2 
    sigmavect <- matrix(out$sigma,J,n*sigmawidth,byrow=T)
    mumat <- array(NA, dim=c(J,k,n))
    sigmamat <- array(NA, dim=c(J,sigmawidth,n))
    for (i in 1:n){
        mumat[,,i] <- muvect[,1:k+(i-1)*k]
        sigmamat[,,i] <- sigmavect[,1:sigmawidth+(i-1)*sigmawidth]
    }

    ## return C output to R program
    if (out.type == 1){
        indmat <- matrix(out$ind,J,N,byrow=T)
        return(list(index=indmat, phi=phimat, mu=mumat, sigma=sigmamat))
    }
    if (out.type > 1){
        J <- dim(mumat)[1]
        d <- dim(mumat)[2]
        n <- dim(mumat)[3]

        parmat <- array(NA, dim=c(1+d+d^2,n,J))
        pmat <- phi2p.fun(phimat)
        parmat[1,,] <- t(pmat)
        parmat[1+1:d,,] <- sapply(1:J, function(i){mumat[i,,]})
        sm <- diag(d)
        sm[upper.tri(sm, diag=TRUE)] <- 1:(d*(d+1)/2)
        sm[sm == 0] <- t(sm)[sm == 0]
        parmat[-(1:(d+1)),,] <- sapply(1:J, function(i){sapply(1:n, function(k){(sigmamat[i,,k])})[as.double(sm),]})
    }
    if (out.type == 2) return(parmat)
    if (out.type > 2){
        pred <- matrix(NA, nrow(Xtest), J)
        if (printlevel > 1) cat("Making predictions ")

        SN <- ceiling(nrow(Xtest) / 5000)
        for (segno in 1:SN){ 
            seg <- seq(ceiling((segno-1)*nrow(Xtest)/SN+1e-9), segno*nrow(Xtest)/SN, by=1)
            pred[seg,] <- getpred.fun(parmat, mat=Xtest[seg,], d=d, n=n, yind=yind)
            if (printlevel > 1) cat("* ")
        }
        if (printlevel > 1) cat("\n")
    }
    if (out.type == 3) return(pred)
    if (out.type == 4) return(rowMeans(pred))
}

# KMixInit <- function(X, S, wl, initval, n=NULL, Burn=1000, J=100, Thin=10, rho=NULL, a=1, b=1, seed="patrick",
#                  iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind, out.type=4){
#     ## Runs the K-Mix algorithm (in C code)
#     ## X = N-by-d data matrix
#     ## S = d-by-d-by-L prior covariance matrices
#     ## wl = weights on prior covariance matrices
#     ## Burn = burn-in iterations
#     ## J = samples
#     ## Thin = thinning iterations per sample
#     ## rho = Wishart prior d.f. parameter (default=d+2)
#     ## a,b = beta prior parameters
#     ## seed = starting seed phrase
#     ## iterprint = print progress statements if printlevel > 0
#     ## printlevel = {0 NONE, 1, 2 DETAILED} for level of printout information
#     ## so.dir = character string for path of kmix.so file
#     ## so.file = the kmix.so file
#     ## Xtest = data matrix for the test set
#     ## yind = index of response (to be predicted)
#     ## out.type = output level:
#     ##   1-K-Mix output parameters as list
#     ##   2-K-Mix output parameters as array
#     ##   3-prediction matrix for Y
#     ##   4-prediction vector for Y (rowMeans of prediction matrix)
#
#     ## error checking
#     if (!(out.type %in% 1:4)) stop("Error: out.type must be in 1:4")
#     if (is.null(n)) n <- ncol(X)*2
#     if (is.null(rho)) rho <- ncol(X)+2
#     if (is.null(so.dir)) so.dir <- paste(getwd(),"/code",sep="")
#     if (so.dir == "pwd") so.dir <- getwd()
#     if (substr(so.dir,nchar(so.dir),nchar(so.dir)) != "/")
#         so.dir <-paste(so.dir, "/", sep="")
#     if (is.null(so.file)) so.file <- "kmixinit.so"
#     if (is.null(dim(Xtest))) Xtest <- matrix(Xtest, ncol=1)
#
#     ##-set up parameters for C program
#     N <- nrow(X)
#     k <- ncol(X)
#     L <- dim(S)[3]
#
#     ## format covariance matrix and weight input for C program
#     Smat <- matrix(NA, k*(k+1)/2+1, L)
#     for (i in 1:k)
# 	for (j in 1:i)
#             Smat[sum(1:i)-i+j,] <- S[i,j,]
#     Smat[k*(k+1)/2+1,] <- wl
#
#     ## load C program
#     dyn.load(paste(so.dir, so.file, sep=""))
#
#     ## call C program
#     out <- .C("KMixInit",
#               as.integer(c(N,k,n,rho,L,Burn,Thin,J,iterprint)),
#               as.double(c(a,b)),
#               as.double(stack(as.data.frame(t(X)))$values),
#               as.double(stack(as.data.frame(Smat))$values),
#               as.double(stack(as.data.frame(initval))$values),
#               ind = integer(N*J),
#               phi = double((n-1)*J),
#               mu = double(n*k*J),
#               sigma = double(n*k*(k+1)/2*J),
#               as.character(seed),
#               as.integer(printlevel)
#               )
#
#     ## unload C program
#     dyn.unload(paste(so.dir, so.file, sep=""))
#
#     ## reformat C output for easier use (this will eventually be done in C itself)
#     ## ind is fine; split par (parameter values) into phi, mu, sigma
#     phimat <- matrix(out$phi,J,n-1,byrow=T)
#     muvect <- matrix(out$mu,J,n*k,byrow=T)
#     sigmawidth <- k*(k+1)/2
#     sigmavect <- matrix(out$sigma,J,n*sigmawidth,byrow=T)
#     mumat <- array(NA, dim=c(J,k,n))
#     sigmamat <- array(NA, dim=c(J,sigmawidth,n))
#     for (i in 1:n){
#         mumat[,,i] <- muvect[,1:k+(i-1)*k]
#         sigmamat[,,i] <- sigmavect[,1:sigmawidth+(i-1)*sigmawidth]
#     }
#
#     ## return C output to R program
#     if (out.type == 1){
#         indmat <- matrix(out$ind,J,N,byrow=T)
#         return(list(index=indmat, phi=phimat, mu=mumat, sigma=sigmamat))
#     }
#     if (out.type > 1){
#         J <- dim(mumat)[1]
#         d <- dim(mumat)[2]
#         n <- dim(mumat)[3]
#
#         parmat <- array(NA, dim=c(1+d+d^2,n,J))
#         pmat <- phi2p.fun(phimat)
#         parmat[1,,] <- t(pmat)
#         parmat[1+1:d,,] <- sapply(1:J, function(i){mumat[i,,]})
#         sm <- diag(d)
#         sm[upper.tri(sm, diag=TRUE)] <- 1:(d*(d+1)/2)
#         sm[sm == 0] <- t(sm)[sm == 0]
#         parmat[-(1:(d+1)),,] <- sapply(1:J, function(i){sapply(1:n, function(k){(sigmamat[i,,k])})[as.double(sm),]})
#     }
#     if (out.type == 2) return(parmat)
#     if (out.type > 2){
#         pred <- matrix(NA, nrow(Xtest), J)
#         if (printlevel > 1) cat("Making predictions ")
#
#         SN <- ceiling(nrow(Xtest) / 5000)
#         for (segno in 1:SN){
#             seg <- seq(ceiling((segno-1)*nrow(Xtest)/SN+1e-9), segno*nrow(Xtest)/SN, by=1)
#             pred[seg,] <- getpred.fun(parmat, mat=Xtest[seg,], d=d, n=n, yind=yind)
#             if (printlevel > 1) cat("* ")
#         }
#         if (printlevel > 1) cat("\n")
#     }
#     if (out.type == 3) return(pred)
#     if (out.type == 4) return(rowMeans(pred))
# }
# KMixInitA <- function(X, S, wl, initval, n=NULL, Burn=1000, J=100, Thin=10, rho=NULL, a=1, b=1, seed="patrick",
#                  iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind, out.type=4){
#     ## InitA has initialization but does not use mean kernels, mean kernel centered at mean(x).
#     ## Runs the K-Mix algorithm (in C code)
#     ## X = N-by-d data matrix
#     ## S = d-by-d-by-L prior covariance matrices
#     ## wl = weights on prior covariance matrices
#     ## Burn = burn-in iterations
#     ## J = samples
#     ## Thin = thinning iterations per sample
#     ## rho = Wishart prior d.f. parameter (default=d+2)
#     ## a,b = beta prior parameters
#     ## seed = starting seed phrase
#     ## iterprint = print progress statements if printlevel > 0
#     ## printlevel = {0 NONE, 1, 2 DETAILED} for level of printout information
#     ## so.dir = character string for path of kmix.so file
#     ## so.file = the kmix.so file
#     ## Xtest = data matrix for the test set
#     ## yind = index of response (to be predicted)
#     ## out.type = output level:
#     ##   1-K-Mix output parameters as list
#     ##   2-K-Mix output parameters as array
#     ##   3-prediction matrix for Y
#     ##   4-prediction vector for Y (rowMeans of prediction matrix)
#
#     ## error checking
#     if (!(out.type %in% 1:4)) stop("Error: out.type must be in 1:4")
#     if (is.null(n)) n <- ncol(X)*2
#     if (is.null(rho)) rho <- ncol(X)+2
#     if (is.null(so.dir)) so.dir <- paste(getwd(),"/code",sep="")
#     if (so.dir == "pwd") so.dir <- getwd()
#     if (substr(so.dir,nchar(so.dir),nchar(so.dir)) != "/")
#         so.dir <-paste(so.dir, "/", sep="")
#     if (is.null(so.file)) so.file <- "kmixinit1.so"
#     if (is.null(dim(Xtest))) Xtest <- matrix(Xtest, ncol=1)
#
#     ##-set up parameters for C program
#     N <- nrow(X)
#     k <- ncol(X)
#     L <- dim(S)[3]
#
#     ## format covariance matrix and weight input for C program
#     Smat <- matrix(NA, k*(k+1)/2+1, L)
#     for (i in 1:k)
# 	for (j in 1:i)
#             Smat[sum(1:i)-i+j,] <- S[i,j,]
#     Smat[k*(k+1)/2+1,] <- wl
#
#     ## load C program
#     dyn.load(paste(so.dir, so.file, sep=""))
#
#     ## call C program
#     out <- .C("KMixInitA",
#               as.integer(c(N,k,n,rho,L,Burn,Thin,J,iterprint)),
#               as.double(c(a,b)),
#               as.double(stack(as.data.frame(t(X)))$values),
#               as.double(stack(as.data.frame(Smat))$values),
#               as.double(stack(as.data.frame(initval))$values),
#               ind = integer(N*J),
#               phi = double((n-1)*J),
#               mu = double(n*k*J),
#               sigma = double(n*k*(k+1)/2*J),
#               as.character(seed),
#               as.integer(printlevel)
#               )
#
#     ## unload C program
#     dyn.unload(paste(so.dir, so.file, sep=""))
#
#     ## reformat C output for easier use (this will eventually be done in C itself)
#     ## ind is fine; split par (parameter values) into phi, mu, sigma
#     phimat <- matrix(out$phi,J,n-1,byrow=T)
#     muvect <- matrix(out$mu,J,n*k,byrow=T)
#     sigmawidth <- k*(k+1)/2
#     sigmavect <- matrix(out$sigma,J,n*sigmawidth,byrow=T)
#     mumat <- array(NA, dim=c(J,k,n))
#     sigmamat <- array(NA, dim=c(J,sigmawidth,n))
#     for (i in 1:n){
#         mumat[,,i] <- muvect[,1:k+(i-1)*k]
#         sigmamat[,,i] <- sigmavect[,1:sigmawidth+(i-1)*sigmawidth]
#     }
#
#     ## return C output to R program
#     if (out.type == 1){
#         indmat <- matrix(out$ind,J,N,byrow=T)
#         return(list(index=indmat, phi=phimat, mu=mumat, sigma=sigmamat))
#     }
#     if (out.type > 1){
#         J <- dim(mumat)[1]
#         d <- dim(mumat)[2]
#         n <- dim(mumat)[3]
#
#         parmat <- array(NA, dim=c(1+d+d^2,n,J))
#         pmat <- phi2p.fun(phimat)
#         parmat[1,,] <- t(pmat)
#         parmat[1+1:d,,] <- sapply(1:J, function(i){mumat[i,,]})
#         sm <- diag(d)
#         sm[upper.tri(sm, diag=TRUE)] <- 1:(d*(d+1)/2)
#         sm[sm == 0] <- t(sm)[sm == 0]
#         parmat[-(1:(d+1)),,] <- sapply(1:J, function(i){sapply(1:n, function(k){(sigmamat[i,,k])})[as.double(sm),]})
#     }
#     if (out.type == 2) return(parmat)
#     if (out.type > 2){
#         pred <- matrix(NA, nrow(Xtest), J)
#         if (printlevel > 1) cat("Making predictions ")
#
#         SN <- ceiling(nrow(Xtest) / 5000)
#         for (segno in 1:SN){
#             seg <- seq(ceiling((segno-1)*nrow(Xtest)/SN+1e-9), segno*nrow(Xtest)/SN, by=1)
#             pred[seg,] <- getpred.fun(parmat, mat=Xtest[seg,], d=d, n=n, yind=yind)
#             if (printlevel > 1) cat("* ")
#         }
#         if (printlevel > 1) cat("\n")
#     }
#     if (out.type == 3) return(pred)
#     if (out.type == 4) return(rowMeans(pred))
# }
# KMixInitB <- function(X, S, wl, initval, n=NULL, Burn=1000, J=100, Thin=10, rho=NULL, a=1, b=1, seed="patrick",
#                  iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind, out.type=4){
#     ## InitB has initialization but does not use mean kernels, and uses a mean far away.
#     ## Runs the K-Mix algorithm (in C code)
#     ## X = N-by-d data matrix
#     ## S = d-by-d-by-L prior covariance matrices
#     ## wl = weights on prior covariance matrices
#     ## Burn = burn-in iterations
#     ## J = samples
#     ## Thin = thinning iterations per sample
#     ## rho = Wishart prior d.f. parameter (default=d+2)
#     ## a,b = beta prior parameters
#     ## seed = starting seed phrase
#     ## iterprint = print progress statements if printlevel > 0
#     ## printlevel = {0 NONE, 1, 2 DETAILED} for level of printout information
#     ## so.dir = character string for path of kmix.so file
#     ## so.file = the kmix.so file
#     ## Xtest = data matrix for the test set
#     ## yind = index of response (to be predicted)
#     ## out.type = output level:
#     ##   1-K-Mix output parameters as list
#     ##   2-K-Mix output parameters as array
#     ##   3-prediction matrix for Y
#     ##   4-prediction vector for Y (rowMeans of prediction matrix)
#
#     ## error checking
#     if (!(out.type %in% 1:4)) stop("Error: out.type must be in 1:4")
#     if (is.null(n)) n <- ncol(X)*2
#     if (is.null(rho)) rho <- ncol(X)+2
#     if (is.null(so.dir)) so.dir <- paste(getwd(),"/code",sep="")
#     if (so.dir == "pwd") so.dir <- getwd()
#     if (substr(so.dir,nchar(so.dir),nchar(so.dir)) != "/")
#         so.dir <-paste(so.dir, "/", sep="")
#     if (is.null(so.file)) so.file <- "kmixinit2.so"
#     if (is.null(dim(Xtest))) Xtest <- matrix(Xtest, ncol=1)
#
#     ##-set up parameters for C program
#     N <- nrow(X)
#     k <- ncol(X)
#     L <- dim(S)[3]
#
#     ## format covariance matrix and weight input for C program
#     Smat <- matrix(NA, k*(k+1)/2+1, L)
#     for (i in 1:k)
# 	for (j in 1:i)
#             Smat[sum(1:i)-i+j,] <- S[i,j,]
#     Smat[k*(k+1)/2+1,] <- wl
#
#     ## load C program
#     dyn.load(paste(so.dir, so.file, sep=""))
#
#     ## call C program
#     out <- .C("KMixInitB",
#               as.integer(c(N,k,n,rho,L,Burn,Thin,J,iterprint)),
#               as.double(c(a,b)),
#               as.double(stack(as.data.frame(t(X)))$values),
#               as.double(stack(as.data.frame(Smat))$values),
#               as.double(stack(as.data.frame(initval))$values),
#               ind = integer(N*J),
#               phi = double((n-1)*J),
#               mu = double(n*k*J),
#               sigma = double(n*k*(k+1)/2*J),
#               as.character(seed),
#               as.integer(printlevel)
#               )
#
#     ## unload C program
#     dyn.unload(paste(so.dir, so.file, sep=""))
#
#     ## reformat C output for easier use (this will eventually be done in C itself)
#     ## ind is fine; split par (parameter values) into phi, mu, sigma
#     phimat <- matrix(out$phi,J,n-1,byrow=T)
#     muvect <- matrix(out$mu,J,n*k,byrow=T)
#     sigmawidth <- k*(k+1)/2
#     sigmavect <- matrix(out$sigma,J,n*sigmawidth,byrow=T)
#     mumat <- array(NA, dim=c(J,k,n))
#     sigmamat <- array(NA, dim=c(J,sigmawidth,n))
#     for (i in 1:n){
#         mumat[,,i] <- muvect[,1:k+(i-1)*k]
#         sigmamat[,,i] <- sigmavect[,1:sigmawidth+(i-1)*sigmawidth]
#     }
#
#     ## return C output to R program
#     if (out.type == 1){
#         indmat <- matrix(out$ind,J,N,byrow=T)
#         return(list(index=indmat, phi=phimat, mu=mumat, sigma=sigmamat))
#     }
#     if (out.type > 1){
#         J <- dim(mumat)[1]
#         d <- dim(mumat)[2]
#         n <- dim(mumat)[3]
#
#         parmat <- array(NA, dim=c(1+d+d^2,n,J))
#         pmat <- phi2p.fun(phimat)
#         parmat[1,,] <- t(pmat)
#         parmat[1+1:d,,] <- sapply(1:J, function(i){mumat[i,,]})
#         sm <- diag(d)
#         sm[upper.tri(sm, diag=TRUE)] <- 1:(d*(d+1)/2)
#         sm[sm == 0] <- t(sm)[sm == 0]
#         parmat[-(1:(d+1)),,] <- sapply(1:J, function(i){sapply(1:n, function(k){(sigmamat[i,,k])})[as.double(sm),]})
#     }
#     if (out.type == 2) return(parmat)
#     if (out.type > 2){
#         pred <- matrix(NA, nrow(Xtest), J)
#         if (printlevel > 1) cat("Making predictions ")
#
#         SN <- ceiling(nrow(Xtest) / 5000)
#         for (segno in 1:SN){
#             seg <- seq(ceiling((segno-1)*nrow(Xtest)/SN+1e-9), segno*nrow(Xtest)/SN, by=1)
#             pred[seg,] <- getpred.fun(parmat, mat=Xtest[seg,], d=d, n=n, yind=yind)
#             if (printlevel > 1) cat("* ")
#         }
#         if (printlevel > 1) cat("\n")
#     }
#     if (out.type == 3) return(pred)
#     if (out.type == 4) return(rowMeans(pred))
# }
# BMix <- function(X, S, wl, n=NULL, Burn=1000, J=100, Thin=10, rho=NULL, a=1, b=1, seed="patrick",
#                  iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind=1, out.type=4){
#     ## Runs the K-Mix algorithm but recoded with only a single mean and covariance kernel
#     ## X = N-by-d data matrix
#     ## S = d-by-d-by-L prior covariance matrices
#     ## wl = weights on prior covariance matrices
#     ## Burn = burn-in iterations
#     ## J = samples
#     ## Thin = thinning iterations per sample
#     ## rho = Wishart prior d.f. parameter (default=d+1)
#     ## a,b = beta prior parameters
#     ## seed = starting seed phrase
#     ## iterprint = print progress statements if printlevel > 0
#     ## printlevel = {0 NONE, 1, 2 DETAILED} for level of printout information
#     ## so.dir = character string for path of bmix.so file
#     ## so.file = the bmix.so file
#     ## Xtest = data matrix for the test set
#     ## yind = index of response (to be predicted)
#     ## out.type = output level:
#     ##   1-K-Mix output parameters as list
#     ##   2-K-Mix output parameters as array
#     ##   3-prediction matrix for Y
#     ##   4-prediction vector for Y (rowMeans of prediction matrix)
#
#     ## error checking
#     if (!(out.type %in% 1:4)) stop("Error: out.type must be in 1:4")
#     if (is.null(n)) n <- ncol(X)*2
#     if (is.null(rho)) rho <- ncol(X)+2
#     if (is.null(so.dir)) so.dir <- paste(getwd(),"/code",sep="")
#     if (so.dir == "pwd") so.dir <- getwd()
#     if (substr(so.dir,nchar(so.dir),nchar(so.dir)) != "/")
#         so.dir <-paste(so.dir, "/", sep="")
#     if (is.null(so.file)) so.file <- "bmix.so"
#     if (is.null(dim(Xtest))) Xtest <- matrix(Xtest, ncol=1)
#
#     ##-set up parameters for C program
#     N <- nrow(X)
#     k <- ncol(X)
#     L <- dim(S)[3]
#
#     ## format covariance matrix and weight input for C program
#     ## DON'T NEED TO ALTER prior covariance matrices, bmix.c code only uses FIRST MATRIX (l=0)
#     Smat <- matrix(NA, k*(k+1)/2+1, 1)
#     for (i in 1:k)
# 	for (j in 1:i)
#             Smat[sum(1:i)-i+j,1] <- S[i,j,1]
#     Smat[k*(k+1)/2+1,1] <- wl[1]
#
#     ## load C program
#     dyn.load(paste(so.dir, so.file, sep=""))
#
#     ## call C program
#     out <- .C("BMix",
#               as.integer(c(N,k,n,rho,1,Burn,Thin,J,iterprint)),
#               as.double(c(a,b)),
#               as.double(stack(as.data.frame(t(X)))$values),
#               as.double(stack(as.data.frame(Smat))$values),
#               as.double(colMeans(X)),
#               ind = integer(N*J),
#               phi = double((n-1)*J),
#               mu = double(n*k*J),
#               sigma = double(n*k*(k+1)/2*J),
#               as.character(seed),
#               as.integer(printlevel)
#               )
#
#     ## unload C program
#     dyn.unload(paste(so.dir, so.file, sep=""))
#
#     ## reformat C output for easier use
#     ## ind is fine; split par (parameter values) into phi, mu, sigma; keep M for checking purposes only
#     phimat <- matrix(out$phi,J,n-1,byrow=T)
#     muvect <- matrix(out$mu,J,n*k,byrow=T)
#     sigmawidth <- k*(k+1)/2
#     sigmavect <- matrix(out$sigma,J,n*sigmawidth,byrow=T)
#     mumat <- array(NA, dim=c(J,k,n))
#     sigmamat <- array(NA, dim=c(J,sigmawidth,n))
#     for (i in 1:n){
#         mumat[,,i] <- muvect[,1:k+(i-1)*k]
#         sigmamat[,,i] <- sigmavect[,1:sigmawidth+(i-1)*sigmawidth]
#     }
#
#     ## return C output to R program
#     if (out.type == 1){
#         indmat <- matrix(out$ind,J,N,byrow=T)
#         return(list(index=indmat, phi=phimat, mu=mumat, sigma=sigmamat))
#     }
#     if (out.type > 1){
#         J <- dim(mumat)[1]
#         d <- dim(mumat)[2]
#         n <- dim(mumat)[3]
#
#         parmat <- array(NA, dim=c(1+d+d^2,n,J))
#         pmat <- phi2p.fun(phimat)
#         parmat[1,,] <- t(pmat)
#         parmat[1+1:d,,] <- sapply(1:J, function(i){mumat[i,,]})
#         sm <- diag(d)
#         sm[upper.tri(sm, diag=TRUE)] <- 1:(d*(d+1)/2)
#         sm[sm == 0] <- t(sm)[sm == 0]
#         parmat[-(1:(d+1)),,] <- sapply(1:J, function(i){sapply(1:n, function(k){(sigmamat[i,,k])})[as.double(sm),]})
#     }
#     if (out.type == 2) return(parmat)
#     if (out.type > 2){
#         pred <- matrix(NA, nrow(Xtest), J)
#         if (printlevel > 1) cat("Making predictions ")
#
#         SN <- ceiling(nrow(Xtest) / 5000)
#         for (segno in 1:SN){
#             seg <- seq(ceiling((segno-1)*nrow(Xtest)/SN+1e-9), segno*nrow(Xtest)/SN, by=1)
#             pred[seg,] <- getpred.fun(parmat, mat=Xtest[seg,], d=d, n=n, yind=yind)
#             if (printlevel > 1) cat("* ")
#         }
#         if (printlevel > 1) cat("\n")
#     }
#     if (out.type == 3) return(pred)
#     if (out.type == 4) return(rowMeans(pred))
# }
# BMixInit <- function(X, S, wl, initval, n=NULL, Burn=1000, J=100, Thin=10, rho=NULL, a=1, b=1, seed="patrick",
#                  iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL, Xtest, yind=1, out.type=4){
#     ## Runs the K-Mix algorithm but recoded with only a single mean and covariance kernel
#     ## X = N-by-d data matrix
#     ## S = d-by-d-by-L prior covariance matrices
#     ## wl = weights on prior covariance matrices
#     ## par = parameter matrix values (columns are component values, rows are p, mu[1:k], sigma[1:kk])
#     ## Burn = burn-in iterations
#     ## J = samples
#     ## Thin = thinning iterations per sample
#     ## rho = Wishart prior d.f. parameter (default=d+1)
#     ## a,b = beta prior parameters
#     ## seed = starting seed phrase
#     ## iterprint = print progress statements if printlevel > 0
#     ## printlevel = {0 NONE, 1, 2 DETAILED} for level of printout information
#     ## so.dir = character string for path of bmix.so file
#     ## so.file = the bmix.so file
#     ## Xtest = data matrix for the test set
#     ## yind = index of response (to be predicted)
#     ## out.type = output level:
#     ##   1-K-Mix output parameters as list
#     ##   2-K-Mix output parameters as array
#     ##   3-prediction matrix for Y
#     ##   4-prediction vector for Y (rowMeans of prediction matrix)
#
#     ## error checking
#     if (!(out.type %in% 1:4)) stop("Error: out.type must be in 1:4")
#     if (is.null(n)) n <- ncol(X)*2
#     if (is.null(rho)) rho <- ncol(X)+2
#     if (is.null(so.dir)) so.dir <- paste(getwd(),"/code",sep="")
#     if (so.dir == "pwd") so.dir <- getwd()
#     if (substr(so.dir,nchar(so.dir),nchar(so.dir)) != "/")
#         so.dir <-paste(so.dir, "/", sep="")
#     if (is.null(so.file)) so.file <- "bmixinit.so"
#     if (is.null(dim(Xtest))) Xtest <- matrix(Xtest, ncol=1)
#
#     ##-set up parameters for C program
#     N <- nrow(X)
#     k <- ncol(X)
#     L <- dim(S)[3]
#
#     ## format covariance matrix and weight input for C program
#     ## DON'T NEED TO ALTER prior covariance matrices, bmix.c code only uses FIRST MATRIX (l=0)
#     Smat <- matrix(NA, k*(k+1)/2+1, 1)
#     for (i in 1:k)
# 	for (j in 1:i)
#             Smat[sum(1:i)-i+j,1] <- S[i,j,1]
#     Smat[k*(k+1)/2+1,1] <- wl[1]
#
#     ## load C program
#     dyn.load(paste(so.dir, so.file, sep=""))
#
#     ## call C program
#     out <- .C("BMixInit",
#               as.integer(c(N,k,n,rho,1,Burn,Thin,J,iterprint)),
#               as.double(c(a,b)),
#               as.double(stack(as.data.frame(t(X)))$values),
#               as.double(stack(as.data.frame(Smat))$values),
#               as.double(colMeans(X)),
#               as.double(stack(as.data.frame(initval))$values),
#               ind = integer(N*J),
#               phi = double((n-1)*J),
#               mu = double(n*k*J),
#               sigma = double(n*k*(k+1)/2*J),
#               as.character(seed),
#               as.integer(printlevel)
#               )
#
#     ## unload C program
#     dyn.unload(paste(so.dir, so.file, sep=""))
#
#     ## reformat C output for easier use
#     ## ind is fine; split par (parameter values) into phi, mu, sigma; keep M for checking purposes only
#     phimat <- matrix(out$phi,J,n-1,byrow=T)
#     muvect <- matrix(out$mu,J,n*k,byrow=T)
#     sigmawidth <- k*(k+1)/2
#     sigmavect <- matrix(out$sigma,J,n*sigmawidth,byrow=T)
#     mumat <- array(NA, dim=c(J,k,n))
#     sigmamat <- array(NA, dim=c(J,sigmawidth,n))
#     for (i in 1:n){
#         mumat[,,i] <- muvect[,1:k+(i-1)*k]
#         sigmamat[,,i] <- sigmavect[,1:sigmawidth+(i-1)*sigmawidth]
#     }
#
#     ## return C output to R program
#     if (out.type == 1){
#         indmat <- matrix(out$ind,J,N,byrow=T)
#         return(list(index=indmat, phi=phimat, mu=mumat, sigma=sigmamat))
#     }
#     if (out.type > 1){
#         J <- dim(mumat)[1]
#         d <- dim(mumat)[2]
#         n <- dim(mumat)[3]
#
#         parmat <- array(NA, dim=c(1+d+d^2,n,J))
#         pmat <- phi2p.fun(phimat)
#         parmat[1,,] <- t(pmat)
#         parmat[1+1:d,,] <- sapply(1:J, function(i){mumat[i,,]})
#         sm <- diag(d)
#         sm[upper.tri(sm, diag=TRUE)] <- 1:(d*(d+1)/2)
#         sm[sm == 0] <- t(sm)[sm == 0]
#         parmat[-(1:(d+1)),,] <- sapply(1:J, function(i){sapply(1:n, function(k){(sigmamat[i,,k])})[as.double(sm),]})
#     }
#     if (out.type == 2) return(parmat)
#     if (out.type > 2){
#         pred <- matrix(NA, nrow(Xtest), J)
#         if (printlevel > 1) cat("Making predictions ")
#
#         SN <- ceiling(nrow(Xtest) / 5000)
#         for (segno in 1:SN){
#             seg <- seq(ceiling((segno-1)*nrow(Xtest)/SN+1e-9), segno*nrow(Xtest)/SN, by=1)
#             pred[seg,] <- getpred.fun(parmat, mat=Xtest[seg,], d=d, n=n, yind=yind)
#             if (printlevel > 1) cat("* ")
#         }
#         if (printlevel > 1) cat("\n")
#     }
#     if (out.type == 3) return(pred)
#     if (out.type == 4) return(rowMeans(pred))
# }
# KMix.mod <- function(X, S, wl, n=NULL, Burn=1000, J=100, Thin=10, rho=NULL, a=1, b=1,
#                      seed="patrick", iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL,
#                      Xtest, yind=1, out.type=4, ntry=NULL, cvfolds=4, seedno=2014){
#     set.seed(seedno)
#     if (is.null(ntry)) ntry <- ncol(X)
#     startseed <- seed
#
#     res <- rep(NA, ntry)
#     for (tryno in 1:ntry){
#         out <- KMix(X=X, S=S, wl=wl, n=n, Burn=250, J=50, Thin=5, rho=rho, a=a, b=b,
#                         seed=paste(startseed,tryno,sep=""), iterprint=10000, printlevel=0,
#                         so.dir=so.dir, so.file=so.file, Xtest=X, yind=yind, out.type=4)
#         res[tryno] <- sum((out - X[,yind])^2)
#     }
#     if (ntry > 1) tryvar <- sd(res)/mean(res) else tryvar <- NA
#     out <- KMix(X=X, S=S, wl=wl, n=n, Burn=Burn, J=J, Thin=Thin, rho=rho, a=a, b=b,
#                     seed=paste(startseed,which.min(res),sep=""),
#                     iterprint=iterprint, printlevel=printlevel,
#                     so.dir=so.dir, so.file=so.file, Xtest=Xtest, yind=yind, out.type=out.type)
#     return(list(out=out, var=tryvar, ind=which.min(res)))
# }
# BMix.mod <- function(X, S, wl, n=NULL, Burn=1000, J=100, Thin=10, rho=NULL, a=1, b=1,
#                      seed="patrick", iterprint=1000, printlevel=1, so.dir=NULL, so.file=NULL,
#                      Xtest, yind=1, out.type=4, ntry=NULL, cvfolds=4, seedno=2014){
#     set.seed(seedno)
#     if (is.null(ntry)) ntry <- ncol(X)
#     startseed <- seed
#
#     res <- rep(NA, ntry)
#     for (tryno in 1:ntry){
#         out <- BMix(X=X, S=S, wl=wl, n=n, Burn=250, J=50, Thin=5, rho=rho, a=a, b=b,
#                         seed=paste(startseed,tryno,sep=""), iterprint=10000, printlevel=0,
#                         so.dir=so.dir, so.file=so.file, Xtest=X, yind=yind, out.type=4)
#         res[tryno] <- sum((out - X[,yind])^2)
#     }
#     if (ntry > 1) tryvar <- sd(res)/mean(res) else tryvar <- NA
#     out <- BMix(X=X, S=S, wl=wl, n=n, Burn=Burn, J=J, Thin=Thin, rho=rho, a=a, b=b,
#                     seed=paste(startseed,which.min(res),sep=""),
#                     iterprint=iterprint, printlevel=printlevel,
#                     so.dir=so.dir, so.file=so.file, Xtest=Xtest, yind=yind, out.type=out.type)
#     return(list(out=out, var=tryvar, ind=which.min(res)))
# }
