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
