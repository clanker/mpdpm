#' @title Fit mclust to data.
#' @description Fits mclust and returns fit information needed for MPDPM.
#' @export
#' @return list of covariance matrices (S), weights (wl), and obs. per matrix (obs)
#' @param split_data Initial_split of data.
#' @param mclust_fit mclust run of training data.
#' @examples
fit_mclust <- function(split_data) {
  samp <- rsample::training(split_data)
  outcov <- run.mclust_v3(samp, sampsize = nrow(samp), seedno = 2021, 
                          cmax = 8, cmin = 1, printlevel = 0, noIflag = TRUE)
  S1 = outcov$mat
  wl1 = rep(1, dim(S1)[3])
  obs1 = round(outcov$obs, 2)
  obs1[obs1 == 1] = 0
  return(list(S = S1, wl = wl1, obs = obs1))
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
  
  ranklist <- (1:length(biclist))[order(biclist, decreasing=T)]
  rankn <- min(length(ranklist),8)
  if (printlevel > 0){
    cat("For seedno", seedno, "mclust finished (", length(biclist),
        ")  Best bic groups are: ", ranklist[1:rankn], "\n")
  }
  return(list(mat = priormat, wt = priorwt, bic = biclist, obs = obs, prob=prokeep, n = sampsize, rl=ranklist))
}


## trymclust functions ---------------------------------------------------------------
##
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