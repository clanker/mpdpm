#' @title Fit MPDPM model.
#' @description Fits MPDPM model to data.
#' @export
#' @return Predicted values on testing data.
#' @param split_data Initial_split of data.
#' @param mclust_fit mclust run of training data.
#' @examples
model_mpdpm <- function(split_data, mclust_fit, ...) {
  samp <- rsample::training(split_data)
  N <- nrow(samp)
  p <- ncol(samp)
  rho <- p + 3
  K <- N / (rho + p/2 + 1/2)
  out <- MPDPM(as.matrix(samp), mclust_fit$S, 
               mclust_fit$obs * 1000, mclust_fit$wl, 
               n = ceiling(K), Burn = 1000, J = 400, Thin = 10, rho = rho, 
               a = K ^ (1/3), b = K, seed = "patrick", iterprint = 1000, 
               printlevel = 2, so.dir = './script/c-code/', so.file = 'mpdpm.so', 
               Xtest = as.matrix(rsample::testing(split_data)), yind = 1, out.type = 4)  
  out
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