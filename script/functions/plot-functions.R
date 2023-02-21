# Functions
plotsize <- function(wd,ht){options(repr.plot.width=wd, repr.plot.height=ht)}


# based on ~/Dropbox/Rmodels/mvn_plots.R
out2dplots <- function(train, out, delay=1){
  N=1000
  x3 = c(seq(1,0,length=N/4)^.333,rep(0,3*N/4))
  x2 = c(seq(0,1,length=N/4)^.667,rep(1,N/4),seq(1,0,length=N/2))
  x1 = c(rep(0,N/4),sqrt(seq(0,1,length=N/4)),rep(1,N/2))
  nclust = max(out$index)+1
  maxmem2=nrow(train)
  for (i in 1:nrow(out$index)){
    plot(train[,2], train[,1], xlim=c(-8,8), ylim=c(-8,8))
    #main=paste("sample #", i))
    for (j in 1:nclust){
      mem = sum(out$index[i,] == j-1)
      nc2 = round(N*(j-0.5)/nclust)
      col = rgb(x1[nc2],x2[nc2],x3[nc2])
      if (mem == 0)
        col = 'grey70'
        meanvect = out$mu[i,c(1,2),j]
        sigmat <- diag(dim(out$mu)[2])
        sigmat[upper.tri(sigmat, diag=TRUE)] <- as.double(out$sigma[i,,j])
        sigmat <- sigmat + t(sigmat) - diag(diag(sigmat))
        ellipse.ci(meanvect, sigmat[1:2,1:2], p=.9, col=col, lwd=0.5+sqrt(mem))
        isx = which(out$index[i,] == j-1)
        if (length(isx) > 0){
          points(train[isx,2], train[isx,1], cex=.75, col=col, pch=16)
        }
    }
    Sys.sleep(delay)
  }
}


library(plotrix)
ellipse.ci <- function(dat, mat, p, col=2, lwd=2){
  g=svd(mat)
  d=atan(g$u[2, 1] / g$u[1, 1])*180/pi
  q=qnorm(1-(1-p)/2)
  draw.ellipse(dat[1], dat[2],
               q*sqrt(pi/2*g$d[1]), q*sqrt(pi/2*g$d[2]),
               angle=d, border=col, lwd=lwd)
}
## This function was in sigmatest_v1a.Rmd
# ellipse.ci <- function(dat, mat, p, col=2, lwd=2){
#   g=svd(mat)
#   d=atan(g$u[2,1]/g$u[2,2])*180/pi
#   q=qnorm(1-(1-p)/2)
#   draw.ellipse(dat[2],dat[1],
#                q*sqrt(pi/2*g$d[2]),q*sqrt(pi/2*g$d[1]),
#                angle=d, border=col, lwd=lwd)		
# }


get.rmse <- function(y, yhat, X, xs){
  sst = sum((y-mean(y))^2)
  sse = sum((y-yhat)^2)
  v1 = floor(min(y)); v2 = ceiling(max(y))
  contour(xs,xs,matrix(y,L,L),levels = seq(v1,v2,len=v2-v1+1))
  text(X[,1],X[,2],round(X[,3],1),col='red')
  v1 = floor(min(yhat)); v2 = ceiling(max(yhat))
  contour(xs,xs,matrix(yhat,L,L),levels = seq(v1,v2,len=v2-v1+1),col='blue',lty=2,add = T)
  return(list(rmse=sqrt(sum((y-yhat)^2)/length(y)), r2=1-sse/sst))
}