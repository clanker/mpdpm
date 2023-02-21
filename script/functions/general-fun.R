calc_regression_metrics <- function(pred_mat) {
  # pred_mat has two columns: estimate and truth
  SST <- sum((pred_mat$truth - mean(pred_mat$truth))^2)
  SSE <- sum((pred_mat$truth - pred_mat$estimate)^2)
  R2 <- 1 - SSE / SST
  RMSE <- sqrt(SSE / nrow(pred_mat))
  AAE <- mean(abs(pred_mat$truth - pred_mat$estimate))
  return(list(Rsquared = R2, RMSE = RMSE, AAE = AAE))
}


#' @title Calculate performance
#' @description Return accuracy performance of testing dataset.
#' @export
#' @return Vector of performance characteristics
#' @param split_data Initial_split of data
#' @param predictions Prediction values
# @param method Name of method that generated prediction values.
#' @examples
calc_accuracy <- function(split_data, predictions) {
  test_truth <- split_data |>
    rsample::testing() |>
    dplyr::pull(y)
  estimates_tbl <- tibble::tibble(
    truth = test_truth,
    estimate = predictions
  )
  estimates_tbl |>
    calc_regression_metrics() |>
    tibble::as_tibble()
}


phi2p.fun <- function(phi){
  ## converts phi variables to mixture proportions p
  phi <- cbind(phi, 1)
  k <- ncol(phi)
  phi2 <- 1 - phi
  for (j in 2:k){
    phi2[ , j] <- phi2[ , j - 1] * phi2[ , j]
    phi[ , j] <- phi2[ , j - 1] * phi[ , j]
  }	
  return(phi)
}


indivlik.fun <- function(par, mat, d, xind, yind){
  ## likelihood calculation
  return(par[1] * mvtnorm::dmvnorm(mat, 
                                   mean = par[1+xind],
                                   sigma = as.matrix(matrix(par[-(1:(d + 1))], d, d)[-yind, -yind])))
}


indivpred.fun <- function(par, mat, d, xind, yind){
  ## predicted value calculation
  sig <- matrix(par[-(1:(d + 1))], d, d)
  return(sig[yind, -yind] %*% solve(sig[-yind, -yind]) %*% (t(mat) - par[1 + xind]) + par[1 + yind])
}


iterpred.fun <- function(par2, mat, d, xind, yind){
  ## get prediction per iteration
  lik <- apply(par2, 2, indivlik.fun, mat = as.matrix(mat), d = d, xind = xind, yind = yind)
  lik[, 1] <- lik[, 1] + 1e-300
  pred <- apply(par2, 2, indivpred.fun, mat = mat, d = d, xind = xind, yind = yind)
  rowSums( pred * lik * (1 / rowSums(lik)) )
}


getpred.fun <- function(parmat, mat, d = NULL, n = NULL, yind = NULL){
  ## get overall prediction: the mean of all iteration predicted values
  if (is.null(d)) {d <- ncol(mat)}
  if (is.null(yind)) {yind <- d}
  xind <- setdiff(1:d, yind)
  apply(parmat[ , 1:n, ], 3, iterpred.fun, mat = mat[, xind], d = d, xind = xind, yind = yind)
}
