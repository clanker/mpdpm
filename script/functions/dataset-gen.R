#' @title Generate data.
#' @description Generates Gaussian mixture data for DPMG testing.
#' @export
#' @return A tibble of three columns (Y, X1, X2)
#' @param
#' @examples
#' main_data <- generate_data()
generate_data <- function(seed) {
  N <- 250  # number of points per group
  out <- matrix(0, 0, 3)
  cov_mat <- diag(3) / 4
  cov_mat[1, 1] <- 1/16
  set.seed(seed)
  for (x1 in c(-1, 1)) {
    for (x2 in c(-1, 1)) {
      mu_vec <- c(x1 * x2, x1, x2)
      out <- rbind(out, mvtnorm::rmvnorm(N, mu_vec, cov_mat))
    }
  }
  colnames(out) <- c("y", "x1", "x2")
  return(tibble::as_tibble(out))
}


plot_data_raw <- function(dat) {
  dat |>
    ggplot(aes(x = x1, y = x2, color = y)) +
    geom_point(alpha = 1/3, size = 1)
}


#' @title Split data.
#' @description Split data into training and testing.
#' @export
#' @return A rsample initial_split object.
#' @param data Tibble of data
#' @examples
split_data <- function(data, seed = 0) {
  if (!is.null(seed)) { set.seed(seed) }
  return(data |> rsample::initial_split(prop = 0.1))
}


fit_model <- function(split_data, mclust_fit, param, seed, method) {
  fxn <- eval(parse(text = paste0("model_", method)))
  predictions <- fxn(split_data, mclust_fit, numerator_diff = param)
  calc_accuracy(split_data, predictions) |>
    dplyr::mutate(numerator_diff = param, seed = seed, method = method)
}
