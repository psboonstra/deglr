#' Calculate the Cox partial log-likelihood
#'
#' @param y A vector of class Surv()
#' @param x A vector with length equal to that of y or matrix with number of rows
#' equal to that of y. Each element corresponds to a linear predictor.
#'
#' @return If x is a vector, then a single number corresponding to the partial
#' log-likelihood. If x is a matrix, then a vector corresponding to the partial
#' log-likelihoods for each column in x.
#'
#' @importFrom survival coxph
#' @importFrom survival coxph.control
#' @importFrom purrr map_dbl
#'
#' @export
#'
#' @examples
#'
coxph_ploglik <- function(y, x) {

  if("numeric" %in% class(x) && length(x) == length(y)) {
    sum(coxph(y ~ x, init = c(1), control = coxph.control(iter.max = 0))$loglik)
  } else if("matrix" %in% class(x) && nrow(x) == length(y)) {
    map_dbl(asplit(x, 2), ~sum(coxph(y ~ .x, init = 1, control = coxph.control(iter.max = 0))$loglik))
  } else {
    stop("'x' must either be a numeric vector of equal length as 'y' or 'x' must be a matrix with number of rows equal to the length of 'y'")
  }
}

