#' Internal helper function for creating default L matrix if not specified
#'
#' @param xtar Predictor matrix from the target population
#' @param family One of 'gaussian', 'binomial', or 'cox'
#'
#'
#' @return A square matrix with dimension equal to ncol(xtar)
#' @export
#'
#' @examples

make_L <- function(xtar,
                   family) {

  if(family != "cox") {
    L <- cbind(0, rbind(0, solve(base::chol(crossprod(xtar)))))
    L[1,1] <- 1
  } else if(family == "cox") {
    L <- solve(base::chol(crossprod(xtar)))
  }

  L
}
