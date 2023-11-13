#' Safe function for accurately calculating log(1 + exp(x)) even for very large and small x's
#'
#' @param x Vector of one or more numbers
#'
#' @return A vector of equal length
#' @export
#'
#' @examples

log1plex = function(x) {
  - stats::plogis(-x, log.p = TRUE)
}
