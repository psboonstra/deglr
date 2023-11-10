#' Internal helper function for reformatting external and target data
#'
#' @param xtar Predictor matrix from the target population
#' @param ytar Outcome vector from the target population
#' @param xext Predictor matrix from the external population
#' @param yext Outcome vector from the external population
#' @param family One of 'gaussian', 'binomial', or 'cox'
#' @param L Square matrix that scales the bias parameter gamma. The default
#'   value is to use inverse of the cholesky decomposition of crossprod(xtar)
#' @param standardize If TRUE, then xext will be standardized but only in the
#'   part of the augmented data matrix corresponding to gamma.
#'
#' @importFrom glmnet stratifySurv
#'
#' @return A named list to be used by deglr
#'
#' @examples
augment_data <- function(xtar,
                         ytar,
                         xext,
                         yext,
                         family,
                         L = NULL,
                         standardize = TRUE){

  if(is.null(L)) {
    L <- make_L(xtar, family)
  }

  if(family != "cox") {

    xext_star <- cbind(1, xext) %*% L

    if (standardize) {
      xext_star_scale = c(1, glmnet:::weighted_mean_sd(xext_star[,-1])$sd)
    } else {
      xext_star_scale = rep(1, ncol(xext_star))
    }

    xext_star <- scale(xext_star, center = FALSE, scale = xext_star_scale)

    x_aug <- cbind(rbind(xtar, xext),
                   rbind(matrix(0, nrow = nrow(xtar), ncol = ncol(xtar) + 1), xext_star))
    y_aug = c(ytar, yext)

  } else {

    xext_star <- xext %*% L

    if (standardize) {
      xext_star_scale = glmnet:::weighted_mean_sd(xext_star)$sd
    } else {
      xext_star_scale = rep(1, ncol(xext_star))
    }

    xext_star <- scale(xext_star, center = FALSE, scale = xext_star_scale)

    x_aug <- cbind(rbind(xtar, xext),
                   rbind(matrix(0, nrow = nrow(xtar), ncol = ncol(xtar)), xext_star))
    y_aug = stratifySurv(c(ytar, yext), c(rep(1, length(ytar)), rep(2, length(yext))))

  }

  return( list(y_aug = y_aug,
               x_aug = x_aug,
               xext_star_scale = xext_star_scale,
               L = L) )
}
