#' Internal helper function for reformatting external and target data
#'
#' @param xtar Predictor matrix from the target population
#' @param ytar Outcome vector from the target population
#' @param xext Predictor matrix from the external population
#' @param yext Outcome vector from the external population
#' @param L Square matrix that scales the bias parameter gamma. The default
#'   value is to use inverse of the cholesky decomposition of crossprod(xtar)
#' @param standardize If TRUE, then xext will be standardized but only in the
#'   part of the augmented data matrix corresponding to gamma.
#'
#' @return A named list to be used by deglr
#'
#' @examples
augment_data <- function(xtar,
                         ytar,
                         xext,
                         yext,
                         L = NULL,
                         standardize = TRUE){

  stopifnot(ncol(xtar) == ncol(xext))
  if(is.null(L)) {
    L <- cbind(0, rbind(0, solve(base::chol(crossprod(xtar)))))
    L[1,1] <- 1
  }

  xext_star <- cbind(1, xext) %*% L

  if (standardize) {
    xext_star_scale = c(1, glmnet:::weighted_mean_sd(xext_star[,-1])$sd)
  } else {
    xext_star_scale = rep(1, ncol(xext_star))
  }

  xext_star <- scale(xext_star, center = FALSE, scale = xext_star_scale)

  x_aug <- cbind(rbind(xtar, xext),
                 rbind(matrix(0, nrow = nrow(xtar), ncol = ncol(xtar) + 1), xext_star))

  return( list(y_aug = c(ytar, yext),
               x_aug = x_aug,
               xext_star_scale = xext_star_scale,
               L = L) )
}
