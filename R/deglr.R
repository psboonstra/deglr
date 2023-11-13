#' Fit data-enriched generalized linear regression model
#'
#' @param xtar Predictor matrix from the target population
#' @param ytar Outcome vector from the target population
#' @param xext Predictor matrix from the external population
#' @param yext Outcome vector from the external population
#' @param family One of 'gaussian', 'binomial', or 'cox'
#' @param alpha See `glmnet::cv.glmnet()`
#' @param nlambda See glmnet::glmnet() but note that in contrast to glmnet() the
#'   sequence of lambdas will always be exactly this long
#' @param lambda.min.ratio See glmnet::glmnet()
#' @param lambda User can provide a fixed sequence of lambdas instead
#'   of an automatically constructed sequence. If provided, this will override
#'   `nlambda` and `lambda.min.ratio`
#' @param return_what This argument allows the user to separately indicate that
#'   they only want results returned corresponding to specific values of lambda
#'   or that they only want to calculate an appropriate sequence of lambdas.
#'   This can take on one of three values: `NULL` (its default), or a numeric
#'   vector of non-negative numbers, or anything else. If `NULL`, then the
#'   function will return the entire solution path. If a numeric vector, then
#'   the function will return the solution path only for these values of lambda.
#'   If anything else, then the function will return a sequence of lambdas and
#'   nothing else.
#' @param L Square matrix that scales the bias parameter gamma. If `NULL`, then
#'   the `make_L()` function will be called to set `L` as the cholesky
#'   decomposition of `crossprod(xtar)`
#' @param standardize If TRUE, then xext will be standardized but only in the
#'   part of the augmented data matrix corresponding to gamma.
#' @param thresh See glmnet::glmnet()
#' @param maxit See glmnet::glmnet()
#'
#' @importFrom glmnet glmnet
#' @importFrom glmnet coef.glmnet
#' @return A named list
#' @export
#'
#' @examples

deglr <- function(xtar,
                  ytar,
                  xext,
                  yext,
                  family,
                  alpha = 0,
                  nlambda = 20,
                  lambda.min.ratio = ifelse(nrow(xtar) < d/2, 1e-3, 1e-6),
                  lambda = NULL,
                  return_what = NULL,
                  L = NULL,
                  standardize = TRUE,
                  thresh = 1e-14,
                  maxit = 1e5){

  d <- ncol(xtar)
  beta_index <- 1:(d + I(family != "cox"))
  if(ncol(xtar) != ncol(xext)) {
    stop("'xtar' and 'xext' must have same numbers of columns")
  }

  if(!family %in% c("gaussian", "binomial", "cox")) {
    stop("'family' must be one of 'gaussian', 'binomial', or 'cox'")
  }

  if(is.null(L)) {
    L <- make_L(xtar, family)
  } else if(any(dim(L) != (ncol(xtar) + (family != "cox")))) {
    stop("'L' must be a square matrix with dimension equal to the number of columns in 'xtar' and 'xext' plus 1 (if 'family' is not 'cox') or plus 0 (if 'family' is 'cox')")
  }

  foo <-
    augment_data(xtar = xtar,
                 ytar = ytar,
                 xext = xext,
                 yext = yext,
                 family = family,
                 L = L,
                 standardize = standardize)
  y_aug = foo$y_aug
  x_aug = foo$x_aug
  xext_star_scale = foo$xext_star_scale

  if(is.null(lambda)) {
    # Run this solely to get the max lambda:
    mod <- glmnet(x = x_aug,
                  y = y_aug,
                  family = family,
                  alpha = alpha,
                  nlambda = nlambda,
                  #Note: lambda.min.ratio = 0.5 is not a typo
                  lambda.min.ratio = 0.5,
                  #Note: lambda = NULL is not a typo
                  lambda = NULL,
                  penalty.factor = c(rep(0, d), rep(1, d + (family != "cox"))),
                  # Standardization taken care of in parent
                  standardize = FALSE)
    # We don't use the default lambda sequence from glmnet because it stops
    # too soon, i.e. doesn't give small enough values of lambda
    max_log_lambda <- log(max(mod$lambda))
    lambda_seq <- exp(seq(from = max_log_lambda, to = max_log_lambda + log(lambda.min.ratio), length = nlambda))
  } else {
    lambda_seq <- lambda
  }

  if(is.null(return_what)){
    return_what <- numeric(0)
    which_lambda <- seq_along(lambda_seq);
  } else if(is.numeric(return_what)) {
    lambda_seq <- rev(sort(unique(c(lambda_seq, return_what))))
    which_lambda <- match(return_what, lambda_seq)
  }

  if(is.numeric(return_what)) {
    mod2 <- glmnet(x = x_aug,
                   y = y_aug,
                   family = family,
                   alpha = alpha,
                   lambda = lambda_seq,
                   penalty.factor = c(rep(0, d), rep(1, d + (family != "cox"))),
                   # Standardization taken care of above
                   standardize = FALSE,
                   thresh = thresh,
                   maxit = maxit)


    beta_hat <- coef.glmnet(mod2)[beta_index, which_lambda]
    gamma_hat_star <- coef.glmnet(mod2)[-beta_index, which_lambda]

    gamma_hat <- drop(L %*% diag(1 / xext_star_scale) %*% gamma_hat_star)
    return( list( lambda = lambda_seq[which_lambda], beta_hat = beta_hat, gamma_hat = gamma_hat, glmnet = mod2) )

  } else {
    return(list(lambda = lambda_seq))
  }

}








