#' Fit data-enriched generalized linear regression model
#'
#' @param xtar Predictor matrix from the target population
#' @param ytar Outcome vector from the target population
#' @param xext Predictor matrix from the external population
#' @param yext Outcome vector from the external population
#' @param family One of 'gaussian', 'binomial', or 'cox'
#' @param nlambda See glmnet::glmnet() but note that in contrast to glmnet() the
#'   sequence of lambdas will always be exactly this long
#' @param lambda.min.ratio See glmnet::glmnet()
#' @param fixed_lambda_seq User can provide a fixed sequence of lambdas instead
#'   of an automatically constructed sequence. If provided, this will override
#'   `nlambda` and `lambda.min.ratio`
#' @param return_only_at A vector of lambdas that may or may not have overlap
#'   with `fixed_lambda_seq.` This argument allows the user to separately indicate
#'   that they only want results returned corresponding to specific values of
#'   lambda.
#' @param L Square matrix that scales the bias parameter gamma. The default
#'   value is to use inverse of the cholesky decomposition of crossprod(xtar)
#' @param standardize If TRUE, then xext will be standardized but only in the
#'   part of the augmented data matrix corresponding to gamma.
#' @param thresh See glmnet::glmnet()
#' @param maxit See glmnet::glmnet()
#'
#' @importFrom glmnet glmnet
#' @importFrom glmnet coef.glmnet
#' @return
#' @export
#'
#' @examples
deglr <- function(xtar,
                  ytar,
                  xext,
                  yext,
                  family,
                  nlambda = 50,
                  lambda.min.ratio = 1e-10,
                  fixed_lambda_seq = NULL,
                  return_only_at = NULL,
                  L = NULL,
                  standardize = TRUE,
                  thresh = 1e-14,
                  maxit = 1e5){

  d <- ncol(xtar)
  beta_index <- 1:(d + I(family != "cox"))
  if(!family %in% c("gaussian", "binomial", "cox")) {
    stop("'family' must be one of 'gaussian', 'binomial', or 'cox'")
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
  L = foo$L
  xext_star_scale = foo$xext_star_scale

  if(is.null(fixed_lambda_seq)) {
    # Run this solely to get the max lambda:
    mod <- glmnet(x = x_aug,
                  y = y_aug,
                  family = family,
                  alpha = 0.0,
                  nlambda = nlambda,
                  lambda.min.ratio = lambda.min.ratio,
                  #Note: lambda = NULL is not a typo
                  lambda = NULL,
                  penalty.factor = c(rep(0, d), rep(1, d + (family != "cox"))),
                  # Standardization taken care of in parent
                  standardize = FALSE)

    # Note that we scale up the max lambda by a factor of 100:
    max_log_lambda = log(100) + log(max(mod$lambda))
    lambda_seq <- exp(seq(from = max_log_lambda, to = max_log_lambda + log(lambda.min.ratio), length = nlambda))
  } else {
    lambda_seq <- fixed_lambda_seq
  }

  if(is.null(return_only_at)){
    which_lambda <- seq_along(lambda_seq);
  } else {
    lambda_seq <- rev(sort(unique(c(lambda_seq, return_only_at))))
    which_lambda <- match(return_only_at, lambda_seq)
  }

  mod2 <- glmnet(x = x_aug,
                 y = y_aug,
                 family = family,
                 alpha = 0.0,
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
}








