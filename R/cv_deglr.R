#' Use cross-validation to select appropriate value of the tuning parameter in data-enriched generalized linear regression
#'
#' @param xtar Predictor matrix from the target population
#' @param ytar Outcome vector from the target population
#' @param xext Predictor matrix from the external population
#' @param yext Outcome vector from the external population
#' @param family One of 'gaussian', 'binomial', or 'cox'
#' @param alpha See `glmnet::cv.glmnet()`
#' @param nfolds See `glmnet::cv.glmnet()`
#' @param ncvreps Number of partitions to create. Choosing > 1 improves
#'   stability of selection but takes longer
#' @param nlambda See `glmnet::glmnet()` but note that in contrast to glmnet() the
#'   sequence of lambdas will always be exactly this long
#' @param lambda.min.ratio See `glmnet::glmnet()`
#' @param lambda User can provide a fixed sequence of lambdas instead
#'   of an automaticaly constructed sequence. If provided, this will override
#'   `nlambda` and `lambda.min.ratio`
#' @param L Square matrix that scales the bias parameter gamma. If `NULL`, then
#'   the `make_L()` function will be called to set `L` as the cholesky
#'   decomposition of `crossprod(xtar)`
#' @param standardize If TRUE, then `xext` will be standardized but only in the
#'   part of the augmented data matrix corresponding to gamma.
#' @param thresh See glmnet::glmnet()
#' @param maxit See glmnet::glmnet()
#'
#' @return A named list
#'
#' @export
#'
#' @examples

cv_deglr <- function(xtar,
                     ytar,
                     xext,
                     yext,
                     family,
                     alpha = 0,
                     nfolds = 5,
                     ncvreps = 2,
                     nlambda = 20,
                     lambda.min.ratio = ifelse(nrow(xtar) < d/2, 1e-3, 1e-6),
                     lambda = NULL,
                     L = NULL,
                     standardize = TRUE,
                     thresh = 1e-14,
                     maxit = 1e5) {

  n_tar <- nrow(xtar)
  n_ext <- nrow(xext)
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

  if(is.null(lambda)) {
    mod_full <- deglr(xtar = xtar,
                      ytar = ytar,
                      xext = xext,
                      yext = yext,
                      family = family,
                      alpha = alpha,
                      nlambda = nlambda,
                      lambda.min.ratio = lambda.min.ratio,
                      lambda = NULL,
                      return_what = "lambda_seq",
                      L = L,
                      standardize = standardize)
    lambda_seq <- mod_full$lambda;
  } else {
    lambda_seq <- lambda;
  }

  n_tar.k <- floor(n_tar / nfolds) * nfolds
  neg_log_likelihood <- numeric(length(lambda_seq))
  for(reps in 1:ncvreps) {
    i.e <- sample(x = n_ext, size = round(n_ext * (nfolds - 1) / nfolds ), replace = FALSE)
    ind <- matrix(sample(x = n_tar, size = n_tar.k, replace = FALSE), ncol = nfolds, nrow = n_tar.k/nfolds)
    for(k in 1:nfolds) {


      mod_k <-  deglr(xtar = xtar[-ind[,k],,drop = FALSE],
                      ytar = ytar[-ind[,k]],
                      xext = xext[i.e,,drop = FALSE],
                      yext = yext[i.e],
                      family = family,
                      alpha = alpha,
                      lambda = lambda_seq,
                      L = L,
                      standardize = standardize,
                      thresh = thresh,
                      maxit = maxit)

      ytar.k <- ytar[ind[,k]]
      if(family != "cox") {
        xtar.beta.k_hat <- cbind(1, xtar[ind[,k],,drop = FALSE]) %*% as.matrix(mod_k$beta_hat)
      } else {
        xtar.beta.k_hat <- xtar[ind[,k],,drop = FALSE] %*% as.matrix(mod_k$beta_hat)
      }

      if(family == "gaussian") {
        neg_log_likelihood <-
          neg_log_likelihood +
          colSums((ytar.k - xtar.beta.k_hat)^2)  / ncvreps;
      } else if(family == "binomial") {
        neg_log_likelihood <-
          neg_log_likelihood -
          drop(t(ytar.k) %*% (xtar.beta.k_hat)) / ncvreps +
          colSums(log1plex(xtar.beta.k_hat)) / ncvreps;
      } else if(family == "cox") {
        neg_log_likelihood <-
          neg_log_likelihood -
          coxph_ploglik(ytar.k, xtar.beta.k_hat);
      }
    }

  }

  lambda_opt <- lambda_seq[which.min(neg_log_likelihood)]

  mod_full <- deglr(xtar = xtar,
                    ytar = ytar,
                    xext = xext,
                    yext = yext,
                    family = family,
                    alpha = alpha,
                    lambda = lambda_seq,
                    return_what = lambda_opt,
                    L = L,
                    standardize = standardize,
                    thresh = thresh,
                    maxit = maxit)

  return(list(neg_log_likelihood = neg_log_likelihood,
              lambda_opt = lambda_opt,
              beta_hat =  mod_full$beta_hat,
              gamma_hat = mod_full$gamma_hat,
              lambda_seq = lambda_seq))

}

