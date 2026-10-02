# Standalone Vecchia GP: needs GpGp. Test responses are never supplied.
# Fit and predict on the original inputs. The GpGp covariance name
# "matern25_scaledim" denotes an ARD kernel, not an input transformation.
fit_vecchia <- function(X, y, X_test, settings, ncores = 1L) {
  set.seed(settings$seed + 10L)
  intercept <- matrix(1, nrow(X), 1)
  started <- proc.time()[3]
  fit <- GpGp::fit_model(
    y, X,
    X = intercept, covfun_name = "matern25_scaledim",
    m_seq = c(10L, settings$m), max_iter = settings$vecchia_maxit, silent = TRUE
  )
  fit_seconds <- unname(proc.time()[3] - started)

  # Conditional response simulations include the fitted observation nugget.
  started <- proc.time()[3]
  intercept_test <- matrix(1, nrow(X_test), 1)
  prediction_mean <- GpGp::predictions(
    fit = fit, locs_pred = X_test,
    X_pred = intercept_test, m = settings$m
  )
  draws <- GpGp::cond_sim(
    fit = fit, locs_pred = X_test,
    X_pred = intercept_test, m = settings$m, nsims = settings$vecchia_simulations
  )
  prediction_variance <- apply(draws, 1, stats::var)
  prediction_seconds <- unname(proc.time()[3] - started)

  list(
    mean = as.numeric(prediction_mean), variance = prediction_variance,
    fit_seconds = fit_seconds, prediction_seconds = prediction_seconds,
    total_seconds = fit_seconds + prediction_seconds,
    initialization_seconds = 0, particles = NA_integer_, convergence = NA_integer_,
    estimator = "Vecchia maximum likelihood",
    diagnostics = list(
      covparms = fit$covparms,
      intercept = fit$betahat, simulations = ncol(draws)
    )
  )
}
