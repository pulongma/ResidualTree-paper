# Standalone Vecchia GP: needs GpGp. Test responses are never supplied.
# Raw: fit once. Scaled: refit twice after dividing inputs by fitted ARD ranges.
# The final fit, prediction and simulations all use the same input scaling.
fit_vecchia <- function(X, y, X_test, settings, ncores = 1L) {
  s <- settings
  set.seed(s$seed + 10L)
  d <- ncol(X)
  range_columns <- 1L + seq_len(d)
  intercept <- matrix(1, nrow(X), 1)
  started <- proc.time()[3]
  fit <- GpGp::fit_model(
    y, X, X = intercept, covfun_name = "matern25_scaledim",
    m_seq = c(10L, s$m), max_iter = s$vecchia_maxit, silent = TRUE
  )
  range_divisor <- rep(1, d)  # cumulative division of the inputs by fitted ranges
  for (step in seq_len(s$vecchia_scaling_rounds)) {
    range_divisor <- range_divisor * fit$covparms[range_columns]
    fit <- GpGp::fit_model(
      y, sweep(X, 2, range_divisor, "/"), X = intercept, covfun_name = "matern25_scaledim",
      m_seq = c(10L, s$m), max_iter = s$vecchia_maxit, silent = TRUE,
      start_parms = c(fit$covparms[1], rep(1, d), fit$covparms[d + 2])
    )
  }
  covparms <- fit$covparms
  covparms[range_columns] <- covparms[range_columns] * range_divisor  # ranges on the original inputs
  fit_seconds <- unname(proc.time()[3] - started)

  # Conditional response simulations include the fitted observation nugget.
  started <- proc.time()[3]
  X_test_model <- X_test  # raw jobs predict on the original held-out inputs
  if (s$vecchia_scaling_rounds > 0L) {
    X_test_model <- sweep(X_test, 2, range_divisor, "/")
  }
  intercept_test <- matrix(1, nrow(X_test), 1)
  mu <- GpGp::predictions(fit = fit, locs_pred = X_test_model,
    X_pred = intercept_test, m = s$m)
  draws <- GpGp::cond_sim(fit = fit, locs_pred = X_test_model,
    X_pred = intercept_test, m = s$m, nsims = s$vecchia_simulations)
  variance <- apply(draws, 1, stats::var)
  prediction_seconds <- unname(proc.time()[3] - started)

  list(mean = as.numeric(mu), variance = variance,
    fit_seconds = fit_seconds, prediction_seconds = prediction_seconds,
    total_seconds = fit_seconds + prediction_seconds,
    initialization_seconds = 0, particles = NA_integer_, convergence = NA_integer_,
    estimator = if (s$vecchia_scaling_rounds > 0) "Scaled Vecchia maximum likelihood" else "Vecchia maximum likelihood",
    diagnostics = list(covparms = covparms, covparms_scaled = fit$covparms,
                       input_scale = range_divisor, scaling_rounds = s$vecchia_scaling_rounds,
                       intercept = fit$betahat, simulations = ncol(draws)))
}
