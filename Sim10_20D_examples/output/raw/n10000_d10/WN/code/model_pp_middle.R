# Standalone model file: needs only ResTree and GpGp.
# PP always uses deterministic middle cuts, for raw and scaled inputs.
# Each job fits its own training-only initializer; no other job is needed.
fit_pp_middle <- function(X, y, X_test, settings, ncores = 1L) {
  s <- settings
  set.seed(s$seed + 10L)
  started <- proc.time()[3]

  # 1. Fit the Vecchia initializer using training data only.
  # Raw: one fit. Scaled: two additional fits after dividing by fitted ranges.
  d <- ncol(X)
  range_columns <- 1L + seq_len(d)
  intercept <- matrix(1, nrow(X), 1)
  initial <- GpGp::fit_model(
    y, X, X = intercept, covfun_name = "matern25_scaledim",
    m_seq = c(10L, s$m), max_iter = s$vecchia_maxit, silent = TRUE
  )
  range_divisor <- rep(1, d)
  for (step in seq_len(s$vecchia_scaling_rounds)) {
    range_divisor <- range_divisor * initial$covparms[range_columns]
    initial <- GpGp::fit_model(
      y, sweep(X, 2, range_divisor, "/"), X = intercept, covfun_name = "matern25_scaledim",
      m_seq = c(10L, s$m), max_iter = s$vecchia_maxit, silent = TRUE,
      start_parms = c(initial$covparms[1], rep(1, d), initial$covparms[d + 2])
    )
  }
  # Convert the final initializer ranges back to the original input units.
  initial_covparms <- initial$covparms
  initial_covparms[range_columns] <- initial_covparms[range_columns] * range_divisor
  initial_ranges <- initial_covparms[range_columns]

  # Main (raw) study: use the original training and test inputs directly.
  input_scale <- rep(1, d)
  X_model <- X
  X_test_model <- X_test
  # Separate scaled study: use the same training-derived scale for both designs.
  if (s$scale_inputs) {
    input_scale <- min(initial_ranges) / initial_ranges
    X_model <- sweep(X, 2, input_scale, "*")
    X_test_model <- sweep(X_test, 2, input_scale, "*")
  }
  theta <- ResTree::restree_theta(
    sig2 = initial_covparms[1], range = initial_ranges * input_scale,
    nugget = initial_covparms[d + 2],
    nu = s$nu, form = "ARD", family = "matern"
  )
  initialization_seconds <- unname(proc.time()[3] - started)

  # 2. Construct and fit the residual-tree model.
  model <- ResTree::restree_model(
    X_model, y, depth = s$depth, r = s$r, leaf_model = "PP",
    design = "maximin", cut_method = "middle",
    prior_rho = s$prior_rho, prior_beta = s$prior_beta, seed = s$seed + 20L
  )
  # Empirical Bayes maximizes tree-integrated evidence, with no parameter
  # prior. The package uses this same particle budget for EB and final SMC.
  fit <- ResTree::restree_fit(model, theta, method = "ebayes",
    nparticles = s$particles_per_depth * s$depth,
    seed = s$seed + 20L, ncores = ncores,
    resampling = "stratified", temper_alpha = 1,
    control = list(maxit = s$eb_maxit, prior = "none", sd_reps = 0L))
  fit_seconds <- unname(proc.time()[3] - started)

  # 3. Predict noisy responses at held-out inputs, without using test responses.
  started <- proc.time()[3]
  prediction <- ResTree::restree_predict(fit, X_test_model,
                                         ncores = ncores, keep_history = FALSE)
  mu <- as.numeric(prediction$mean)
  variance <- as.numeric(prediction$var)
  prediction_seconds <- unname(proc.time()[3] - started)

  # 4. Diagnostics and saving are outside the fit/prediction timers.
  diagnostics <- ResTree::restree_diagnostics(fit)
  fitted_theta <- list(
    sig2 = fit$theta$sig2, range = fit$theta$range,
    nugget = fit$theta$nugget, nu = fit$theta$nu,
    family = fit$theta$family, form = fit$theta$form,
    range_original = fit$theta$range / input_scale
  )
  list(mean = mu, variance = variance, prediction = prediction,
    fit_seconds = fit_seconds, prediction_seconds = prediction_seconds,
    total_seconds = fit_seconds + prediction_seconds,
    initialization_seconds = initialization_seconds,
    particles = as.integer(fit$nparticles), convergence = as.integer(diagnostics$convergence),
    estimator = "Empirical Bayes (maximum evidence)",
    diagnostics = list(fit = diagnostics, theta = fitted_theta,
      input_scale = input_scale, scale_inputs = s$scale_inputs,
      initial_covparms = initial_covparms,
      vecchia_scaling_rounds = s$vecchia_scaling_rounds, actual_method = fit$method,
      design = "maximin", cut_method = "middle"))
}
