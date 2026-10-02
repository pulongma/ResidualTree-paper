# Standalone model file: needs only ResTree and GpGp.
# Each job fits its own training-only initializer; no other job is needed.
fit_wn_pmmh <- function(X, y, X_test, settings, ncores = 1L) {
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
    X_model, y, depth = s$depth, r = s$r, leaf_model = "WN",
    design = "maximin", cut_method = "uniform",
    prior_rho = s$prior_rho, prior_beta = s$prior_beta, seed = s$seed + 20L
  )
  # PMMH uses proper parameter priors. These match current package defaults:
  # sig2 ~ inverse-gamma(2, var(y)); ranges and nugget ~ half-Cauchy(1).
  # var(y) uses training responses only. All settings are saved with the chain.
  prior <- list(sig2 = c(shape = 2, rate = stats::var(y)), range = 1, nugget = 1)
  fit <- ResTree::restree_fit(model, theta, method = "pmmh",
    nparticles = s$particles_per_depth * s$depth, n_iter = s$pmmh_iter,
    seed = s$seed + 20L, ncores = ncores,
    resampling = "stratified", temper_alpha = 1,
    control = list(burnin = s$pmmh_burnin,
      prediction = "mcmc", point_estimate = "median",
      prior = prior, prop_sd = s$pmmh_prop_sd,
      sig2_moves = 0L,  # cheap variance-only moves require deterministic cuts
      adapt = s$pmmh_adapt, target_accept = s$pmmh_target_accept,
      sd_reps = s$pmmh_sd_reps))
  elapsed_fit <- unname(proc.time()[3] - started)
  # The package times its extra evidence-noise checks separately. Exclude
  # those diagnostics from fitting time; the initializer remains included.
  evidence_check_seconds <- unname(fit$fit$sd_reps_time_s)
  fit_seconds <- elapsed_fit - evidence_check_seconds

  # 3. Predict from the saved PMMH trees; this does not run SMC again.
  # ResTree averages every post-burn-in state, including rejection repeats.
  # Batch test points because component matrices grow with the number of test points and trees.
  started <- proc.time()[3]
  mu <- variance <- numeric(nrow(X_test))
  for (first in seq.int(1L, nrow(X_test), by = s$pmmh_prediction_batch)) {
    rows <- seq.int(first, min(first + s$pmmh_prediction_batch - 1L, nrow(X_test)))
    prediction <- ResTree::restree_predict(fit, X_test_model[rows, , drop = FALSE],
      prediction = "mcmc", ncores = ncores, keep_history = FALSE)
    mu[rows] <- prediction$mean
    variance[rows] <- prediction$var  # within-tree + between-tree/iteration variance
    rm(prediction)
    invisible(gc())  # release component matrices before the next batch
    cat("PMMH prediction:", max(rows), "of", nrow(X_test), "test points\n")
    flush.console()
  }
  prediction_seconds <- unname(proc.time()[3] - started)

  # 4. Diagnostics and saving are outside the fit/prediction timers.
  diagnostics <- ResTree::restree_diagnostics(fit)
  diagnostics$ESS_final <- utils::tail(fit$ESS, 1L)  # particle ESS, not chain ESS
  diagnostics$evidence_check_seconds <- evidence_check_seconds
  mcmc <- list(chain = fit$chain, burnin = fit$burnin,
    prediction_index = fit$prediction_index, accepted = fit$accepted,
    prediction_method = fit$prediction_method,
    prediction_source = "retained_smc_ensembles",
    prediction_batch = s$pmmh_prediction_batch,
    state = fit$mcmc$trees$state, state_seed = fit$mcmc$trees$seed,
    point_estimate = fit$point_estimate,
    logZ_trace = fit$logZ_trace, prop_scale_trace = fit$prop_scale_trace,
    n_requested = s$pmmh_iter, interrupted = fit$diagnostics$interrupted)
  fitted_theta <- list(
    sig2 = fit$theta$sig2, range = fit$theta$range,
    nugget = fit$theta$nugget, nu = fit$theta$nu,
    family = fit$theta$family, form = fit$theta$form,
    range_original = fit$theta$range / input_scale
  )
  list(mean = mu, variance = variance,
    fit_seconds = fit_seconds, prediction_seconds = prediction_seconds,
    total_seconds = fit_seconds + prediction_seconds,
    initialization_seconds = initialization_seconds,
    particles = as.integer(fit$nparticles), convergence = NA_integer_,
    estimator = "PMMH (posterior predictive average)",
    diagnostics = list(fit = diagnostics, theta = fitted_theta,
      input_scale = input_scale, scale_inputs = s$scale_inputs,
      initial_covparms = initial_covparms,
      vecchia_scaling_rounds = s$vecchia_scaling_rounds,
      initial_covfun = initial$covfun_name,
      actual_method = fit$method,
      prior = fit$fit$prior, mcmc = mcmc,
      design = "maximin", cut_method = "uniform"))
}
