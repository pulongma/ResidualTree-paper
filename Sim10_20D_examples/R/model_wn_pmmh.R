# Standalone model file: needs only ResTree and GpGp.
# Each job fits its own training-only initializer; no other job is needed.
fit_wn_pmmh <- function(X, y, X_test, settings, ncores = 1L) {
  set.seed(settings$seed + 10L)
  started <- proc.time()[3]

  # 1. Initialize covariance parameters with a training-only Vecchia fit.
  dimension <- ncol(X)
  intercept <- matrix(1, nrow(X), 1)
  initial_fit <- GpGp::fit_model(
    y, X,
    X = intercept, covfun_name = "matern25_scaledim",
    m_seq = c(10L, settings$m), max_iter = settings$vecchia_maxit, silent = TRUE
  )
  initial_covparms <- initial_fit$covparms

  theta <- ResTree::restree_theta(
    sig2 = initial_covparms[1], range = initial_covparms[1L + seq_len(dimension)],
    nugget = initial_covparms[dimension + 2],
    nu = settings$nu, form = "ARD", family = "matern"
  )
  initialization_seconds <- unname(proc.time()[3] - started)

  # 2. Fit the model. The timer includes the Vecchia initialization.
  model <- ResTree::restree_model(
    X, y,
    depth = settings$depth, r = settings$r, leaf_model = "WN",
    design = "maximin", cut_method = "uniform",
    prior_rho = settings$prior_rho, prior_beta = settings$prior_beta, seed = settings$seed + 20L
  )
  # PMMH uses proper parameter priors. These match current package defaults:
  # sig2 ~ inverse-gamma(2, var(y)); ranges and nugget ~ half-Cauchy(1).
  # var(y) uses training responses only. All settings are saved with the chain.
  prior <- list(sig2 = c(shape = 2, rate = stats::var(y)), range = 1, nugget = 1)
  fit <- ResTree::restree_fit(model, theta,
    method = "pmmh",
    nparticles = settings$particles_per_depth * settings$depth, n_iter = settings$pmmh_iter,
    seed = settings$seed + 20L, ncores = ncores,
    resampling = "stratified", temper_alpha = 1,
    control = list(
      burnin = settings$pmmh_burnin,
      prediction = "mcmc", point_estimate = "median",
      prior = prior, prop_sd = settings$pmmh_prop_sd,
      sig2_moves = 0L, # cheap variance-only moves require deterministic cuts
      adapt = settings$pmmh_adapt, target_accept = settings$pmmh_target_accept,
      sd_reps = settings$pmmh_sd_reps
    )
  )
  elapsed_fit <- unname(proc.time()[3] - started)
  # The package times its extra evidence-noise checks separately. Exclude
  # those diagnostics from fitting time; the initializer remains included.
  evidence_check_seconds <- unname(fit$fit$sd_reps_time_s)
  fit_seconds <- elapsed_fit - evidence_check_seconds

  # 3. Predict from the saved PMMH trees; this does not run SMC again.
  # ResTree averages every post-burn-in state, including rejection repeats.
  # Batch test points because component matrices grow with the number of test points and trees.
  started <- proc.time()[3]
  prediction_mean <- prediction_variance <- numeric(nrow(X_test))
  for (first in seq.int(1L, nrow(X_test), by = settings$pmmh_prediction_batch)) {
    rows <- seq.int(first, min(first + settings$pmmh_prediction_batch - 1L, nrow(X_test)))
    prediction <- ResTree::restree_predict(fit, X_test[rows, , drop = FALSE],
      prediction = "mcmc", ncores = ncores, keep_history = FALSE
    )
    prediction_mean[rows] <- prediction$mean
    prediction_variance[rows] <- prediction$var # within-tree + between-tree/iteration variance
    rm(prediction)
    invisible(gc()) # release component matrices before the next batch
    cat("PMMH prediction:", max(rows), "of", nrow(X_test), "test points\n")
    flush.console()
  }
  prediction_seconds <- unname(proc.time()[3] - started)

  # 4. Return predictions, timing and parameter estimates.
  diagnostics <- ResTree::restree_diagnostics(fit)
  diagnostics$ESS_final <- utils::tail(fit$ESS, 1L) # particle ESS, not chain ESS
  diagnostics$evidence_check_seconds <- evidence_check_seconds
  mcmc <- list(
    chain = fit$chain, burnin = fit$burnin,
    prediction_index = fit$prediction_index, accepted = fit$accepted,
    prediction_method = fit$prediction_method,
    prediction_source = "retained_smc_ensembles",
    prediction_batch = settings$pmmh_prediction_batch,
    state = fit$mcmc$trees$state, state_seed = fit$mcmc$trees$seed,
    point_estimate = fit$point_estimate,
    logZ_trace = fit$logZ_trace, prop_scale_trace = fit$prop_scale_trace,
    n_requested = settings$pmmh_iter, interrupted = fit$diagnostics$interrupted
  )
  fitted_theta <- list(
    sig2 = fit$theta$sig2, range = fit$theta$range,
    nugget = fit$theta$nugget, nu = fit$theta$nu,
    family = fit$theta$family, form = fit$theta$form,
    range_original = fit$theta$range
  )
  list(
    mean = prediction_mean, variance = prediction_variance,
    fit_seconds = fit_seconds, prediction_seconds = prediction_seconds,
    total_seconds = fit_seconds + prediction_seconds,
    initialization_seconds = initialization_seconds,
    particles = as.integer(fit$nparticles), convergence = NA_integer_,
    estimator = "PMMH (posterior predictive average)",
    diagnostics = list(
      fit = diagnostics, theta = fitted_theta,
      initial_covparms = initial_covparms,
      initial_covfun = initial_fit$covfun_name,
      actual_method = fit$method,
      prior = fit$fit$prior, mcmc = mcmc,
      design = "maximin", cut_method = "uniform"
    )
  )
}
