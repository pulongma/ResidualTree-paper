# Standalone model file: needs only ResTree and GpGp.
# PP uses deterministic middle cuts on the original inputs.
# Each job fits its own training-only initializer; no other job is needed.
fit_pp_middle <- function(X, y, X_test, settings, ncores = 1L) {
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
    depth = settings$depth, r = settings$r, leaf_model = "PP",
    design = "maximin", cut_method = "middle",
    prior_rho = settings$prior_rho, prior_beta = settings$prior_beta, seed = settings$seed + 20L
  )
  # Empirical Bayes uses the same particle budget for fitting and final SMC.
  fit <- ResTree::restree_fit(model, theta,
    method = "ebayes",
    nparticles = settings$particles_per_depth * settings$depth,
    seed = settings$seed + 20L, ncores = ncores,
    resampling = "stratified", temper_alpha = 1,
    control = list(maxit = settings$eb_maxit, prior = "none", sd_reps = 0L)
  )
  fit_seconds <- unname(proc.time()[3] - started)

  # 3. Predict the held-out responses; their observed values are not used.
  started <- proc.time()[3]
  prediction <- ResTree::restree_predict(fit, X_test,
    ncores = ncores, keep_history = FALSE
  )
  prediction_mean <- as.numeric(prediction$mean)
  prediction_variance <- as.numeric(prediction$var)
  prediction_seconds <- unname(proc.time()[3] - started)

  # 4. Return predictions, timing and parameter estimates.
  diagnostics <- ResTree::restree_diagnostics(fit)
  fitted_theta <- list(
    sig2 = fit$theta$sig2, range = fit$theta$range,
    nugget = fit$theta$nugget, nu = fit$theta$nu,
    family = fit$theta$family, form = fit$theta$form,
    range_original = fit$theta$range
  )
  list(
    mean = prediction_mean, variance = prediction_variance, prediction = prediction,
    fit_seconds = fit_seconds, prediction_seconds = prediction_seconds,
    total_seconds = fit_seconds + prediction_seconds,
    initialization_seconds = initialization_seconds,
    particles = as.integer(fit$nparticles), convergence = as.integer(diagnostics$convergence),
    estimator = "Empirical Bayes (maximum evidence)",
    diagnostics = list(
      fit = diagnostics, theta = fitted_theta,
      initial_covparms = initial_covparms,
      actual_method = fit$method,
      design = "maximin", cut_method = "middle"
    )
  )
}
