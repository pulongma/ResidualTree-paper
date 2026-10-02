# Standalone model file: needs only BayesTree. BayesTree uses a constant noise SD.
fit_bart <- function(X, y, X_test, settings, ncores = 1L) {
  set.seed(settings$seed + 30L)
  started <- proc.time()[3]

  # Fitting and test prediction happen in the same MCMC call.
  fit <- BayesTree::bart(
    x.train = X, y.train = as.numeric(y), x.test = X_test,
    ntree = settings$bart_trees, nskip = settings$bart_burn, ndpost = settings$bart_draws,
    keepevery = 1L, keeptrainfits = FALSE,
    usequants = FALSE, numcut = 100L, printevery = 1000L, verbose = FALSE
  )
  prediction_mean <- as.numeric(fit$yhat.test.mean)
  latent_variance <- apply(fit$yhat.test, 2, stats::var)
  # Add residual variance: yhat.test contains f(x) draws, not new noisy y draws.
  prediction_variance <- latent_variance + mean(fit$sigma^2)
  latent_interval <- apply(fit$yhat.test, 2, stats::quantile, c(0.025, 0.975))
  total_seconds <- unname(proc.time()[3] - started)

  list(
    mean = prediction_mean, variance = prediction_variance,
    latent_lower = latent_interval[1, ], latent_upper = latent_interval[2, ],
    fit_seconds = NA_real_, prediction_seconds = NA_real_,
    total_seconds = total_seconds, initialization_seconds = 0,
    particles = NA_integer_, convergence = NA_integer_,
    estimator = "BART posterior mean",
    diagnostics = list(
      sigma_draws = fit$sigma,
      latent_variance = latent_variance, noise_variance = mean(fit$sigma^2),
      retained_draws = nrow(fit$yhat.test), burn = settings$bart_burn, trees = settings$bart_trees
    )
  )
}
