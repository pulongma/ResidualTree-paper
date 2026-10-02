# Study-wide settings; method-specific choices are in the model files.
# All methods use the original inputs.
study_settings <- function(n, d) {
  stopifnot(n %in% c(5000L, 10000L, 50000L), d %in% c(10L, 20L))

  list(
    # Simulation design and reproducible seeds.
    study = "Multiscale-2026-09-23-100-depth", variant = "raw",
    n = as.integer(n), n_test = as.integer(n / 5L), d = as.integer(d),
    replicate = 1L, seed = as.integer(20260923L + n + 1000L * d),

    # GP covariance, tree and optimization settings.
    r = 60L, # ResTree knot budget / Full leaf size
    m = 60L, # Vecchia neighbors
    depth = as.integer(ceiling(log2(n / 60))),
    particles_per_depth = 100L, # study budget; independent of package defaults
    nu = 2.5,
    eb_maxit = 100L,
    full_starts = 2L,
    vecchia_maxit = 100L,
    vecchia_simulations = 200L,

    # PMMH: retain every post-burn-in state for prediction.
    pmmh_iter = 3000L, # total, including burn-in
    pmmh_burnin = 1000L,
    # Batching splits test inputs; it does not thin the chain.
    pmmh_prediction_batch = 32L,
    pmmh_sd_reps = 10L,
    pmmh_prop_sd = 0.15,
    pmmh_adapt = "burnin", pmmh_target_accept = 0.15,
    prior_rho = 0.05, prior_beta = 0,

    # BART: burn-in followed by retained posterior draws.
    bart_trees = 200L,
    bart_burn = 3000L,
    bart_draws = 7000L
  )
}

METHODS <- c("Full", "PP_middle", "WN", "Vecchia", "BART", "WN_PMMH")
