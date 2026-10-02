# Main cluster jobs use raw inputs. Scaled jobs are submitted separately.
# The WN PMMH batch job is provided for n = 10,000.
# variant = "raw":    GpGp orders and picks neighbours on the raw inputs and
#                     the trees place maximin knots on the raw inputs
#                     (the Michalewicz study's code path; PP uses middle cuts).
# variant = "scaled": scaled Vecchia (refits on range-scaled inputs) and trees
#                     on inputs multiplied by min(range)/range_j.
# Both variants use middle cuts for PP and uniform cuts for WN (EB and PMMH).
VARIANTS <- c("raw", "scaled")
study_settings <- function(n, d, variant, smoke = FALSE) {
  stopifnot(d %in% c(10L, 20L), n >= 120L, variant %in% VARIANTS)
  if (!smoke) stopifnot(n %in% c(5000L, 10000L, 50000L))

  list(
    study = "Multiscale-2026-09-23-100-depth", smoke = smoke, variant = variant,
    n = as.integer(n), n_test = as.integer(n / 5L), d = as.integer(d),
    replicate = 1L, seed = as.integer(20260923L + n + 1000L * d),
    r = 60L,                          # ResTree knot budget / Full leaf size
    m = 60L,                          # Vecchia neighbors
    depth = as.integer(ceiling(log2(n / 60))),
    particles_per_depth = 100L,        # explicit PP/WN budget, including final SMC
    nu = 2.5,
    eb_maxit = if (smoke) 3L else 100L,
    full_starts = if (smoke) 1L else 2L,
    vecchia_maxit = if (smoke) 3L else 100L,
    vecchia_scaling_rounds = if (variant == "scaled") 2L else 0L,
    scale_inputs = variant == "scaled",
    vecchia_simulations = 200L,
    pmmh_iter = if (smoke) 20L else 3000L,  # total, including burn-in
    pmmh_burnin = if (smoke) 10L else 1000L,
    # Every retained state is averaged; batching only splits test points.
    pmmh_prediction_batch = if (smoke) 7L else 32L,
    pmmh_sd_reps = if (smoke) 2L else 10L,
    pmmh_prop_sd = 0.15,
    pmmh_adapt = "burnin", pmmh_target_accept = 0.15,
    prior_rho = 0.05, prior_beta = 0,
    bart_trees = 200L,
    bart_burn = if (smoke) 10L else 3000L,
    bart_draws = if (smoke) 20L else 7000L
  )
}

METHODS <- c("Full", "PP_middle", "WN", "WN_middle", "Vecchia", "BART", "WN_PMMH")
