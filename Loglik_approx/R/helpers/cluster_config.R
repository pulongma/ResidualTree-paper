cluster.config = function(sample_sizes, dimensions, minimum_ranges, nugget_values,
                          knot_counts, pp_seed_count, accuracy_tolerances, output_suffix) {
  ncores = 10L
  local_cores = Sys.getenv("RESTREE_NCORES", "")
  slurm_cores = Sys.getenv("SLURM_CPUS_PER_TASK", "")
  if (nzchar(local_cores)) {
    ncores = as.integer(local_cores)
  }
  if (nzchar(slurm_cores)) {
    ncores = as.integer(slurm_cores)
  }
  exact_memory_limit_gb = as.numeric(Sys.getenv("RESTREE_HD_EXACT_MAX_GB", "50"))
  Sys.setenv(OMP_NUM_THREADS = ncores)

  regimes = expand.grid(
    d = dimensions,
    rho_min = minimum_ranges,
    nugget = nugget_values
  )

  config = list(
    mode = "cluster",
    base_seed = 2023L,
    output_suffix = output_suffix,
    min_package_version = "0.1.8",
    ncores = ncores,
    tol = accuracy_tolerances,
    N = list(
      n = sample_sizes,
      d = dimensions,
      regimes = regimes,
      r_grid = knot_counts,
      pp_seeds = pp_seed_count,
      exact_max_gb = exact_memory_limit_gb
    )
  )

  config
}
