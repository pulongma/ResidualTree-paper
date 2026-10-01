run.cell = function(n, d, nugget, rho_min, case_config, results) {
  section = "N"
  k_eff = NA_integer_
  range_multiplier = NA_real_
  seed = case.seed(base.seed, n, d, nugget)
  ranges = geometric.ranges(d, s = rho_min)
  regime = regime.id(rho_min, nugget)

  exact_case = exact.gp.case(
    n, d, seed, ranges, nugget, known.sig2,
    cache.dir, case_config$exact_max_gb, rho_min,
    context = benchmark.context
  )
  if (exact_case$status != "ok") {
    warning(sprintf(
      "[%s] d=%d n=%d %s: exact case %s: %s",
      section, d, n, regime, exact_case$status, exact_case$message
    ))
    return(results)
  }

  X = exact_case$X
  y = exact_case$y
  loglik_exact = exact_case$exact_loglik
  theta = restree_theta(
    sig2 = known.sig2, range = ranges, nugget = nugget, nu = known.nu, form = "ARD"
  )
  reference_source = if (exact_case$from_cache) "cached" else "dense"
  cat(sprintf(
    "[%s] d=%2d n=%6d %-14s exact %.3f (%s data; likelihood %.3fs) rho %.4f..%.4f | %s\n",
    section, d, n, regime, loglik_exact, reference_source, exact_case$seconds,
    min(ranges), max(ranges), format(Sys.time(), "%H:%M:%S")
  ))

  loglik_package = NA_real_
  if (n <= crosscheck.n.max) {
    loglik_package = exact.gp.package(X, y, theta, ncores)
    allowed_difference = crosscheck.rtol * max(1, abs(loglik_exact))
    if (abs(loglik_package - loglik_exact) > allowed_difference) {
      stop(sprintf(
        "ResTree r=n likelihood %.6f != dense reference %.6f (d=%d n=%d %s)",
        loglik_package, loglik_exact, d, n, regime
      ))
    }
  }

  budgets = r.grid.for(case_config$r_grid, n)
  full_jobs = data.frame(
    method_id = "restgp_full", arm = "main", knob = budgets, pp_seed = NA_integer_
  )
  pp_jobs = expand.grid(
    method_id = "restgp_pp", arm = "main", knob = budgets,
    pp_seed = seq_len(case_config$pp_seeds), stringsAsFactors = FALSE
  )
  vecchia_jobs = data.frame(
    method_id = "ssgv", arm = "main", knob = m.grid.for(budgets, n),
    pp_seed = NA_integer_
  )
  jobs = rbind(full_jobs, pp_jobs, vecchia_jobs)

  for (job_index in seq_len(nrow(jobs))) {
    job = jobs[job_index, ]

    key = paste(
      section, paste0("nu=", known.nu), n, d, regime, job$method_id, job$knob, job$pp_seed,
      cut.method, knot.design, range.order, temper.alpha, resampling, restree.api,
      paste0("data_seed=", seed),
      sep = "|"
    )
    if (!is.null(results) && any(
      results$key == key & results$status == "ok" &
        results$timing_protocol == timing.protocol &
        results$timing_context == benchmark.context &
        (job$method_id != "restgp_pp" |
          results$nparticles == nparticles.for(tree.depth(n, job$knob))),
      na.rm = TRUE
    )) {
      next
    }

    gc(verbose = FALSE)
    evaluation = tryCatch(
      {
        value = evaluate(job, X, y, theta, ranges, nugget, n, seed)
        if (!is.finite(value$loglik)) {
          stop("Non-finite log likelihood")
        }
        value
      },
      error = function(error) error
    )

    failed = inherits(evaluation, "error")
    error_message = NA_character_
    if (failed) {
      error_message = conditionMessage(evaluation)
      warning(sprintf(
        "[%s] d=%d n=%d %s %s knob=%d: %s",
        section, d, n, regime, job$method_id, job$knob, error_message
      ))
      evaluation = list(
        loglik = NA_real_, depth = NA_integer_, nparticles = NA_integer_,
        n_knot = NA_real_, seconds_model = NA_real_, seconds = NA_real_,
        seconds_ordering = NA_real_, seconds_neighbors = NA_real_
      )
    }

    gap = evaluation$loglik - loglik_exact
    result_row = data.frame(
      key = key, section = section, nu = known.nu, d = d, n = n, seed = seed,
      k_eff = k_eff, range_multiplier = range_multiplier, nugget = nugget, regime = regime,
      range_order = range.order, rho_min = min(ranges), rho_max = max(ranges),
      method_id = job$method_id, method = unname(method.labels[job$method_id]),
      arm = job$arm, knob_name = unname(knob.names[job$method_id]),
      knob = job$knob, pp_seed = job$pp_seed,
      depth = evaluation$depth, nparticles = evaluation$nparticles,
      n_knot = evaluation$n_knot,
      loglik = evaluation$loglik, exact = loglik_exact, gap = gap,
      gap_per_obs = gap / n, abs_gap_per_obs = abs(gap) / n,
      seconds_model = evaluation$seconds_model, seconds = evaluation$seconds,
      seconds_total = evaluation$seconds_model + evaluation$seconds,
      seconds_ordering = evaluation$seconds_ordering,
      seconds_neighbors = evaluation$seconds_neighbors,
      timing_protocol = timing.protocol, timing_context = benchmark.context,
      timing_host = unname(Sys.info()["nodename"]),
      exact_seconds = exact_case$seconds, exact_package = loglik_package,
      exact_seconds_covariance = exact_case$seconds_covariance,
      exact_seconds_cholesky = exact_case$seconds_cholesky,
      exact_seconds_solve = exact_case$seconds_solve,
      exact_simulation_seconds = exact_case$seconds_simulation,
      exact_legacy_seconds = if (is.null(exact_case$legacy_generation_seconds)) NA_real_ else exact_case$legacy_generation_seconds,
      exact_timing_protocol = exact_case$timing_protocol,
      exact_timing_context = exact_case$timing_context,
      exact_timing_host = exact_case$timing_host,
      exact_timing_utc = exact_case$timing_utc,
      exact_timing_reused = exact_case$timing_reused,
      cut_method = cut.method, design = knot.design,
      prior_rho = prior.rho, prior_beta = prior.beta,
      resampling = resampling, temper_alpha = temper.alpha, ncores = ncores,
      package_version = as.character(packageVersion("ResTree")),
      package_api = restree.api, data_file = exact_case$file,
      status = if (failed) "failed" else "ok", message = error_message,
      time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      stringsAsFactors = FALSE
    )

    if (!is.null(results)) {
      results = results[results$key != key, ]
    }
    results = rbind(results, result_row)
    write.summary(results, section)
    cat(sprintf(
      "[%s] d=%2d n=%6d %-14s %-18s %s=%4d  gap/n %+.3e  (%.3fs)\n",
      section, d, n, regime, job$method_id, result_row$knob_name,
      job$knob, result_row$gap_per_obs, result_row$seconds_total
    ))
  }
  results
}
