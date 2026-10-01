required.knob = function(results, tols) {
  results = results[results$status == "ok", ]
  group_columns = c(
    "section", "nu", "d", "n", "k_eff", "range_multiplier", "nugget", "regime", "method_id", "method", "pp_seed"
  )

  groups = results[group_columns]
  groups$pp_seed[is.na(groups$pp_seed)] = 0L

  groups$k_eff[is.na(groups$k_eff)] = 0L
  groups$range_multiplier[is.na(groups$range_multiplier)] = 0
  grouped_results = split(results, groups, drop = TRUE)

  summaries = lapply(grouped_results, function(method_results) {
    method_results = method_results[order(method_results$knob), ]
    best_row = which.min(method_results$abs_gap_per_obs)

    tolerance_rows = lapply(tols, function(tol) {
      first_hit = which(method_results$abs_gap_per_obs <= tol)[1L]
      cbind(
        method_results[1L, group_columns],
        tol = tol,
        knob_name = unname(knob.names[method_results$method_id[1L]]),
        knob_required = method_results$knob[first_hit],
        seconds_at_required = method_results$seconds_total[first_hit],
        knob_max_tried = max(method_results$knob),
        best_abs_gap_per_obs = method_results$abs_gap_per_obs[best_row],
        best_gap_per_obs = method_results$gap_per_obs[best_row],
        best_knob = method_results$knob[best_row]
      )
    })
    do.call(rbind, tolerance_rows)
  })

  summary = do.call(rbind, summaries)
  rownames(summary) = NULL
  summary[order(
    summary$section, summary$d, summary$n, summary$k_eff,
    summary$nugget, summary$method_id, summary$tol
  ), ]
}

accuracy.comparison = function(results, tols) {
  case_columns = c("section", "nu", "n", "d", "k_eff", "range_multiplier", "nugget", "regime", "seed")
  method_columns = c(case_columns, "knob", "loglik", "status")

  full = results[results$method_id == "restgp_full", method_columns]
  names(full) = c(case_columns, "r", "loglik_full", "status_full")

  pp = results[results$method_id == "restgp_pp", c(method_columns, "pp_seed")]
  names(pp) = c(case_columns, "r", "loglik_pp", "status_pp", "pp_seed")

  vecchia = results[results$method_id == "ssgv", method_columns]
  names(vecchia) = c(case_columns, "m", "loglik_vecchia", "status_vecchia")

  exact = unique(results[c(case_columns, "exact")])
  comparison = merge(full, pp, by = c(case_columns, "r"), all = TRUE)
  comparison = merge(comparison, exact, by = case_columns, all.x = TRUE)
  comparison$m = comparison$r
  comparison = merge(comparison, vecchia, by = c(case_columns, "m"), all.x = TRUE)
  comparison$m_over_r = 1L

  comparison$complete = with(
    comparison,
    is.finite(loglik_full) & is.finite(loglik_pp) &
      is.finite(loglik_vecchia) & is.finite(exact) &
      status_full == "ok" & status_pp == "ok" & status_vecchia == "ok"
  )
  comparison$full_abs_gap_per_obs = with(comparison, abs(loglik_full - exact) / n)
  comparison$pp_abs_gap_per_obs = with(comparison, abs(loglik_pp - exact) / n)
  comparison$vecchia_abs_gap_per_obs = with(comparison, abs(loglik_vecchia - exact) / n)
  comparison$both_higher_than_vecchia = with(
    comparison,
    loglik_full > loglik_vecchia & loglik_pp > loglik_vecchia
  )
  comparison$both_closer_than_vecchia = with(
    comparison,
    full_abs_gap_per_obs < vecchia_abs_gap_per_obs &
      pp_abs_gap_per_obs < vecchia_abs_gap_per_obs
  )

  summaries = lapply(tols, function(tol) {
    comparison$tol = tol
    comparison$all_within_tol = with(
      comparison,
      complete & full_abs_gap_per_obs <= tol &
        pp_abs_gap_per_obs <= tol & vecchia_abs_gap_per_obs <= tol
    )
    comparison$requested_condition = with(
      comparison,
      all_within_tol & both_closer_than_vecchia
    )
    comparison
  })
  summary = do.call(rbind, summaries)
  summary[order(
    summary$n, summary$d, summary$k_eff, summary$nugget,
    summary$r, summary$m, summary$pp_seed, summary$tol
  ), ]
}
