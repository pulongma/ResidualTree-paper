vecchia.ref = function(X, y, rho, nugget, m, sig2 = 1) {
  m = min(as.integer(m), nrow(X) - 1L)
  ordering = timed({
    ord = GpGp::order_maxmin(X)
    list(X = X[ord, , drop = FALSE], y = y[ord])
  })
  ordered = ordering$value
  neighbours = timed(GpGp::find_ordered_nn(ordered$X, m = m))
  likelihood = timed(GpGp::vecchia_meanzero_loglik(
    c(sig2, rho, nugget), "matern15_scaledim", ordered$y, ordered$X, neighbours$value
  )$loglik)
  list(
    loglik = as.numeric(likelihood$value), m = m,
    seconds_ordering = ordering$seconds, seconds_neighbors = neighbours$seconds,
    seconds_model = ordering$seconds + neighbours$seconds,
    seconds = likelihood$seconds
  )
}

evaluate = function(job, X, y, theta, ranges, nugget, n, seed) {
  r = job$knob
  if (job$method_id == "restgp_full") {
    model_build = timed(restree_model(
      X, y,
      depth = tree.depth(n, r),
      r = r,
      leaf_model = "Full",
      design = knot.design,
      cut_method = cut.method
    ))
    likelihood = timed(restree_loglik(
      model_build$value, theta,
      ncores = ncores, diagnostics = TRUE
    ))

    return(list(
      loglik = as.numeric(likelihood$value),
      depth = model_build$value$depth,
      nparticles = NA_integer_,
      n_knot = as.numeric(attr(likelihood$value, "n_knot")),
      seconds_model = model_build$seconds,
      seconds = likelihood$seconds,
      seconds_ordering = NA_real_, seconds_neighbors = NA_real_
    ))
  }

  if (job$method_id == "restgp_pp") {
    depth = tree.depth(n, r)
    particle_count = nparticles.for(depth)
    sampler_seed = seed + 1000L * job$pp_seed + r

    model_build = timed(restree_model(
      X, y,
      depth = depth,
      r = r,
      leaf_model = "PP",
      design = knot.design,
      cut_method = cut.method,
      prior_rho = prior.rho,
      prior_beta = prior.beta
    ))
    likelihood = timed(restree_loglik(
      model_build$value, theta,
      nparticles = particle_count,
      resampling = resampling,
      temper_alpha = temper.alpha,
      seed = sampler_seed,
      ncores = ncores
    ))

    return(list(
      loglik = as.numeric(likelihood$value),
      depth = depth,
      nparticles = particle_count,
      n_knot = NA_real_,
      seconds_model = model_build$seconds,
      seconds = likelihood$seconds,
      seconds_ordering = NA_real_, seconds_neighbors = NA_real_
    ))
  }

  if (job$method_id == "ssgv") {
    set.seed(as.integer((as.double(seed) + 200000L + r) %% (.Machine$integer.max - 1)))
    likelihood = vecchia.ref(X, y, ranges, nugget, m = r, sig2 = known.sig2)
    return(list(
      loglik = likelihood$loglik,
      depth = NA_integer_,
      nparticles = NA_integer_,
      n_knot = NA_real_,
      seconds_model = likelihood$seconds_model,
      seconds = likelihood$seconds,
      seconds_ordering = likelihood$seconds_ordering,
      seconds_neighbors = likelihood$seconds_neighbors
    ))
  }
}
