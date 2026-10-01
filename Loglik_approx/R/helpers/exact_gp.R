dense.matern15.ard = function(X, rho, nugget, sig2 = 1, block.rows = 256L) {
  X = as.matrix(X)
  n = nrow(X)
  Z = sweep(X, 2, rho, "/")
  zz = rowSums(Z^2)
  C = matrix(0, n, n)
  for (lo in seq.int(1L, n, by = block.rows)) {
    ii = lo:min(n, lo + block.rows - 1L)
    D2 = pmax(outer(zz[ii], zz, "+") - 2 * tcrossprod(Z[ii, , drop = FALSE], Z), 0)
    DD = sqrt(D2)
    C[ii, ] = sig2 * (1 + DD) * exp(-DD)
  }
  diag(C) = diag(C) + sig2 * nugget
  C
}

exact.gp.likelihood = function(X, y, rho, nugget, sig2 = 1) {
  covariance = timed(dense.matern15.ard(X, rho, nugget, sig2))
  factorisation = timed(chol(covariance$value))
  evaluation = timed({
    U = factorisation$value
    whitened_y = forwardsolve(U, y, upper.tri = TRUE, transpose = TRUE)
    -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(whitened_y^2))
  })
  list(
    loglik = as.numeric(evaluation$value),
    seconds_covariance = covariance$seconds,
    seconds_cholesky = factorisation$seconds,
    seconds_solve = evaluation$seconds,
    seconds = covariance$seconds + factorisation$seconds + evaluation$seconds
  )
}

benchmark.exact.case = function(case, context) {
  gc(verbose = FALSE)
  value = exact.gp.likelihood(
    case$X, case$y, case$theta$range,
    case$theta$nugget, case$theta$sig2
  )
  if (abs(value$loglik - case$exact_loglik) > 1e-7 * max(1, abs(case$exact_loglik))) {
    stop("Timed exact likelihood disagrees with the simulation reference")
  }
  modifyList(case, list(
    seconds = value$seconds,
    seconds_covariance = value$seconds_covariance,
    seconds_cholesky = value$seconds_cholesky,
    seconds_solve = value$seconds_solve,
    evaluated_loglik = value$loglik,
    timing_protocol = timing.protocol, timing_context = context,
    timing_host = unname(Sys.info()["nodename"]),
    timing_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    timing_reused = FALSE
  ))
}

exact.gp.case = function(n, d, seed, rho, nugget, sig2, cache.dir, max.gb, rho_min,
                          context = timing.context()) {
  path = file.path(cache.dir, sprintf(
    "exact_nu15_n%d_d%d_%s_%s_seed%d.rds",
    n, d, regime.id(rho_min, nugget), range.order, seed
  ))
  theta = list(
    sig2 = sig2, range = as.numeric(rho), nugget = nugget, nu = known.nu, form = "ARD",
    k_eff = NA_integer_, range_order = range.order
  )
  if (file.exists(path)) {
    old = readRDS(path)
    if (identical(old$status, "ok") && isTRUE(all.equal(old$theta, theta, tolerance = 0))) {
      old$from_cache = TRUE
      same_timing = identical(old$timing_protocol, timing.protocol) &&
        identical(old$timing_context, context)
      if (same_timing) {
        old$timing_reused = TRUE
      } else {
        if (is.null(old$timing_protocol)) {
          old$legacy_generation_seconds = old$seconds
          old$seconds_simulation = NA_real_
        }
        old = benchmark.exact.case(old, context)
        saveRDS(old, path, compress = "xz")
      }
      return(old)
    }
  }
  projected.gb = 2.35 * 8 * as.double(n)^2 / 1024^3
  case = list(
    status = "skipped_memory", message = NA_character_, file = path, n = n, d = d,
    seed = seed, theta = theta, projected_peak_gb = projected.gb, max_gb = max.gb,
    from_cache = FALSE, seconds = NA_real_, seconds_simulation = NA_real_,
    blas = unname(extSoftVersion()["BLAS"]), lapack = La_library(),
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  if (projected.gb > max.gb) {
    case$message = sprintf("projected dense peak %.1f GB exceeds the guard %.1f GB", projected.gb, max.gb)
  } else {
    case = tryCatch(
      {
        simulation = timed({
          X = make.design(n, d, seed)
          C = dense.matern15.ard(X, rho, nugget, sig2)
          U = chol(C)
          rm(C)
          set.seed(seed + 1L)
          z = rnorm(n)
          y = as.numeric(crossprod(U, z))
          ll = -0.5 * (n * log(2 * pi) + 2 * sum(log(diag(U))) + sum(z^2))
          rm(U, z)
          list(X = X, y = y, exact_loglik = ll)
        })
        case = modifyList(case, c(
          simulation$value,
          list(status = "ok", seconds_simulation = simulation$seconds)
        ))
        benchmark.exact.case(case, context)
      },
      error = function(e) modifyList(case, list(status = "failed", message = conditionMessage(e)))
    )
  }
  saveRDS(case, path, compress = "xz")
  case
}
