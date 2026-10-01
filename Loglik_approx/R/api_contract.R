known.nu = 1.5

restree.api = "matern15_sizeconstrained_full_2026-09-18"
restree.min.version = "0.1.8"

cut.method = "middle"
knot.design = "maximin"
prior.rho = 0.05
prior.beta = 0

resampling = "stratified"
temper.alpha = 1
nparticles.for = function(depth) 100L * max(1L, as.integer(depth))
tree.depth = function(n, r) max(1L, min(20L, as.integer(ceiling(log2(n / r)))))

r.max = 200L

range.order = "tight_first"

exact.gp.package = function(X, y, theta, ncores = 1L) {
  X = as.matrix(X)
  as.numeric(restree_loglik(
    restree_model(X, y,
      depth = 0L, r = nrow(X), leaf_model = "Full",
      design = knot.design, cut_method = cut.method
    ),
    theta,
    ncores = as.integer(ncores)
  ))
}
crosscheck.n.max = 2000L
crosscheck.rtol = 1e-7

dense.loglik = function(X, y, sig2, rho, nugget) {
  Z = sweep(as.matrix(X), 2, rho, "/")
  D = as.matrix(stats::dist(Z))
  C = sig2 * (1 + D) * exp(-D)
  diag(C) = diag(C) + sig2 * nugget
  U = chol(C)
  z = forwardsolve(t(U), y)
  -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(z^2))
}

require.restree.api = function() {
  fail = function(...) stop(..., ". Install the current ResTree source with R CMD INSTALL --preclean /path/to/ResTree into this R library.",
    call. = FALSE
  )
  v = utils::packageVersion("ResTree")
  if (v < restree.min.version) fail("ResTree ", v, " is older than ", restree.min.version)
  need.l = c("model", "theta", "tree", "nparticles", "resampling", "temper_alpha", "seed", "ncores", "diagnostics")
  need.m = c("x", "y", "depth", "r", "leaf_model", "design", "cut_method", "prior_rho", "prior_beta")
  miss = c(
    setdiff(need.l, names(formals(ResTree::restree_loglik))),
    setdiff(need.m, names(formals(ResTree::restree_model)))
  )
  if (length(miss)) fail("missing argument(s) in ResTree: ", paste(miss, collapse = ", "))

  set.seed(1L)
  n0 = 150L
  rho0 = c(0.15, 0.3, 0.6)
  nug0 = 0.05
  X0 = matrix(runif(n0 * 3L), n0, 3L)
  Z0 = sweep(X0, 2, rho0, "/")
  D0 = as.matrix(stats::dist(Z0))
  C0 = (1 + D0) * exp(-D0)
  diag(C0) = diag(C0) + nug0
  y0 = as.numeric(t(chol(C0)) %*% rnorm(n0))
  th0 = restree_theta(sig2 = 1, range = rho0, nugget = nug0, nu = known.nu, form = "ARD")
  m.dense = restree_model(X0, y0,
    depth = 0L, r = n0, leaf_model = "Full",
    design = knot.design, cut_method = cut.method
  )
  m.split = restree_model(X0, y0,
    depth = 0L, r = 20L, leaf_model = "Full",
    design = knot.design, cut_method = cut.method
  )
  if (m.dense$depth != 0L) fail("a Full model with r = n did not stay at depth 0")
  if (m.split$depth < 1L) fail(
    "a Full model with r = 20 < n = 150 did not split at depth = 0; ",
    "the installed package predates the 0.1.6 size-constrained Full rule"
  )
  ll.r = dense.loglik(X0, y0, 1, rho0, nug0)
  ll.pkg = as.numeric(restree_loglik(m.dense, th0, ncores = 1L))
  ord = GpGp::order_maxmin(X0)
  NN = GpGp::find_ordered_nn(X0[ord, ], m = n0 - 1L)
  ll.gp = GpGp::vecchia_meanzero_loglik(
    c(1, rho0, nug0), "matern15_scaledim",
    y0[ord], X0[ord, ], NN
  )$loglik
  tol = 1e-6 * max(1, abs(ll.r))
  if (abs(ll.pkg - ll.r) > tol) fail(sprintf("ResTree likelihood %.8f != dense R %.8f", ll.pkg, ll.r))
  if (abs(ll.gp - ll.r) > tol) fail(sprintf("GpGp likelihood at m = n - 1 %.8f != dense R %.8f", ll.gp, ll.r))
  invisible(TRUE)
}

api.line = function()
  sprintf(
    "ResTree %s | nu=1.5 | %s | cut_method=%s design=%s prior_rho=%g prior_beta=%g | resampling=%s temper_alpha=%g | particles=100*depth | r<=%d | range order %s",
    utils::packageVersion("ResTree"), restree.api, cut.method, knot.design,
    prior.rho, prior.beta, resampling, temper.alpha, r.max, range.order
  )
