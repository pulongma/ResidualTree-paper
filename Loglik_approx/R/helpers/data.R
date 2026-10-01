regime.id = function(rho_min, nugget) {
  sprintf("rho_min%.17g_nug%s", rho_min, format(nugget, trim = TRUE))
}

geometric.ranges = function(d, s) {
  s * 4^seq(0, 1, length.out = d)
}

make.design = function(n, d, seed) {
  X = DiceDesign::lhsDesign(n, d, seed = seed)$design
  pmin(pmax(X, 1e-4), 1 - 1e-4)
}

prepare.regimes = function(regimes) {
  regimes$k_eff = NA_integer_
  regimes$range_multiplier = NA_real_
  regimes$rho_max = 4 * regimes$rho_min
  regimes
}

case.seed = function(base.seed, n, d, nugget) {
  as.integer((base.seed + n + 100003 * d + 1009 * 3200L +
    7919 * round(1e6 * nugget)) %% (.Machine$integer.max - 1))
}

r.grid.for = function(r.grid, n) {
  budgets = sort(unique(as.integer(r.grid)))
  budgets[budgets <= r.max & budgets <= n / 2]
}
m.grid.for = function(r.grid, n) {
  budgets = sort(unique(as.integer(r.grid)))
  budgets[budgets < n]
}
