# Deterministic multiscale response on the unit cube with a regime change on
# x1 and input-dependent noise on x2. Designed 2026-09-23 (see README.md).
#
#   r(x)     = || (x3, x4, x5) - (0.3, 0.6, 0.4) ||           3-D radial coordinate
#   g(t)     = plogis((t - 0.5) / 0.02)                       sharp gate on x1
#   f(x)     = A * [ (1 - g(x1)) sin(2 pi nu0 r) + g(x1) sin(2 pi nu1 r + phi) ]
#              + J * g(x1) + b * sum_j cos(pi x_j) / sqrt(d)
#   sigma(x) = s_lo + (s_hi - s_lo) * plogis((x2 - 0.5) / 0.05)
#
# with A = 2, nu0 = 2, nu1 = 3, phi = 1, J = 2, b = 0.5, s_lo = 0.05, s_hi = 1.
# The radial wave is a smooth three-way interaction; every coordinate enters
# the mean through the cosine background; x1 switches the wave frequency and
# shifts the level; x2 switches the noise standard deviation.
multiscale_constants <- list(
  wave_coords = 3:5, wave_center = c(0.3, 0.6, 0.4),
  amplitude = 2, nu0 = 2, nu1 = 3, phase = 1, jump = 2, background = 0.5,
  gate_width = 0.02, noise_lo = 0.05, noise_hi = 1, noise_width = 0.05
)

multiscale_mean <- function(U, k = multiscale_constants) {
  stopifnot(ncol(U) >= max(k$wave_coords))
  r <- sqrt(rowSums(sweep(U[, k$wave_coords, drop = FALSE], 2, k$wave_center)^2))
  g <- plogis((U[, 1] - 0.5) / k$gate_width)
  wave <- (1 - g) * sin(2 * pi * k$nu0 * r) + g * sin(2 * pi * k$nu1 * r + k$phase)
  k$amplitude * wave + k$jump * g + k$background * rowSums(cos(pi * U)) / sqrt(ncol(U))
}

noise_sd <- function(U, k = multiscale_constants) {
  k$noise_lo + (k$noise_hi - k$noise_lo) * plogis((U[, 2] - 0.5) / k$noise_width)
}

make_data <- function(n, n_test, d, seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  X <- DiceDesign::lhsDesign(n, dimension = d, randomized = TRUE,
                             seed = seed)$design

  # Independent LHS for testing. A disjoint first coordinate guarantees no
  # repeated training row. Redraw the whole design in the unlikely collision.
  test_seed <- seed + 1000000L
  repeat {
    X_test <- DiceDesign::lhsDesign(n_test, dimension = d, randomized = TRUE,
                                    seed = test_seed)$design
    if (!any(X_test[, 1] %in% X[, 1])) break
    test_seed <- test_seed + 1L
  }
  colnames(X) <- colnames(X_test) <- paste0("x", seq_len(d))
  f <- multiscale_mean(X)
  f_test <- multiscale_mean(X_test)
  sigma <- noise_sd(X)
  sigma_test <- noise_sd(X_test)

  # Separate noise seeds keep both designs independent of noise generation.
  set.seed(seed + 2000000L)
  y <- f + rnorm(n, sd = sigma)
  set.seed(seed + 3000000L)
  y_test <- f_test + rnorm(n_test, sd = sigma_test)
  list(X = X, y = y, f = f, sigma = sigma,
       X_test = X_test, y_test = y_test, f_test = f_test, sigma_test = sigma_test,
       seed = seed, test_seed = test_seed)
}
