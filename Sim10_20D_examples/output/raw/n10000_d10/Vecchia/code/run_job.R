# Main study (original inputs): Rscript R/run_job.R 10000 10 WN raw
# Separate rescaled study:     Rscript R/run_job.R 10000 10 WN scaled
# Short validation only:       Rscript R/run_job.R 120 10 WN raw smoke
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
cat("Study directory:", root, "\n")
cat("Working directory:", getwd(), "\n")
flush.console()
source(file.path(root, "R", "settings.R"))
source(file.path(root, "R", "data.R"))
source(file.path(root, "R", "score_predictions.R"))

args <- commandArgs(TRUE)
if (!length(args) %in% c(4L, 5L)) {
  stop("Usage: Rscript R/run_job.R n d method raw|scaled [smoke]")
}
method <- match.arg(args[3], METHODS)
variant <- match.arg(args[4], VARIANTS)
smoke <- length(args) == 5L && identical(args[5], "smoke")
if (length(args) == 5L && !smoke) stop("The only optional argument is 'smoke'.")
s <- study_settings(as.integer(args[1]), as.integer(args[2]), variant, smoke)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
Sys.setenv(OMP_NUM_THREADS = ncores, OPENBLAS_NUM_THREADS = 1,
           MKL_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1)
source(file.path(root, "R", paste0("model_", tolower(method), ".R")))
fit_method <- switch(method, Full = fit_full, PP_middle = fit_pp_middle,
                     WN = fit_wn, WN_middle = fit_wn_middle,
                     Vecchia = fit_vecchia, BART = fit_bart,
                     WN_PMMH = fit_wn_pmmh)

# Each task owns its directory. No shared data file, prerequisite job or writer.
output_root <- if (smoke) file.path(root, "validation", "smoke") else file.path(root, "output")
parent <- file.path(output_root, variant, sprintf("n%d_d%d", s$n, s$d))
folder <- file.path(parent, method)
cat("Output directory:", folder, "\n")
flush.console()
if (!dir.exists(folder) && !dir.create(folder, recursive = TRUE)) {
  stop("Could not create output directory: ", folder)
}
# Reuse this method's folder. Keep the last complete result as a backup, but
# do not let the notebook mistake it for the result of the new run.
result_file <- file.path(folder, "result.rds")
if (file.exists(result_file)) {
  previous <- file.path(folder, "result.previous.rds")
  if (!file.rename(result_file, previous)) {
    stop("Could not preserve the previous result: ", previous)
  }
  cat("Previous completed result:", previous, "\n")
}
# Remove generated files from earlier attempts; other methods are untouched.
unlink(file.path(folder, c("result.pending", "predictions.csv", "metrics.csv",
  "diagnostics.rds", "chain.csv", "parameter_summary.csv", "pmmh_incomplete.rds")))
# ResTree supplies the common scores, including for Vecchia and BART.
packages <- c("DiceDesign", "ResTree", if (method == "BART") "BayesTree" else "GpGp")
package_info <- data.frame(package = packages,
  version = vapply(packages, function(p) as.character(utils::packageVersion(p)), ""),
  library = vapply(packages, find.package, ""))
print(package_info)
saveRDS(s, file.path(folder, "settings.rds"))
code_files <- list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)
code_hashes <- tools::md5sum(code_files)
dir.create(file.path(folder, "code"), showWarnings = FALSE)
invisible(file.copy(code_files, file.path(folder, "code"), overwrite = TRUE))

# Data generation and response centering are outside the model timers.
cat("Generating and saving training/test data...\n")
flush.console()
dat <- make_data(s$n, s$n_test, s$d, s$seed)
saveRDS(dat, file.path(folder, "data.rds"))
y_center <- mean(dat$y)
cat(method, variant, "n =", s$n, "d =", s$d, "started", format(Sys.time()), "\n")
flush.console()
out <- fit_method(dat$X, dat$y - y_center, dat$X_test, s, ncores)
# Preserve a shortened chain without marking the requested run complete.
if (!is.null(out$diagnostics$mcmc) &&
    nrow(out$diagnostics$mcmc$chain) != s$pmmh_iter) {
  saveRDS(out$diagnostics, file.path(folder, "pmmh_incomplete.rds"))
  stop("PMMH ended early; saved pmmh_incomplete.rds.")
}
stopifnot(length(out$mean) == s$n_test, length(out$variance) == s$n_test,
          all(is.finite(out$mean)), all(is.finite(out$variance)), all(out$variance > 0))
cat("Fitting and prediction finished; saving this method's predictions...\n")
flush.console()

# All exported predictions are back in the original response units.
predictions <- data.frame(test_id = seq_len(s$n_test), observed = dat$y_test,
  truth = dat$f_test, noise_sd = dat$sigma_test,
  mean = out$mean + y_center, sd = sqrt(out$variance), dat$X_test,
  check.names = FALSE)
predictions$lower <- predictions$mean - qnorm(0.975) * predictions$sd
predictions$upper <- predictions$mean + qnorm(0.975) * predictions$sd
predictions$latent_lower <- if (method == "BART") out$latent_lower + y_center else NA_real_
predictions$latent_upper <- if (method == "BART") out$latent_upper + y_center else NA_real_

# Save this method immediately, before scoring. Other jobs can still be running.
write.csv(predictions, file.path(folder, "predictions.csv"), row.names = FALSE)
saveRDS(out$diagnostics, file.path(folder, "diagnostics.rds"))
writeLines(capture.output(sessionInfo()), file.path(folder, "sessionInfo.txt"))
if (!is.null(out$diagnostics$mcmc)) {
  chain <- out$diagnostics$mcmc$chain
  write.csv(data.frame(iteration = seq_len(nrow(chain)), chain),
            file.path(folder, "chain.csv"), row.names = FALSE)
  write.csv(out$diagnostics$fit$summary,
            file.path(folder, "parameter_summary.csv"), row.names = FALSE)
}
cat("Predictions saved:", file.path(folder, "predictions.csv"), "\n")
cat("Computing scores for", method, "...\n")
flush.console()

# Evaluation and file I/O are not charged to fitting or prediction.
# Gaussian scores for every method; exact mixture scores where the prediction
# object is available (Full, PP_middle, WN).
# Exact lookup: `$prediction` would partially match `prediction_seconds`
# for methods that return only predictive means and variances.
scores <- score_predictions(predictions, out[["prediction"]], y_center)
metrics <- cbind(data.frame(method = method, variant = variant, n = s$n, n_test = s$n_test,
  d = s$d, replicate = s$replicate, seed = s$seed, smoke = smoke,
  fit_seconds = out$fit_seconds, prediction_seconds = out$prediction_seconds,
  total_seconds = out$total_seconds,
  initialization_seconds = out$initialization_seconds,
  particles = out$particles, convergence = out$convergence,
  ESS_final = if (is.null(out$diagnostics$fit$ESS_final)) NA_real_ else out$diagnostics$fit$ESS_final,
  estimator = out$estimator), scores)
result <- list(settings = s, metrics = metrics, predictions = predictions,
  diagnostics = out$diagnostics, y_center = y_center,
  packages = package_info, code_hashes = code_hashes,
  data_hash = unname(tools::md5sum(file.path(folder, "data.rds"))),
  session = sessionInfo())
cat("Saving scores and completed result for", method, "...\n")
flush.console()
write.csv(metrics, file.path(folder, "metrics.csv"), row.names = FALSE)
# Publish the completion marker only after the result has been fully written.
pending <- file.path(folder, "result.pending")
saveRDS(result, pending)
if (!file.rename(pending, result_file)) {
  stop("Could not publish the completed result: ", folder)
}
print(metrics)
cat("Saved", folder, "\n")
