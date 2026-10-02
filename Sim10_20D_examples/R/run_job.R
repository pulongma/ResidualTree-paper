# Run one method: Rscript R/run_job.R 10000 10 WN

# 1. Read the command and load the selected model.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(dirname(normalizePath(script, mustWork = TRUE)))
source(file.path(root, "R", "settings.R"))
source(file.path(root, "R", "data.R"))
source(file.path(root, "R", "score_predictions.R"))

args <- commandArgs(TRUE)
# Accept older batch commands with a final "raw" argument.
if (length(args) == 4L && args[4] == "raw") args <- args[1:3]
if (length(args) != 3L) stop("Usage: Rscript R/run_job.R n d method")
method <- match.arg(args[3], METHODS)
settings <- study_settings(as.integer(args[1]), as.integer(args[2]))
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
Sys.setenv(
  OMP_NUM_THREADS = ncores, OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1, VECLIB_MAXIMUM_THREADS = 1
)
source(file.path(root, "R", paste0("model_", tolower(method), ".R")))
fit_method <- switch(method,
  Full = fit_full,
  PP_middle = fit_pp_middle,
  WN = fit_wn,
  Vecchia = fit_vecchia,
  BART = fit_bart,
  WN_PMMH = fit_wn_pmmh
)

# 2. Prepare this method's folder and save the settings and code.
method_dir <- file.path(
  root, "output", "raw",
  sprintf("n%d_d%d", settings$n, settings$d), method
)
cat("Output directory:", method_dir, "\n")
flush.console()
if (!dir.exists(method_dir) && !dir.create(method_dir, recursive = TRUE)) {
  stop("Could not create output directory: ", method_dir)
}
# Keep the previous completed result when reusing this folder.
result_file <- file.path(method_dir, "result.rds")
if (file.exists(result_file) &&
  !file.rename(result_file, file.path(method_dir, "result.previous.rds"))) {
  stop("Could not preserve the previous result: ", result_file)
}
unlink(file.path(method_dir, c(
  "result.pending", "predictions.csv", "metrics.csv", "diagnostics.rds",
  "chain.csv", "parameter_summary.csv", "pmmh_incomplete.rds"
)))
packages <- c("DiceDesign", "ResTree", if (method == "BART") "BayesTree" else "GpGp")
package_info <- data.frame(
  package = packages,
  version = vapply(packages, function(p) as.character(utils::packageVersion(p)), ""),
  library = vapply(packages, find.package, "")
)
print(package_info)
saveRDS(settings, file.path(method_dir, "settings.rds"))
code_files <- list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)
code_hashes <- tools::md5sum(code_files)
dir.create(file.path(method_dir, "code"), showWarnings = FALSE)
invisible(file.copy(code_files, file.path(method_dir, "code"), overwrite = TRUE))

# 3. Generate data, center the training responses, and fit the model.
# Each model times only fitting (including initialization) and prediction.
data <- make_data(settings$n, settings$n_test, settings$d, settings$seed)
saveRDS(data, file.path(method_dir, "data.rds"))
y_center <- mean(data$y)
cat("Starting", method, "at", format(Sys.time()), "\n")
flush.console()
model_result <- fit_method(data$X, data$y - y_center, data$X_test, settings, ncores)
mcmc <- model_result$diagnostics$mcmc
if (!is.null(mcmc) && nrow(mcmc$chain) != settings$pmmh_iter) {
  saveRDS(model_result$diagnostics, file.path(method_dir, "pmmh_incomplete.rds"))
  stop("PMMH ended early; saved pmmh_incomplete.rds.")
}
stopifnot(
  length(model_result$mean) == settings$n_test,
  length(model_result$variance) == settings$n_test,
  all(is.finite(model_result$mean)), all(is.finite(model_result$variance)),
  all(model_result$variance > 0)
)

# 4. Save predictions immediately, before scoring or waiting for other jobs.
predictions <- data.frame(
  test_id = seq_len(settings$n_test), observed = data$y_test,
  truth = data$f_test, noise_sd = data$sigma_test,
  mean = model_result$mean + y_center, sd = sqrt(model_result$variance),
  data$X_test, check.names = FALSE
)
predictions$lower <- predictions$mean - qnorm(0.975) * predictions$sd
predictions$upper <- predictions$mean + qnorm(0.975) * predictions$sd
predictions$latent_lower <- if (method == "BART") model_result$latent_lower + y_center else NA_real_
predictions$latent_upper <- if (method == "BART") model_result$latent_upper + y_center else NA_real_
write.csv(predictions, file.path(method_dir, "predictions.csv"), row.names = FALSE)
saveRDS(model_result$diagnostics, file.path(method_dir, "diagnostics.rds"))
writeLines(capture.output(sessionInfo()), file.path(method_dir, "sessionInfo.txt"))
if (!is.null(mcmc)) {
  write.csv(data.frame(iteration = seq_len(nrow(mcmc$chain)), mcmc$chain),
    file.path(method_dir, "chain.csv"),
    row.names = FALSE
  )
  write.csv(model_result$diagnostics$fit$summary,
    file.path(method_dir, "parameter_summary.csv"),
    row.names = FALSE
  )
}
cat("Predictions saved; computing scores.\n")
flush.console()

# 5. Save scores and mark this run complete.
# Exact lookup avoids partially matching the name "prediction_seconds".
scores <- score_predictions(predictions, model_result[["prediction"]], y_center)
ess <- model_result$diagnostics$fit$ESS_final
metrics <- cbind(data.frame(
  method = method, variant = settings$variant, n = settings$n, n_test = settings$n_test,
  d = settings$d, replicate = settings$replicate, seed = settings$seed,
  fit_seconds = model_result$fit_seconds, prediction_seconds = model_result$prediction_seconds,
  total_seconds = model_result$total_seconds,
  initialization_seconds = model_result$initialization_seconds,
  particles = model_result$particles, convergence = model_result$convergence,
  ESS_final = if (is.null(ess)) NA_real_ else ess,
  estimator = model_result$estimator
), scores)
result <- list(
  settings = settings, metrics = metrics, predictions = predictions,
  diagnostics = model_result$diagnostics, y_center = y_center,
  packages = package_info, code_hashes = code_hashes,
  data_hash = unname(tools::md5sum(file.path(method_dir, "data.rds"))),
  session = sessionInfo()
)
write.csv(metrics, file.path(method_dir, "metrics.csv"), row.names = FALSE)
# Only a fully written result receives the filename read by the notebook.
pending <- file.path(method_dir, "result.pending")
saveRDS(result, pending)
if (!file.rename(pending, result_file)) stop("Could not publish: ", result_file)
print(metrics)
cat("Saved", result_file, "\n")
