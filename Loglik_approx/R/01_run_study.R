source(file.path(script.dir, "00_setup.R"))

config = getOption("hd.config")
if (is.null(config)) {
  stop("Run 01_run_cluster.R or 01_run_cluster_r100.R to set the production configuration.")
}

mode = config$mode
base.seed = config$base_seed
known.sig2 = 1
ncores = as.integer(config$ncores)
tols = config$tol
benchmark.context = timing.context(ncores)

out.dir = config$output_dir
if (is.null(out.dir)) {
  out.dir = file.path(script.dir, "..", "output")
}
cache.dir = file.path(out.dir, "exact_gp_cache")
dir.create(cache.dir, showWarnings = FALSE, recursive = TRUE)
output.suffix = config$output_suffix

stem = function(name) {
  file.path(out.dir, paste0(name, output.suffix, ".rds"))
}

load.table = function(section) {
  result_file = stem(section)
  if (file.exists(result_file)) readRDS(result_file) else NULL
}

write.summary = function(results, name) {
  result_file = stem(name)
  saveRDS(results, result_file)
  write.csv(results, sub("\\.rds$", ".csv", result_file), row.names = FALSE)
}

restree.min.version = config$min_package_version
require.restree.api()
cat("[api]", api.line(), "\n")
cat("[library]", find.package("ResTree"), "\n")
cat(
  "[env]", R.version.string, "| BLAS", unname(extSoftVersion()["BLAS"]),
  "| ncores", ncores, "\n"
)
cat("[timing]", timing.protocol, "| elapsed setup + likelihood; simulation excluded\n")

study_config = config$N
regimes = prepare.regimes(study_config$regimes)
write.summary(regimes, "N_regimes")
cat("[N] regimes:\n")
print(regimes)

results = load.table("N")
for (sample_size in sort(study_config$n)) {
  for (regime_index in seq_len(nrow(regimes))) {
    setting = regimes[regime_index, ]
    results = run.cell(
      sample_size, setting$d, setting$nugget, setting$rho_min,
      study_config, results
    )
  }
}

write.summary(required.knob(results, tols), "N_required_knob")
comparison = accuracy.comparison(results, tols)
write.summary(comparison, "N_accuracy_comparison")
cat("[N] Accuracy comparisons (each PP seed reported separately):\n")
print(aggregate(
  cbind(complete, all_within_tol, requested_condition) ~ tol + m_over_r,
  comparison, sum
))

session_file = file.path(out.dir, paste0("session_", mode, output.suffix, ".txt"))
writeLines(
  c(
    api.line(), paste("Timing protocol:", timing.protocol),
    paste("Timing environment:", benchmark.context),
    capture.output(str(config)), capture.output(sessionInfo())
  ),
  session_file
)
cat("[done]", format(Sys.time()), "\n")
