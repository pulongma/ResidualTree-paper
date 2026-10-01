sample_sizes = c(2000L, 5000L, 10000L)
dimensions = c(2L, 5L, 10L, 20L)
minimum_ranges = 1
nugget_values = 0.1
knot_counts = 100L
pp_seed_count = 3L
accuracy_tolerances = c(1e-3, 1e-4)
output_suffix = "_s1_r100"

file_argument = grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_argument) > 0) {
  script_path = sub("^--file=", "", file_argument[1])
  script.dir = dirname(normalizePath(script_path))
} else {
  script.dir = getwd()
}

source(file.path(script.dir, "helpers", "cluster_config.R"))
config = cluster.config(
  sample_sizes, dimensions, minimum_ranges, nugget_values,
  knot_counts, pp_seed_count, accuracy_tolerances, output_suffix
)
options(hd.config = config)
source(file.path(script.dir, "01_run_study.R"))
