timing.protocol = "elapsed_setup_likelihood_v2"

timing.context = function(ncores = 1L) {
  threads = Sys.getenv(c(
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"
  ))
  paste(R.version.string, Sys.info()["machine"], extSoftVersion()["BLAS"],
    La_library(), paste0("ncores=", ncores),
    paste(names(threads), threads, sep = "=", collapse = ";"),
    paste0("ResTree=", packageVersion("ResTree")),
    paste0("GpGp=", packageVersion("GpGp")),
    sep = " | "
  )
}
timed = function(expr) {
  start_time = proc.time()[["elapsed"]]
  value = force(expr)
  elapsed_seconds = as.numeric(proc.time()[["elapsed"]] - start_time)
  list(value = value, seconds = elapsed_seconds)
}
