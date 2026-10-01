# Likelihood approximation: production reproduction

This release contains the two production runs used by the plotting notebooks.
Run both from this folder. No preliminary experiments are needed.

## Software and resources

The recorded production environment used R 4.4.3, ResTree 0.1.8, GpGp 0.5.1
and DiceDesign 1.10 on Linux with OpenBLAS 0.3.27. Use these versions to
reproduce the numerical results. The figures additionally require dplyr,
ggplot2, scales and an R notebook kernel (IRkernel). The simulation records
its package versions and numerical environment in `output/session_cluster_*.txt`.

The original allocation used 10 CPUs and 64 GB of memory. Exact-GP generation
and likelihood evaluation use dense matrices; the configured memory guard
remains 50 GB. Timing results depend on hardware and the BLAS configuration.

## Run the production simulations

```bash
cd /path/to/Loglik_approx
export RESTREE_NCORES=10 OMP_NUM_THREADS=10
export OPENBLAS_NUM_THREADS=10 MKL_NUM_THREADS=10
export BLIS_NUM_THREADS=10 VECLIB_MAXIMUM_THREADS=10
export OMP_DYNAMIC=FALSE OMP_PLACES=cores OMP_PROC_BIND=spread
Rscript R/01_run_cluster.R
Rscript R/01_run_cluster_r100.R
```

Each command evaluates n = 2,000, 5,000, 10,000 and d = 2, 5, 10, 20:
12 data settings and 60 approximation evaluations, including all three PP
seeds. The first command uses r = m = 30; the second uses r = m = 100.
They share the same generated observations and exact-data caches.

All generating and method settings are unchanged: mean 0, process variance 1,
Matérn smoothness 1.5, nugget variance 0.1, coordinate ranges
rho_j = 4^((j - 1)/(d - 1)) from 1 to 4, and base seed 2023. Full and PP use
middle cuts and sequential maximin knots. PP uses stratified resampling,
alpha = 1, prior_rho = 0.05, prior_beta = 0, H = ceiling(log2(n/r)) and 100H
particles. Vecchia uses m = r, maximin ordering and Euclidean nearest
predecessors in the original coordinates. All likelihoods are evaluated at
the true covariance parameters.

The production seed formula, result keys, data-cache names and numerical
calculations are preserved. RDS/CSV checkpoints are saved after each method
evaluation; rerunning a command resumes its completed evaluations.

## Reproduce the figures

| Budget | Result file | Notebook |
|---|---|---|
| r = m = 30 | `output/N_s1_r30.rds` | `plot_likelihood.ipynb` |
| r = m = 100 | `output/N_s1_r100.rds` | `plot_likelihood_r100.ipynb` |

When running simulations on a cluster, copy the matching RDS file into the
local `output/` directory. Open the notebook with this folder as its working
directory, then restart the kernel and run all cells. Each notebook displays
and saves three figures: likelihood, relative log-likelihood error and
approximation time against dimension. The PDFs are written to `figure/`.
Selections, labels, titles and output names are editable at the top.

The error axis displays abs(approximate - exact) / abs(exact) as a fraction:
0.01 means 1% and 0.1 means 10%. Timing plots show Full, PP and Vecchia.
Exact-GP timings remain in the result files. Full/PP timing includes model
construction and likelihood evaluation; Vecchia includes ordering, neighbor
search and likelihood evaluation. Exact timing includes fresh covariance
construction, Cholesky factorization and the likelihood solve. Simulation and
file I/O are excluded.

This folder contains code and notebook previews, not the production result
files. Generate those files with the commands above or supply the saved
production outputs. Embedded previews were generated from the saved results.

## Code layout

- `R/01_run_cluster.R` and `R/01_run_cluster_r100.R`: production grids.
- `R/api_contract.R`: shared method settings and covariance agreement checks.
- `R/01_run_study.R`: data-setting loop, checkpoints and session record.
- `R/00_setup.R`: libraries and shared helpers.
- `R/helpers/`: configuration, data generation, exact likelihood,
  approximation methods, timing, per-setting execution and summaries.
