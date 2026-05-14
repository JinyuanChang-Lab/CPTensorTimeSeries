# CP-Factorization for High Dimensional Tensor Time Series and Double Projection Iterations

This repository contains the R code used to reproduce the simulation studies and the real data analysis for the paper **CP-Factorization for High Dimensional Tensor Time Series and Double Projection Iterations**.

The code is organized into three main folders:

- `Simulation/`: simulation experiments in the main paper and supplementary material.
- `Application_AirPollution/`: code and data resources for the Beijing air-pollution application.
- `Application_FamaFrenchReturn/`: code and data resources for the Fama-French 100 portfolios return application.

All scripts write their results to the current working directory. To keep the relative `source("cp_cciso.R")` calls valid, run each script from inside its own folder.

---

## Repository Structure

### Simulation/

This folder contains the simulation code for the numerical experiments in the main paper and supplementary material.

#### Helper code

- `cp_cciso.R`  
  Implements auxiliary methods used for comparison, including cPCA, RP-PCA, CC-ISO and HOPE.

#### Main simulation scripts and outputs

| Script | Output file(s) | Reproduces |
|---|---|---|
| `simulation_baseline.R` | `table_rank_estimation.csv` | Table 1 in the main paper. |
| `simulation_baseline.R` | `table_rho_Ahat_A.csv` | Table 2 in the main paper. |
| `simulation_baseline.R` | `iterative_steps_selected_methods.pdf` | Figure 1 in the main paper. |
| `simulation_baseline.R` | `Pro.iter_estimated_h1_hist.pdf` | Figure 2 in the main paper. |
| `simulation_baseline.R` | `Pro.iter_estimated_h2_hist.pdf` | Figure 3 in the main paper. |
| `simulation_runtime_benchmark.R` | `runtime_benchmark.pdf` | Figure 4 in the main paper. |
| `simulation_variance_estimation.R` | `rmse_iter_variance.csv` | Table T6 in the supplementary material. |
| `simulation_baseline.R` | `table_common_component_error.csv` | Table T7 in the supplementary material. |
| `simulation_rank_misspecification_table.R` | `table_rank_estimation_stage.csv` | Table T8 in the supplementary material. |
| `simulation_ar_error.R` | `table_rho_Ahat_A_ARerror.csv` | Table T9 in the supplementary material. |
| `simulation_heavy_tail_error.R` | `table_rho_Ahat_A_heavytailed.csv` | Table T10 in the supplementary material. |
| `simulation_runtime_benchmark.R` | `ram_benchmark.pdf` | Figure F1 in the supplementary material. |
| `simulation_k_robustness.R` | `Krobust_accuracy.pdf` | Figure F2 in the supplementary material. |
| `simulation_k_robustness.R` | `Krobust_error.pdf` | Figure F3 in the supplementary material. |
| `simulation_rank_misspecification_figure.R` | `rtilde_selected_methods.pdf` | Figure F11 in the supplementary material. |

---

### Application_AirPollution

The `Application_AirPollution/` folder contains the code and data resources for the real data analysis of the Beijing air-pollution tensor time series.

#### Files

- `application_airpollution.R`  
  Runs the full empirical analysis, including model estimation, loading inference, comparison with competing methods, robustness checks using winsorized data, and plotting.

- `cp_cciso.R`  
  Helper implementation of comparison methods used by `application_airpollution.R`.

- `gadm41_CHN_3.json`  
  China administrative-boundary GeoJSON file used for the station-mode map figures.

#### Real-data outputs

| Script | Output file(s) | Reproduces |
|---|---|---|
| `application_airpollution.R` | `table_loading_a2_main.csv` | Table 3 in the main paper. |
| `application_airpollution.R` | `beijing_map_a1.pdf`, `beijing_map_a2.pdf` | Figure 5 in the main paper. |
| `application_airpollution.R` | `mode3_a1_Pro.iter.pdf`, `mode3_a2_Pro.iter.pdf` | Figure 6 in the main paper. |
| `application_airpollution.R` | `timeseries_factor1.pdf`, `timeseries_factor2.pdf` | Figure 7 in the main paper. |
| `application_airpollution.R` | `table_loading_a2_all.csv` | Table T1 in the supplementary material. |
| `application_airpollution.R` | `table_loading_a1_all.csv` | Table T2 in the supplementary material. |
| `application_airpollution.R` | `table_moment_test.csv` | Table T3 in the supplementary material. |
| `application_airpollution.R` | `table_loading_a2_robust.csv` | Table T4 in the supplementary material. |
| `application_airpollution.R` | `Y_time_series_airpollution.pdf` | Figure F4 in the supplementary material. |
| `application_airpollution.R` | `mode3_a1_compare.pdf`, `mode3_a2_compare.pdf` | Figure F5 in the supplementary material. |
| `application_airpollution.R` | `factor1_timeseries_compare.pdf`, `factor2_timeseries_compare.pdf` | Figure F6 in the supplementary material. |
| `application_airpollution.R` | `beijing_map_a1_robust.pdf`, `beijing_map_a2_robust.pdf` | Figure F7 in the supplementary material. |
| `application_airpollution.R` | `mode3_a1_robust.pdf`, `mode3_a2_robust.pdf` | Figure F8 in the supplementary material. |
| `application_airpollution.R` | `timeseries_factor1_robust.pdf`, `timeseries_factor2_robust.pdf` | Figure F9 in the supplementary material. |

---

### Application_FamaFrenchReturn

The `Application_FamaFrenchReturn/` folder contains the code and data resources for the Fama-French 100 portfolios return application.

#### Files

- `application_famafrenchreturn.R`  
  Runs the rolling forecast comparison for the Fama-French 100 portfolios return tensor time series.

- `tensor_cp_functions.R`, `CP_functions.R`, `cp_cciso.R`, `tools4cp.cpp`  
  Helper functions used by `application_famafrenchreturn.R`.

- `100_Portfolios.CSV`  
  Data used in the rolling forecast application.

#### Real-data outputs

| Script | Output file(s) | Reproduces |
|---|---|---|
| `application_famafrenchreturn.R` | `rolling_final_display_table.csv` | Table T5 in the supplementary material. |
| `application_famafrenchreturn.R` | `Y_time_series_100.pdf` | Figure F10 in the supplementary material. |

---

## How to Run

The function `CP_TTS` in the package `HDTSA` implements the proposed method. The scripts also require the following R packages:

```r
install.packages(c(
  "HDTSA", "pracma", "Rcpp", "RcppEigen", "foreach", "doParallel", "jointDiag",
  "R.utils", "tensor", "rTensor", "MASS", "vars", "ggplot2", "cowplot",
  "gridExtra", "dplyr", "tidyr", "sf", "zoo", "xts", "scales",
  "ggrepel", "tibble", "tensorTS", "OLCPM"
))
```

Run simulation scripts from `Simulation/`, for example:

```r
setwd("Simulation")
source("simulation_baseline.R")
```

Run the real data analysis from `Application_AirPollution/`:

```r
setwd("Application_AirPollution")
source("application_airpollution.R")
```

Run the rolling forecast analysis from `Application_FamaFrenchReturn/`:

```r
setwd("Application_FamaFrenchReturn")
source("application_famafrenchreturn.R")
```

For the simulation studies, we recommend running the code on a server with 128 CPU cores and 512 GB RAM. Otherwise, please reset the number of cores used by `makeCluster()` before running the code.

Each simulation script also saves the complete workspace and full intermediate results in an `.RData` file.
