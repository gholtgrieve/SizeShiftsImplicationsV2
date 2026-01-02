# SizeShiftsImplicationsV2

> Fishery management implications of demographic changes in Chinook salmon.

## Overview

**SizeShiftsImplicationsV2** is an R package that implements a management strategy evaluation (MSE) framework for exploring how shifts in age, sex, and length composition propagate through 
stock–recruit dynamics and harvest control rules. This package accompanies and extends analyses related to Ohlberger et al. (2025) and is designed to make those workflows easier to share, 
reproduce, and adapt.

## Relationship to the original project

This package **is based on, builds on, and reuses code from** the original repository by Jan Ohlberger: **`janohlberger/SizeShiftsImplications`**. 

### What’s new in V2 (differences in objective & implementation)

1. **R package structure & documentation** – Formalized as an installable package with a consistent API (`run_scenarios()`, `make_figures()`, etc.), a standard results object (`ssi_run`), and roxygen2 documentation. This makes the project easier to install, share, test, and cite.
2. **Parallel execution** – The simulation runner supports **parallel, scenario‑level execution** (via `future`/`future.apply`) with sensible defaults for worker counts and BLAS threading. Large grids complete much faster on multi‑core machines.
3. **Kuskokwim River analysis** – Adds a ready‑to‑run parameter profile and figure builders tailored to the Kuskokwim River, producing a parallel set of figures (A–E) for that system.
4. **Reproducibility by design** – Seeds and run metadata (e.g., RNG kind, BLAS thread caps, package version, scenario grid hash, timestamps) are recorded inside `ssi_run$parameters` so runs can be diagnosed and compared later.
5. **Unified results object** – All per‑scenario/per‑iteration outputs (e.g., SR parameters, fecundity/egg mass, S\_MSY trajectories, observations, implementation errors) are collated into a single `ssi_run` object and saved to disk alongside a scenario table (Excel). Downstream plotting functions use this object directly.
6. **Figure builders** – Programmatic functions to rebuild the Ohlberger (2–7) and Kuskokwim (A–E) figures, writing both **PDF** plots and the underlying **CSV** data used in each figure.
7. **Legacy conversion** – A helper to convert legacy outputs into the new `ssi_run` structure for side‑by‑side comparisons and regression testing.

> **Attribution**: If you use this package, please also acknowledge the original project and its authors.

## Installation

Requires **R ≥ 4.1**.

```r
# Option 1: with pak
install.packages("pak")
pak::pak("gholtgrieve/SizeShiftsImplicationsV2")

# Option 2: with remotes
install.packages("remotes")
remotes::install_github("gholtgrieve/SizeShiftsImplicationsV2")
```

## Quick start

Run a small, parallel test and build a couple of figures.

```r
library(SizeShiftsImplicationsV2)

# 1) Run a modest grid with reproducible seeds
run <- run_scenarios(
  scenarios  = "all",        # or a dplyr-style filter string, or c(ID1, ID2, ...)
  niter      = 50,
  seed       = "reproducible",  # or "random"
  params     = "Ohlberger",     # or "Kuskokwim" or a full list of parameters
  output_dir = file.path(here::here(), "outputs"),
  parallel   = TRUE,             # each worker runs one scenario
  workers    = NULL              # NULL = auto-detect (leaves 1 core free)
)

# 2) Recreate publication-style figures (PDF + CSV)
make_figures(figures = c("2","3"), data = run, figure_dir = here::here("figures"))
```

### Kuskokwim analysis

To run the Kuskokwim River configuration and produce its figures:

```r
run_kusko <- run_scenarios(
  scenarios = "all",
  niter     = 100,
  seed      = "reproducible",
  params    = "Kuskokwim",
  parallel  = TRUE
)
make_figures("Kuskokwim", data = run_kusko, figure_dir = here::here("figures"))
```

## The `ssi_run` object

A run returns (and saves) a **single artifact**:

- `meta`: timestamp, number of scenarios/iterations, seed mode/list, profile, output\_dir.
- `scenarios`: the exact scenario table that was executed.
- `results`: per‑scenario/per‑iteration lists including: `para`, `sr_sim`, `fec`, `egg`, `S_msy`, `data`, `obs`, `ret_by_age`, `meanSaA`, `propfemale`, `selectivity`, `MSY_Goals`, `impl_errors`.

A serialized copy is written to `outputs/run_<timestamp>.rds`, and the scenario table is also exported to `run_<timestamp>_scen.xlsx` for quick inspection. Figure builders write `FigureX.pdf` + matching `FigureX.csv` under `figures/`.

## Parallel performance & reproducibility tips

- **Seeds**: `seed = "reproducible"` uses common random numbers (CRN) across scenarios so iteration *k* is seeded identically across the grid. `seed = "random"` samples seeds for each iteration.
- **Workers**: If `workers = NULL`, the package auto‑detects cores, caps workers at the number of scenarios, and **leaves one core free** for the OS.
- **BLAS threads**: To avoid oversubscription in parallel runs, limit BLAS threads before running:
  ```r
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
  ```

## Figures

- **Ohlberger**: 2–7 (e.g., S\_MSY distributions with 50%/80% bands; harvest/return comparisons; YPR/DLM percentage differences relative to TRM; etc.).
- **Kuskokwim**: A–E (system‑specific metrics; probability of zero harvest versus harvest; etc.).

Each builder also saves the tidy **CSV** used to draw the figure to make manuscript tables/QA straightforward.

## Migrating legacy outputs

Use the included converter to produce a faithful `ssi_run` for from legacy output derived from the **`janohlberger/SizeShiftsImplications`** version of the model.

```r
ssi <- convert_legacy_to_ssi(path_to_Rdata_or_env)
```

## Contributing

Best efforts to address issues will be made as they arise but recognize this package is made available for reproducability purposes, not as maintained software.  

## Citing

If this package or its figures contributed to your work, please cite:

> Ohlberger, J., Schindler, D. E., & Staton, B. A. (2025). Accounting for salmon body size declines in fishery management can reduce conservation risks. *Fish and Fisheries*, 26, 113–130.[ ](https://doi.org/10.1111/faf.12869)[https://doi.org/10.1111/faf.12869](https://doi.org/10.1111/faf.12869)


## License

GPL-3.0. See `LICENSE`.
