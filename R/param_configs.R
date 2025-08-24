#' Built-in parameter profiles (`param_configs`)
#'
#' A named list of **built-in parameter profiles** used to configure
#' `run_model()` and higher-level helpers. Each element is a list of
#' baseline values; **scenario-specific overrides** (e.g., `futureT`,
#' `harvmgmt`, `factorMSY`, `maxsel`, `sdsel`, and switches for
#' `agetrend`, `sizetrends`, `propFtrend`) are applied at runtime.
#'
#' These profiles are intended to be layered over a generic
#' [`default_config()`] and then modified with
#' [`apply_scenario_overrides()`].
#'
#' @format A named list with two elements:
#' \describe{
#'   \item{`Ohlberger`}{Baseline profile for the Ohlberger-style setup.}
#'   \item{`Kuskokwim`}{Baseline profile for the Kuskokwim setup.}
#' }
#'
#' Each profile list can include (non-exhaustive):
#' \itemize{
#'   \item Core switches: `sim_recruits`, `est_method`, `ricker_type`, `var`, `optimize`
#'   \item Time horizon: `nyi`, `nyh`, `ny`  (note: `futureT` is scenario-driven)
#'   \item Errors: `obserr`, `hobserr`, `harverr`
#'   \item Regimes: `reglength`, `regstr`
#'   \item Stock–recruit: `procerr`, `rho`, `alpha_mean`, `beta_mean`, `sr_corr`, `sr_parms_sd`
#'   \item Ages/structure: `ages`, `meanageini`, `sdage`
#'   \item Trend baselines (multiplied by scenario switches): `agetrend`, `sizetrends`, `propFtrend`
#'   \item Sex: `propF`, `agediff`
#'   \item Growth/size: `vonB_Linf`, `vonB_k`, `ocean0s`, `sdSaA`
#'   \item Allometry: `allometry`
#'   \item Alternative SR params: `alt_sr_param`
#' }
#'
#' @details
#' These are **baselines** only. For a full, ready-to-run config:
#' 1) Start with `default_config()`,
#' 2) Overlay one of `param_configs[[...]]`,
#' 3) Call `apply_scenario_overrides(cfg, scen_row)`, then
#' 4) Set `cfg$j`, `cfg$k`, and `cfg$seednum`.
#'
#' @examples
#' # Inspect available profiles
#' names(param_configs)
#'
#' # Build a config using the Ohlberger profile (pseudo-example):
#' # cfg <- utils::modifyList(default_config(), param_configs[["Ohlberger"]])
#' # cfg <- apply_scenario_overrides(cfg, scen_row)  # scen_row = one row from scenarios
#' # cfg$seednum <- 123; cfg$j <- 1L; cfg$k <- 1L
#' # run_model(cfg)
#'
#' # Use directly in a higher-level runner (if it accepts `params =`):
#' # run_scenarios("all", niter = 10, params = param_configs[["Kuskokwim"]])
#'
#' @seealso [default_config()], [apply_scenario_overrides()], [run_scenarios()], [run_model()]
#' @export
param_configs <- list(

  Ohlberger = list(
    # core switches
    sim_recruits = "eggmass",
    est_method   = "GLS",
    ricker_type  = "const_beta",
    var          = "both",
    optimize     = "logFmax",

    # time horizon
    nyi = 10, nyh = 50, ny = 110,

    # observation / implementation error
    obserr = 0.2, hobserr = 0.1, harverr = 0.15,

    # regime settings
    reglength = 20, regstr = 1,

    # management defaults (scenario sets harvmgmt/factorMSY at runtime)
    harvrate = 0.5, goalfreq = 10,

    # stock–recruit process
    procerr = 0.35,
    rho = 0.4,
    alpha_mean = 5,
    beta_mean  = 5e-5,
    sr_corr    = 0.0,
    sr_parms_sd= 0.2,

    # ages and structure
    ages = 1:9,
    meanageini = 5.5,
    sdage = 0.6,

    # trend baselines (multiplied by scenario switches)
    agetrend   = -0.4,
    sizetrends = c(0, 0, 30, 10, -30, -60, rep(-90, 3)),
    propF      = 0.45,
    propFtrend = -0.4,
    agediff    = 1,

    # growth / size
    vonB_Linf = 1200,
    vonB_k    = 0.325,
    ocean0s   = 150,
    sdSaA     = 0.01,

    # allometry
    allometry = c(
      size_ref = 800,
      fec_ref  = 6600,
      b_fec    = 2.4,
      egg_ref  = 916,
      b_eggs   = 4.8
    ),

    # alternative SR parameters
    alt_sr_param = c(
      alpha_N0   = 5.07,
      alpha_EASL = 0.00177,
      alpha_EMASL= 0.01036,
      beta_N0    = 8.6e-6,
      beta_EASL  = 2.973e-09,
      beta_EMASL = 1.6967e-08
    )
  ),

  Kuskokwim = list(
    # core switches
    sim_recruits = "eggmass",
    est_method   = "GLS",
    ricker_type  = "const_beta",
    var          = "both",
    optimize     = "logFmax",

    # time horizon
    nyi = 10, nyh = 50, ny = 110,

    # observation / implementation error
    obserr = 0.2, hobserr = 0.1, harverr = 0.15,

    # regime settings
    reglength = 20, regstr = 1,

    # management defaults
    harvrate = 0.5, goalfreq = 10,

    # stock–recruit process (profile-specific)
    procerr = 0.38,
    rho = 0.593,
    alpha_mean = 6.463,
    beta_mean  = 1.039e-5,
    sr_corr    = 0.0,
    sr_parms_sd= 0.2,

    # ages / structure
    ages = 1:9,
    meanageini = 5.5,
    sdage = 0.6,

    # trend baselines (multiplied by scenario switches)
    agetrend   = -0.4,
    sizetrends = c(0, 0, 0, 0, 10, -40, rep(-100, 3)),
    propF      = 0.41,
    propFtrend = -0.275,
    agediff    = 1,

    # growth / size
    vonB_Linf = 1200,
    vonB_k    = 0.325,
    ocean0s   = 100,
    sdSaA     = 0.01,

    # allometry
    allometry = c(
      size_ref = 800,
      fec_ref  = 6600,
      b_fec    = 2.4,
      egg_ref  = 916,
      b_eggs   = 4.8
    ),

    # alternative SR parameters
    alt_sr_param = c(
      alpha_N0   = 5.07,
      alpha_EASL = 0.00177,
      alpha_EMASL= 0.01036,
      beta_N0    = 8.6e-6,
      beta_EASL  = 2.973e-09,
      beta_EMASL = 1.6967e-08
    )
  )
)
