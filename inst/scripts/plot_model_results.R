# Quick diagnostic plots of model simulation output.
# Run from project root: source("tests/plot_model_results.R")
# Requires the package to be loaded: pkgload::load_all(".")

library(ggplot2)
library(dplyr)
library(tidyr)

out_dir <- "tests"

# ── 1. Run a small batch of iterations ──────────────────────────────────────

build_cfg <- function(scen_row, seed, harvmgmt_override = NULL) {
  profile <- param_configs[["Kuskokwim"]]
  cfg <- utils::modifyList(SizeShiftsImplicationsV2:::default_config(), profile)
  cfg <- SizeShiftsImplicationsV2:::apply_scenario_overrides(cfg, scen_row)
  if (!is.null(harvmgmt_override)) cfg$harvmgmt <- harvmgmt_override
  cfg$j <- 1L; cfg$k <- 1L
  cfg$seednum  <- as.integer(seed)
  cfg$log_dir  <- tempdir()
  SizeShiftsImplicationsV2:::validate_config(cfg)
}

scens    <- SizeShiftsImplicationsV2:::.select_scenarios("kuskokwim")
base_row <- scens[scens$scen_num == 13, ]   # no trends, smsy_goal
trend_row <- scens[scens$scen_num == 15, ]  # age-sex-length trends, smsy_goal

niter  <- 30
seeds  <- seq_len(niter)

message("Running ", niter, " iterations for base scenario (no trends)...")
runs_base <- lapply(seeds, function(s) {
  run_model(build_cfg(base_row, s))
})

message("Running ", niter, " iterations for trend scenario...")
runs_trend <- lapply(seeds, function(s) {
  run_model(build_cfg(trend_row, s))
})

# ── Helper: extract annual data from a list of runs ─────────────────────────
collect_data <- function(run_list, scenario_label) {
  purrr::map_dfr(seq_along(run_list), function(i) {
    d <- run_list[[i]]$data
    d$iter <- i
    d$scenario <- scenario_label
    d
  })
}

df_base  <- collect_data(runs_base,  "No trends")
df_trend <- collect_data(runs_trend, "Age-sex-length trends")
df_all   <- bind_rows(df_base, df_trend)
df_all$year <- df_all$Year  # alias

# ── 2. Plot A: population dynamics time series (ribbons + mean) ──────────────

dyn_summary <- df_all |>
  group_by(scenario, Year) |>
  summarise(
    Esc_med  = median(Esc,  na.rm = TRUE),
    Esc_lo   = quantile(Esc,  0.1, na.rm = TRUE),
    Esc_hi   = quantile(Esc,  0.9, na.rm = TRUE),
    Harv_med = median(Harv, na.rm = TRUE),
    Harv_lo  = quantile(Harv, 0.1, na.rm = TRUE),
    Harv_hi  = quantile(Harv, 0.9, na.rm = TRUE),
    Rec_med  = median(Rec,  na.rm = TRUE),
    Rec_lo   = quantile(Rec,  0.1, na.rm = TRUE),
    Rec_hi   = quantile(Rec,  0.9, na.rm = TRUE),
    .groups = "drop"
  )

dyn_long <- dyn_summary |>
  pivot_longer(
    cols = -c(scenario, Year),
    names_to = c("metric", "stat"),
    names_pattern = "(.+)_(med|lo|hi)"
  ) |>
  pivot_wider(names_from = stat, values_from = value) |>
  mutate(metric = factor(metric, levels = c("Rec","Esc","Harv"),
                         labels = c("Recruitment","Escapement","Harvest")))

pA <- ggplot(dyn_long, aes(x = Year, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(aes(y = med), linewidth = 0.8) +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c("No trends" = "#2166ac",
                                 "Age-sex-length trends" = "#d6604d"),
                      name = NULL) +
  scale_fill_manual(values  = c("No trends" = "#2166ac",
                                "Age-sex-length trends" = "#d6604d"),
                    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Population dynamics: median ± 80% interval",
       subtitle = sprintf("Kuskokwim profile | smsy_goal | %d iterations", niter),
       x = "Simulation year", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "plot_dynamics.png"), pA, width = 8, height = 7, dpi = 150)
message("Saved plot_dynamics.png")


# ── 3. Plot B: stock-recruit scatter ─────────────────────────────────────────

sr_df <- df_all |>
  filter(!is.na(Esc), !is.na(Rec), Esc > 0, Rec > 0) |>
  mutate(lnRS = log(Rec / Esc))

# Ricker curve from median alpha/beta across runs
fit_ricker <- function(run_list) {
  alphas <- sapply(run_list, function(r) r$para["alpha.low"])
  betas  <- sapply(run_list, function(r) r$para["beta"])
  list(alpha = median(alphas), beta = median(betas))
}

rk_base  <- fit_ricker(runs_base)
rk_trend <- fit_ricker(runs_trend)

s_seq <- seq(0, max(sr_df$Esc, na.rm = TRUE) * 1.05, length.out = 300)
curve_df <- bind_rows(
  data.frame(Esc = s_seq, scenario = "No trends",
             Rec = rk_base$alpha  * s_seq * exp(-rk_base$beta  * s_seq)),
  data.frame(Esc = s_seq, scenario = "Age-sex-length trends",
             Rec = rk_trend$alpha * s_seq * exp(-rk_trend$beta * s_seq))
)

pB <- ggplot(sr_df, aes(x = Esc, y = Rec, colour = scenario)) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_line(data = curve_df, linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("No trends" = "#2166ac",
                                 "Age-sex-length trends" = "#d6604d"),
                      name = NULL) +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Stock-recruit relationship",
       subtitle = "Points = simulated years; curves = median Ricker fit",
       x = "Escapement (S)", y = "Recruitment (R)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "plot_stock_recruit.png"), pB, width = 7, height = 5, dpi = 150)
message("Saved plot_stock_recruit.png")


# ── 4. Plot C: size-at-age trends ───────────────────────────────────────────

ages_keep <- c(3, 4, 5, 6, 7)

# meanSaA loses dimnames after matrix arithmetic; rebuild manually
saa_df <- purrr::map_dfr(seq_along(runs_trend), function(i) {
  mat      <- runs_trend[[i]]$meanSaA   # (nyr x nage); columns may be V1..V9
  age_vals <- seq_len(ncol(mat))        # 1:9
  data.frame(
    Year = rep(seq_len(nrow(mat)), times = ncol(mat)),
    age  = rep(age_vals,           each  = nrow(mat)),
    size = as.vector(mat),
    iter = i
  )
}) |>
  filter(age %in% ages_keep, !is.na(size), size > 0)

saa_sum <- saa_df |>
  group_by(Year, age) |>
  summarise(size_med = median(size, na.rm = TRUE),
            size_lo  = quantile(size, 0.1, na.rm = TRUE),
            size_hi  = quantile(size, 0.9, na.rm = TRUE),
            .groups  = "drop")

pC <- ggplot(saa_sum, aes(x = Year, y = size_med,
                          colour = factor(age), fill = factor(age))) +
  geom_ribbon(aes(ymin = size_lo, ymax = size_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_brewer(palette = "Dark2", name = "Age") +
  scale_fill_brewer(palette  = "Dark2", name = "Age") +
  labs(title    = "Mean size-at-age over time (trend scenario)",
       subtitle = "Median ± 80% interval across iterations",
       x = "Simulation year", y = "Mean length (mm)") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "plot_size_at_age.png"), pC, width = 8, height = 5, dpi = 150)
message("Saved plot_size_at_age.png")


# ── 5. Plot D: regime shift effect on recruitment (Bug 1 was here) ──────────

message("Running regime-shift comparison runs...")
shift_row <- base_row
shift_row$sizeT <- 0; shift_row$ageT <- 0; shift_row$sexT <- 0  # no trends, shifts only

runs_noshift <- lapply(seeds, function(s) {
  cfg <- build_cfg(shift_row, s)
  cfg$reglength <- 0; cfg$regstr <- 1
  run_model(SizeShiftsImplicationsV2:::validate_config(cfg))
})

runs_shift <- lapply(seeds, function(s) {
  cfg <- build_cfg(shift_row, s)
  cfg$reglength <- 20; cfg$regstr <- 2  # default Ohlberger regime settings
  run_model(SizeShiftsImplicationsV2:::validate_config(cfg))
})

reg_df <- bind_rows(
  collect_data(runs_noshift, "No regime shifts"),
  collect_data(runs_shift,   "Regime shifts (regstr=2)")
)

reg_sum <- reg_df |>
  group_by(scenario, Year) |>
  summarise(
    Rec_med = median(Rec, na.rm = TRUE),
    Rec_lo  = quantile(Rec, 0.1, na.rm = TRUE),
    Rec_hi  = quantile(Rec, 0.9, na.rm = TRUE),
    .groups = "drop"
  )

pD <- ggplot(reg_sum, aes(x = Year, colour = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = Rec_lo, ymax = Rec_hi), alpha = 0.15, colour = NA) +
  geom_line(aes(y = Rec_med), linewidth = 0.9) +
  scale_colour_manual(values = c("No regime shifts"       = "#1b7837",
                                 "Regime shifts (regstr=2)" = "#762a83"),
                      name = NULL) +
  scale_fill_manual(values  = c("No regime shifts"       = "#1b7837",
                                "Regime shifts (regstr=2)" = "#762a83"),
                    name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(title    = "Effect of regime shifts on recruitment",
       subtitle = "Median ± 80% interval | Bug 1 fixed: alpha changes in the shift year",
       x = "Simulation year", y = "Recruitment") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "plot_regime_shifts.png"), pD, width = 8, height = 5, dpi = 150)
message("Saved plot_regime_shifts.png")

message("\nAll plots saved to: ", out_dir)
