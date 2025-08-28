# Silence R CMD check “no visible binding for global variable ...” for NSE columns
utils::globalVariables(c(
  ".data", ":=", "value", "metric_label",
  "trends", "factorMSY", "mgmt",
  "escapement", "harvest",
  "nscen", "niter", "nyh", "obs.list",
  "n_zero_harvest", "n_obs", "recRec", "obsEsc",
  "scen_num", "scenarios",
  "alpha_N0", "beta_N0", "alpha_EASL", "beta_EASL", "alpha_EMASL", "beta_EMASL",
  "figure_dir"
))
