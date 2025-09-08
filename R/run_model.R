#' Run Management Strategy Evaluation (MSE) simulation
#'
#' This internal function runs a management strategy evaluation simulation
#' using a specified parameter configuration and scenario. It generates
#' escapement, harvest, recruitment, age composition, and size-at-age outputs.
#'
#' @param config List. List of parameter values
#'
#' @return A named list.
#' @keywords internal

run_model <- function(config) {

  ####################################### unpack configuration parameters
  agediff      <- config$agediff
  agetrend     <- config$agetrend     ## mean age trend (mm/year)
  ages         <- config$ages         ## ages to simulate
  allometry    <- config$allometry    ## allometry parameters
  alpha_mean   <- config$alpha_mean   ## mean alpha parameter
  alt_sr_param <- config$alt_sr_param ## alternative SR parameters
  beta_mean    <- config$beta_mean    ## mean beta parameter
  factorMSY    <- config$factorMSY    ## MSY scaling factor
  futureT      <- config$futureT      ## future years after historical period
  goalfreq     <- config$goalfreq     ## frequency of MSY goal reviews
  harverr      <- config$harverr      ## harvest implementation error
  harvmgmt     <- config$harvmgmt     ## harvest management type
  harvrate     <- config$harvrate     ## harvest rate
  hobserr      <- config$hobserr      ## harvest observation error
  maxsel       <- config$maxsel       ## maximum selectivity
  meanageini   <- config$meanageini   ## initial mean age (mm)
  ny           <- config$ny           ## years to simulate
  nyh          <- config$nyh          ## historical years
  nyi          <- config$nyi          ## initial years
  obserr       <- config$obserr       ## escapement observation error
  ocean0s      <- config$ocean0s      ## initial ocean sizes
  procerr      <- config$procerr      ## process error
  propF        <- config$propF
  propFtrend   <- config$propFtrend
  reglength    <- config$reglength    ## length of regime shift
  regstr       <- config$regstr       ## regime shift strength
  rho          <- config$rho          ## AR(1) correlation
  ricker_type  <- config$ricker_type  ## type of Ricker model
  sdSaA        <- config$sdSaA        ## SD of size-at-age anomalies
  sdage        <- config$sdage        ## SD of age
  sdsel        <- config$sdsel        ## SD of selectivity
  seednum      <- config$seednum
  sim_recruits <- config$sim_recruits ## recruits generation mode
  sizetrends   <- config$sizetrends   ## size-at-age trends
  sr_corr      <- config$sr_corr      ## corr(alpha, beta)
  sr_parms_sd  <- config$sr_parms_sd  ## SD of SR parameters
  var          <- config$var          ## DLM variable to vary
  vonB_k       <- config$vonB_k       ## von Bertalanffy k
  vonB_Linf    <- config$vonB_Linf    ## von Bertalanffy Linf

  ## NOTE: no global set.seed() here; we use withr::with_seed(seednum + OFFSET)

  #################################### generate stock-recruit parameters
  ##====================================## Ricker parameters and regimes
  ##------------------------------------------------------------## alpha
  alpha.low <- withr::with_seed(seednum + 101L,
                                rlnorm(1, meanlog = log(alpha_mean), sdlog = sr_parms_sd)
  )
  alpha.high <- alpha.low * regstr ## high productivity regime
  aset <- 1 ## initial regime is assumed to be low productivity
  if (reglength == 0) aset <- 1 ## also without regime shifts
  alpha <- c(alpha.low, alpha.high)[aset] ## select alpha parameter

  ##--------------------------------------------------------## beta/rmax
  if (sr_corr != 0) { ## in case of correlated parameters
    beta_corr <- beta_mean * (1 + (alpha - alpha_mean) * sr_corr)
    beta <- withr::with_seed(seednum + 102L,
                             rlnorm(1, meanlog = log(beta_corr), sdlog = sr_parms_sd)
    )
  } else {
    beta <- withr::with_seed(seednum + 102L,
                             rlnorm(1, meanlog = log(beta_mean), sdlog = sr_parms_sd)
    )
  }
  maxr <- (alpha / beta) * exp(-1) ## Hilborn & Walters 1992 (max recruitment)
  srparms <- c(alpha.low = alpha.low, alpha.high = alpha.high, beta = beta, maxr = maxr)

  ##------------------------------------------------------## S_msy/U_msy (true)
  S_msy_true <- round((1 - gsl::lambert_W0(exp(1 - log(alpha)))) / beta)
  U_msy_true <- round((1 - gsl::lambert_W0(exp(1 - log(alpha)))), 3)

  ##===========================================## alternative parameters
  ## parameters for fecundity/eggmass models (Staton et al. 2021)
  for (i in 1:6) { assign(names(alt_sr_param)[i], as.numeric(alt_sr_param[i])) }
  ##------------------------------------------------------------## Rmax
  rmax_N0   <- (alpha_N0  / beta_N0)  * exp(-1)
  rmax_EASL <- (alpha_EASL / beta_EASL) * exp(-1)
  rmax_EMASL<- (alpha_EMASL/ beta_EMASL)* exp(-1)
  ##-------------------------------------------------------## fecundity/eggmass
  if (sim_recruits == "fecundity") {
    alpha_fec <- (alpha_EASL / alpha_N0) * alpha
    beta_fec  <- (beta_EASL  / beta_N0)  * beta
    rmax_fec  <- (rmax_EASL  / rmax_N0)  * maxr
  }
  if (sim_recruits == "eggmass") {
    alpha_egg <- (alpha_EMASL / alpha_N0) * alpha
    beta_egg  <- (beta_EMASL  / beta_N0)  * beta
    rmax_egg  <- (rmax_EMASL  / rmax_N0)  * maxr
  }

  #################################### arrays for storing simulated data
  nage <- length(ages)
  sexes <- c("M","F"); nsex <- length(sexes)
  nyr <- nyi + ny + nage + 1 ## additional years dropped during reconstruction
  Yvec <- c(seq(1, nyr))
  ##---------------------------------------------## true population data
  pnames <- c("Year","Ret","Rec","eps","reg","Esc","Harv","RepOut")
  PopDat <- data.frame(array(NA, dim = c(nyr, length(pnames))))
  dimnames(PopDat) <- list(Yvec, c(pnames))
  PopDat$Year <- seq(nyr)
  ##-------------------------------## calendar year age-sex compositions
  ret_by_age_sex <- array(NA, dim = c(nyr, nage, nsex))
  dimnames(ret_by_age_sex) <- list(Yvec, paste0("Ret", ages), sexes)
  harv_by_age_sex <- array(NA, dim = c(nyr, nage, nsex))
  dimnames(harv_by_age_sex) <- list(Yvec, paste0("Harv", ages), sexes)
  esc_by_age_sex <- array(NA, dim = c(nyr, nage, nsex))
  dimnames(esc_by_age_sex) <- list(Yvec, paste0("Esc", ages), sexes)
  ##-----------------------------------------------------## observations
  onames <- c("obsRet","obsEsc","obsHarv","RepOut")
  ObsDat <- data.frame(array(NA, dim = c(nyr, length(onames))))
  dimnames(ObsDat) <- list(c(seq(1, nyr)), onames)
  ##----------------------------------------------## harvest selectivity
  selectivities_by_age <- array(0, dim = c(nyr, nage)) ## selectivity by age

  ################################################### demographic trends
  max_yh <- nyi + nyh          ## last year of historical period
  yrs_h  <- (nyi + 1):max_yh   ## historical years index
  yrs_f  <- (nyi + nyh + 1):nyr## future years index

  ##=========================================================## mean age
  ## deterministic trend but age proportions drawn based on mean by year
  meanage <- NA
  meanage[1:nyi] <- meanageini
  for (y in yrs_h) { meanage[y] <- meanage[y-1] + agetrend / (nyh) }
  if (futureT == "yes") {
    for (y in yrs_f) meanage[y] <- meanage[y-1] + agetrend / (nyh)
  } else {
    meanage[yrs_f] <- meanage[max_yh]
  }

  ##===============================================## proportion female
  ## propF trend given as change over entire period in logit space
  propF_y <- NA; propF_y[1:nyi] <- stats::qlogis(propF)
  n_ref <- 5
  for (y in yrs_h) { propF_y[y] <- propF_y[y-1] + propFtrend / (nyh) }
  if (futureT == "yes") {
    for (y in yrs_f) propF_y[y] <- propF_y[y-1] + propFtrend / (nyh)
  } else {
    propF_y[yrs_f] <- mean(propF_y[(max_yh - n_ref + 1):max_yh])
  }
  propF_y <- propF_y + withr::with_seed(seednum + 201L, rnorm(nyr, 0, 0.05)) ## normal error
  propfemale <- stats::plogis(propF_y)

  ##=================================================## mean size-at-age
  meanSaA <- meanSaAdet <- SaA_anoms <- array(0, dim = c(nyr, nage))
  dimnames(meanSaA) <- list(Yvec, paste0("age", ages))
  ##----------------------------------------## initial mean sizes at age
  oceanage <- ages[1:(length(ages) - 2)] ## ocean ages
  growth <- function(age, Linf, k) { Linf * (1 - exp(-k * age)) } ## growth fn
  meansini <- round((ocean0s + growth(oceanage, Linf = (vonB_Linf - ocean0s), k = vonB_k)) / 10) * 10
  meanSaA_ini <- c(0, ocean0s, meansini) ## fish enter the ocean at age 2 (BY+2)
  ##--------------------------------------## size-at-age trend (mm/year)
  trendsSaA <- sizetrends / (nyh) ## for reconstructed time period
  ##--------------------------------------## deterministic size at age
  for (y in 1:nyi) { meanSaAdet[y,] <- meanSaA_ini } ## same for initial years
  for (y in yrs_h) { for (a in 1:nage) { meanSaAdet[y, a] <- meanSaAdet[y-1, a] + trendsSaA[a] } }
  if (futureT == "yes") {
    for (y in yrs_f) { for (a in 1:nage) { meanSaAdet[y, a] <- meanSaAdet[y-1, a] + trendsSaA[a] } }
  } else {
    for (a in 1:nage) { meanSaAdet[yrs_f, a] <- mean(meanSaAdet[(max_yh - n_ref + 1):max_yh, a]) }
  }
  ##-------------------------------------------## size-at-age anomalies
  for (a in 1:nage) {
    SaA_anoms[, a] <- withr::with_seed(seednum + 300L + a,
                                       rnorm(nyr, mean = 0, sd = meanSaA_ini[a] * sdSaA)
    )
  }
  ##-------------------------------------------## stochastic size at age
  meanSaA <- meanSaAdet + SaA_anoms

  ##==================================================================##
  ##=========================================## generate population data
  ##==================================================================##
  age_comp <- array(0, dim = c(nyr, nage, nsex))

  ##################################### generate first few years of data
  PopDat$Esc[1:nyi] <- round(exp(withr::with_seed(seednum + 400L,
                                                  rnorm(nyi, log(S_msy_true) - 0.5 * 0.25^2, 0.25)
  )))
  for (y in 1:nyi) {
    PopDat$reg[y] <- aset ## regime (no regime shifts yet)

    ##----------------------------------## generate recruits by brood year
    if (sim_recruits == "spawners") {
      sr <- .calc_ricker(spawn = PopDat$Esc[y], sigma = procerr,
                        alpha = alpha, beta = beta, rho = rho,
                        last.eps = ifelse(y == 1, 0, PopDat$eps[y-1]),
                        seed = seednum + 500L + y
      )
      PopDat$RepOut[y] <- PopDat$Esc[y]
    }

    ##----------------------------------## alternative reproductive units
    if (sim_recruits != "spawners") {
      age_comp_y <- .calc_agecomp(
        ages = ages, recruits = round(PopDat$Esc[y] * (1 - propfemale[y])),
        meanage = meanage[y] - agediff / 2, sdage = sdage,
        seed = seednum + 515L + y
      )
      sizes_y <- rep(meanSaA[y,], prop.table(age_comp_y) * PopDat$Esc[y])

      if (sim_recruits == "fecundity") {
        fecundity_total <- round(sum(exp(.calc_reprod_output(sizes_y, allometry)[[1]])))
        if (ricker_type == "const_beta") {
          sr <- .calc_ricker(
                            spawn = fecundity_total, sigma = procerr,
                            alpha = alpha_fec, beta = beta_fec, rho = rho,
                            last.eps = ifelse(y == 1, 0, PopDat$eps[y-1]),
                            seed = seednum + 520L + y
          )
        }
        PopDat$RepOut[y] <- fecundity_total
      }

      if (sim_recruits == "eggmass") {
        eggmass_total <- round(sum(exp(.calc_reprod_output(sizes_y, allometry)[[2]])))
        if (ricker_type == "const_beta") {
          sr <- .calc_ricker(
                            spawn = eggmass_total, sigma = procerr,
                            alpha = alpha_egg, beta = beta_egg, rho = rho,
                            last.eps = ifelse(y == 1, 0, PopDat$eps[y-1]),
                            seed = seednum + 521L + y
          )
        }
        PopDat$RepOut[y] <- eggmass_total
      }
    } ## end if 'sim_recruits' statement

    ##---------------------------------## store recruits from brood year y
    PopDat$Rec[y] <- round(sr$rec)
    ##-------------------## store anomaly (epsilon for AR1 process year y)
    PopDat$eps[y] <- round(sr$eps, 5)

    ##-------------------## generate sex-specific recruit age compositions
    age_comp[y, , 1] <- .calc_agecomp(
                                     ages, round(PopDat$Rec[y] * (1 - propfemale[y])),
                                     meanage[y] - agediff/2, sdage,
                                     seed = seednum + 530L + y
    )
    age_comp[y, , 2] <- .calc_agecomp(
                                     ages, round(PopDat$Rec[y] * propfemale[y]),
                                     meanage[y] + agediff/2, sdage,
                                     seed = seednum + 531L + y
    )
    for (k in 1:nage) {
      ret_by_age_sex[y + k, k, 1] <- age_comp[y, k, 1]
      ret_by_age_sex[y + k, k, 2] <- age_comp[y, k, 2]
    }

    ##-----------------------## generate escapement observations from data
    ObsDat$obsEsc[y] <- PopDat$Esc[y] * exp(withr::with_seed(seednum + 540L + y,
                                                             rnorm(1, -0.5 * obserr^2, obserr)
    ))
  } ## end loop over initial years

  ##################################### years for re-assessing MSY goals
  if (harvmgmt == "fix_harv_rate") goalfreq <- FALSE
  firstrev <- 20
  if (goalfreq) { goalrev <- seq(nyi + firstrev, ny, goalfreq) } else { goalrev <- ny }
  nrev <- length(goalrev)

  ##-------## set initial MSY goal (esc or U) as mean/true values
  if (harvmgmt %in% c("smsy_goal","s_eq_goal","smsy_dlm_goal")) {
    msygoal <- S_msy_true; msygoalini <- msygoal
  }
  if (harvmgmt %in% c("umsy_goal","u_eq_goal","umsy_dlm_goal")) {
    msygoal <- U_msy_true; msygoalini <- msygoal
  }

  ##-------------------------------------## storing errors, bias, goals
  impl_errors <- MSY_Goals <- HarvRates <- S_msy_estimate <- NA

  ######################################## loop through esc goal reviews
  for (i in 1:nrev) {
    ##=========================================## calculate years of data
    y_rev <- goalrev[i]
    if (i == 1) nycnt <- nyi else nycnt <- nyi + goalrev[i - 1]
    yindex <- (nycnt + 1):(nyi + y_rev)

    ##===============================## generate subsequent years of data
    MSY_Goals[yindex] <- msygoal
    for (y in yindex) {
      ## if regime shift switch productivity and change regime indicator
      if (reglength != 0) {
        if (withr::with_seed(seednum + 6000L + y,
                             sample(c(rep(0, reglength - 1), 1L), 1L)
        )) {
          if (aset == 1) alpha <- alpha.low else alpha <- alpha.high
          aset <- abs(aset - 3) ## change regime indicator
        }
      }
      PopDat$reg[y] <- aset ## regime
      ##===============================## total return across ages and sexes
      PopDat$Ret[y] <- sum(ret_by_age_sex[y,,], na.rm = TRUE)

      ##==============================## management using fixed harvest rate
      if (harvmgmt == "fix_harv_rate") {
        impl_err <- withr::with_seed(seednum + 6100L + y, rnorm(1, mean = 0, sd = harverr))
        harv_thresh <- round(maxr * 0.25) ## 25% max recruits as conservation threshold
        if (PopDat$Ret[y] >= harv_thresh) {
          harvrate_det <- harvrate
        } else {
          harvrate_det <- harvrate * (PopDat$Ret[y] / (PopDat$Ret[y] + round(maxr * 0.1)))
        }
        harvrate_real <- stats::plogis(stats::qlogis(harvrate_det) + impl_err)
        if (is.na(harvrate_real)) harvrate_real <- 0
        HarvRates[y] <- harvrate_real
        PopDat$Harv[y] <- round(PopDat$Ret[y] * harvrate_real) ## total harvest
        PopDat$Esc[y]  <- PopDat$Ret[y] - PopDat$Harv[y]       ## escapement
      }

      ##==================================## management by harvest rate goal
      if (harvmgmt %in% c("umsy_goal","u_eq_goal","umsy_dlm_goal")) {
        impl_err <- withr::with_seed(seednum + 6200L + y, rnorm(1, mean = 0, sd = harverr))
        harvrate_real <- stats::plogis(stats::qlogis(msygoal) + impl_err) ## realized harvest rate
        PopDat$Harv[y] <- round(PopDat$Ret[y] * harvrate_real)
        PopDat$Esc[y]  <- PopDat$Ret[y] - PopDat$Harv[y]
      }

      ##===================================## management by escapement goal
      if (harvmgmt %in% c("smsy_goal","s_eq_goal","smsy_dlm_goal")) {
        impl_err <- withr::with_seed(seednum + 6300L + y, rnorm(1, mean = -0.5 * harverr^2, sd = harverr))
        bias <- 0
        escape <- round(msygoal * exp(impl_err + bias))
        PopDat$Esc[y]  <- min(escape, PopDat$Ret[y]) ## no larger than return
        PopDat$Harv[y] <- PopDat$Ret[y] - PopDat$Esc[y]
      }

      ##===========================================## harvest by age and sex
      ## proportions pre-harvest
      ret_prop_by_sex_y <- prop.table(ret_by_age_sex[y,,], 1)
      ret_by_age_all_y  <- rowSums(ret_by_age_sex[y,,]) ## this year

      ## 'sample' approximation
      fac <- 10
      nsample <- round(PopDat$Harv[y] / fac)
      return_by_age <- round(ret_by_age_all_y / fac)
      return_by_age[is.na(return_by_age)] <- 0
      age_vector <- rep(ages, return_by_age)

      ## apply selectivity
      selectivities_by_age[y, ] <- withr::with_seed(seednum + 6400L + y,
                                                    .calc_selectivity(meanSaA[y, ], maxsel, sdsel)
      )
      select <- rep(selectivities_by_age[y, ], return_by_age)

      ## sample harvest without replacement
      harvested <- withr::with_seed(seednum + 6500L + y,
                                    sample(age_vector, size = nsample, prob = select, replace = FALSE)
      )
      tab <- table(harvested)

      ## redistribute harvest to sexes using pre-harvest proportions
      harv_prop_M <- ret_prop_by_sex_y[as.numeric(names(tab)), 1]
      harv_by_age_sex[y, as.numeric(names(tab)), 1] <- round(tab * fac * harv_prop_M)
      harv_prop_F <- ret_prop_by_sex_y[as.numeric(names(tab)), 2]
      harv_by_age_sex[y, as.numeric(names(tab)), 2] <- round(tab * fac * harv_prop_F)
      harv_by_age_sex[is.na(harv_by_age_sex)] <- 0

      ## avoid over-harvesting of age-by-sex groups (approximation fix)
      index <- which(harv_by_age_sex[y,,1] > ret_by_age_sex[y,,1])
      if (length(index) >= 1) { harv_by_age_sex[y, index, 1] <- ret_by_age_sex[y, index, 1] }
      index <- which(harv_by_age_sex[y,,2] > ret_by_age_sex[y,,2])
      if (length(index) >= 1) { harv_by_age_sex[y, index, 2] <- ret_by_age_sex[y, index, 2] }

      ##=======================================## escapement by age and size
      esc_by_age_sex[y,,] <- ret_by_age_sex[y,,] - harv_by_age_sex[y,,]

      ##================================================## generate recruits
      if (sim_recruits == "spawners") {
        sr <- .calc_ricker(
                          spawn = PopDat$Esc[y], sigma = procerr,
                          alpha = alpha, beta = beta, rho = rho,
                          last.eps = ifelse(y == 1, 0, PopDat$eps[y-1]),
                          seed = seednum + 6600L + y
        )
        PopDat$RepOut[y] <- PopDat$Esc[y]
      }
      if (sim_recruits != "spawners") {
        ## all individual sizes for females in escapement this year
        sizes_y <- NA
        for (a in 1:nage) {
          esc_F_a_y <- esc_by_age_sex[y, a, which(names(esc_by_age_sex[1,1,]) == "F")]
          if (is.na(esc_F_a_y)) { size_new <- NA } else { if (esc_F_a_y > 0) { size_new <- meanSaA[y, a] } else { size_new <- NA } }
          sizes_y <- c(sizes_y, rep(size_new, esc_F_a_y))
        }
        sizes_y <- sizes_y[!is.na(sizes_y)]

        if (sim_recruits == "fecundity") {
          fecundity_total <- round(sum(exp(.calc_reprod_output(sizes_y, allometry)[[1]])))
          if (ricker_type == "const_beta") {
            sr <- .calc_ricker(
                              spawn = fecundity_total, sigma = procerr,
                              alpha = alpha_fec, beta = beta_fec, rho = rho,
                              last.eps = ifelse(y == 1, 0, PopDat$eps[y-1]),
                              seed = seednum + 6610L + y
            )
          }
          PopDat$RepOut[y] <- fecundity_total
        }

        if (sim_recruits == "eggmass") {
          eggmass_total <- round(sum(exp(.calc_reprod_output(sizes_y, allometry)[[2]])))
          if (ricker_type == "const_beta") {
            sr <- .calc_ricker(
                              spawn = eggmass_total, sigma = procerr,
                              alpha = alpha_egg, beta = beta_egg, rho = rho,
                              last.eps = ifelse(y == 1, 0, PopDat$eps[y-1]),
                              seed = seednum + 6620L + y
            )
          }
          PopDat$RepOut[y] <- eggmass_total
        }
      } ## end if 'sim_recruits' statement

      ##---------------------------------## store recruits from brood year y
      PopDat$Rec[y] <- round(sr$rec)
      ##-------------------## store anomaly (epsilon for AR1 process year y)
      PopDat$eps[y] <- round(sr$eps, 5)

      ##--------------------## generate sex-specific recruit age composition
      age_comp[y, , 1] <- .calc_agecomp(
                                       ages, round(PopDat$Rec[y] * (1 - propfemale[y])),
                                       meanage[y] - agediff/2, sdage,
                                       seed = seednum + 6700L + y
      )
      age_comp[y, , 2] <- .calc_agecomp(
                                       ages, round(PopDat$Rec[y] * propfemale[y]),
                                       meanage[y] + agediff/2, sdage,
                                       seed = seednum + 6710L + y
      )
      for (k in 1:nage) {
        ret_by_age_sex[y + k, k, 1] <- age_comp[y, k, 1]
        ret_by_age_sex[y + k, k, 2] <- age_comp[y, k, 2]
      }

      ##=================================## generate observations from data
      ObsDat$obsEsc[y] <- PopDat$Esc[y] * exp(withr::with_seed(seednum + 6800L + y,
                                                               rnorm(1, -0.5 * obserr^2, obserr)
      ))
      ObsDat$obsHarv[y] <- PopDat$Harv[y] * exp(withr::with_seed(seednum + 6810L + y,
                                                                 rnorm(1, -0.5 * hobserr^2, hobserr)
      ))
      ## store observed return (harvest + escapement)
      ObsDat$obsRet[y] <- ObsDat$obsEsc[y] + ObsDat$obsHarv[y]
      ## reproductive output observation (deterministic bookkeeping)
      ObsDat$RepOut[y] <- PopDat$RepOut[y]
    } ## end year loop

    ##==========================## round population numbers to individuals
    roundit <- c("Ret","Rec","Esc","Harv","RepOut")
    PopDat[, names(PopDat) %in% roundit] <- round(PopDat[, names(PopDat) %in% roundit])
    ObsDat <- round(ObsDat)

    ##========================================## simulated population data
    data_all <- PopDat
    data_all[is.na(data_all)] <- 0
    year_index <- (nyi + 1):(nyi + y_rev - nage - 1) ## match reconstructed data
    nyrec <- length(year_index) ## number of reconstructed years
    data <- data.frame(data_all[year_index, ])

    ##################################################### assess MSY goals
    ##=====================================## Ricker fit to simulated data
    data$lnRS <- log(data$Rec / data$Esc)
    ##--------------------------------------------------------------## GLS
    mod_sim <- nlme::gls(lnRS ~ Esc, data = data, correlation = nlme::corAR1())
    sig_sim <- summary(mod_sim)$sigma ## sigma to correct for log-normal error
    phi_sim <- coef(mod_sim$modelStruct$corStruct, unconstrained = FALSE)
    log_a_sim <- summary(mod_sim)$coefficients[1] ## log-productivity
    a_sim <- exp(log_a_sim) ## productivity (alpha)
    b_sim <- -summary(mod_sim)$coefficients[2] ## density-dependence (beta)
    alpha_sim <- exp(log_a_sim + ((0.5 * sig_sim^2) / (1 - phi_sim^2))) ## corrected
    beta_sim  <- b_sim * (alpha_sim / a_sim) ## corrected beta value
    S_msy_sim <- (1 - gsl::lambert_W0(exp(1 - log(alpha_sim)))) / beta_sim
    U_msy_sim <- (1 - gsl::lambert_W0(exp(1 - log(alpha_sim))))
    S_max_sim <- 1 / beta_sim
    ##----------------------------------------------------------## results
    sr_sim <- data.frame(rbind(c(alpha_sim, beta_sim, phi_sim, sig_sim)))
    colnames(sr_sim) <- c("alpha","beta","rho","sigma")

    ##===================================================## reconstruction
    ## reconstruct recruitment from year y using observed escapement and harvest
    dataObs <- ObsDat ## all years simulated
    dataObs[is.na(dataObs)] <- 0
    ##--------------------------------------## escapement age proportions
    age_comp_esc <- rowSums(esc_by_age_sex, dims = 2)
    age_comp_esc[is.na(age_comp_esc)] <- 0
    age_prop_esc <- prop.table(as.matrix(age_comp_esc), margin = 1)
    age_comp_esc_obs <- array(0, dim = dim(age_prop_esc))
    for (y in 1:nyr) { age_comp_esc_obs[y, ] <- age_prop_esc[y, ] * dataObs$obsEsc[y] }
    age_comp_esc_obs[is.na(age_comp_esc_obs)] <- 0
    age_comp_esc_obs <- round(age_comp_esc_obs)
    ##------------------------------------------## harvest age proportions
    age_comp_harv <- rowSums(harv_by_age_sex, dims = 2)
    age_comp_harv[is.na(age_comp_harv)] <- 0
    age_prop_harv <- prop.table(as.matrix(age_comp_harv), margin = 1)
    age_comp_harv_obs <- array(0, dim = dim(age_prop_harv))
    for (y in 1:nyr) { age_comp_harv_obs[y, ] <- age_prop_harv[y, ] * dataObs$obsHarv[y] }
    age_comp_harv_obs[is.na(age_comp_harv_obs)] <- 0
    age_comp_harv_obs <- round(age_comp_harv_obs)

    ##----------------------------------------## recruitment by brood year
    nage <- length(ages)
    dataObs$recRec <- NA
    age_comp_ret_obs_byBY <- array(NA, dim = dim(age_comp_harv_obs))
    for (yy in (nyi + 1):(nyi + ny - nage)) {
      sum_year_ages <- 0
      for (a in 1:nage) {
        sum_broodyear <- age_comp_esc_obs[yy + a, a] + age_comp_harv_obs[yy + a, a]
        age_comp_ret_obs_byBY[yy, a] <- sum_broodyear
        sum_year_ages <- sum_year_ages + sum_broodyear
      }
      dataObs$recRec[yy] <- sum_year_ages
    }
    dataObs <- data.frame(dataObs[year_index, ])
    dataObs <- dataObs[!is.na(dataObs$recRec), ]
    dataObs <- dataObs[dataObs$recRec > 0, ]

    ##================================## Ricker fits to reconstructed data
    dataObs$obsRpS  <- dataObs$recRec / dataObs$obsEsc
    dataObs$obslnRS <- log(dataObs$recRec / dataObs$obsEsc)

    ##---------------------------------------------------------## MSY goals
    if (harvmgmt %in% c("smsy_goal","umsy_goal")) {
      mod_obs <- nlme::gls(obslnRS ~ obsEsc, data = dataObs, correlation = nlme::corAR1())
      sig_obs <- summary(mod_obs)$sigma
      phi_obs <- coef(mod_obs$modelStruct$corStruct, unconstrained = FALSE)
      dataObs$mod_obs_resid <- residuals(mod_obs)
      log_a_obs <- summary(mod_obs)$coefficients[1]
      a_obs <- exp(log_a_obs)
      b_obs <- -summary(mod_obs)$coefficients[2]
      alpha_obs <- exp(log_a_obs + ((0.5 * sig_obs^2) / (1 - phi_obs^2)))
      beta_obs  <- b_obs * (alpha_obs / a_obs)
      S_msy_obs <- (1 - gsl::lambert_W0(exp(1 - log(alpha_obs)))) / beta_obs
      U_msy_obs <- (1 - gsl::lambert_W0(exp(1 - log(alpha_obs))))
      S_max_obs <- 1 / beta_obs
      sr_obs <- data.frame(rbind(c(alpha_obs, beta_obs, phi_obs, sig_obs)))
      colnames(sr_obs) <- c("alpha","beta","rho","sigma")
    }

    ##==============================================================## DLM
    if (harvmgmt %in% c("smsy_dlm_goal","umsy_dlm_goal")) {
      dataDLM <- dplyr::select(dataObs, Rec = recRec, Esc = obsEsc)
      nrd <- dim(dataDLM)[1]
      if (var == "alpha") mod_dlm <- suppressWarnings(.calc_DLMfit(data = dataDLM, var_alpha = TRUE,  var_beta = FALSE))
      if (var == "beta")  mod_dlm <- suppressWarnings(.calc_DLMfit(data = dataDLM, var_alpha = FALSE, var_beta = TRUE))
      if (var == "both")  mod_dlm <- suppressWarnings(.calc_DLMfit(data = dataDLM, var_alpha = TRUE,  var_beta = TRUE))
      sig_dlm <- mod_dlm$sigma
      log_a_dlm <- mod_dlm$results$alpha_y
      a_dlm <- exp(log_a_dlm)
      b_dlm <- -mod_dlm$results$beta_y
      alpha_dlm <- signif(exp(log_a_dlm + ((0.5 * sig_dlm^2))), 4)
      beta_dlm  <- signif(b_dlm * (alpha_dlm / a_dlm), 4)
      S_msy_dlm <- round((1 - gsl::lambert_W0(exp(1 - log(alpha_dlm)))) / beta_dlm)
      U_msy_dlm <- round((1 - gsl::lambert_W0(exp(1 - log(alpha_dlm)))), 3)
      dlm_out <- data.frame(cbind(alpha_dlm, beta_dlm, S_msy_dlm, U_msy_dlm))
      ## recent years Smsy/Umsy
      ndlm <- 5
      alpha.dlm <- mean(alpha_dlm[(nrd - ndlm + 1):nrd])
      beta.dlm  <- mean(beta_dlm[(nrd - ndlm + 1):nrd])
      Smsy_dlm  <- round((1 - gsl::lambert_W0(exp(1 - log(alpha.dlm)))) / beta.dlm)
      Umsy_dlm  <- round((1 - gsl::lambert_W0(exp(1 - log(alpha.dlm)))), 3)
      Smax_dlm  <- 1 / beta.dlm
    }

    ########################################## yield-per-recruit analysis
    if (harvmgmt %in% c("s_eq_goal","u_eq_goal")) {
      dataObs$lnRperRU <- log(dataObs$recRec / dataObs$RepOut)
      mod_rep_out <- nlme::gls(lnRperRU ~ RepOut, data = dataObs, correlation = nlme::corAR1())
      sig_rep_out <- summary(mod_rep_out)$sigma
      phi_rep_out <- coef(mod_rep_out$modelStruct$corStruct, unconstrained = FALSE)
      log_a_rep_out <- summary(mod_rep_out)$coefficients[1]
      a_rep_out <- exp(log_a_rep_out)
      b_rep_out <- -summary(mod_rep_out)$coefficients[2]
      alpha_rep_out <- exp(log_a_rep_out + ((0.5 * sig_rep_out^2) / (1 - phi_rep_out^2)))
      beta_rep_out  <- b_rep_out * (alpha_rep_out / a_rep_out)
      sr_rep_out <- data.frame(rbind(c(alpha_rep_out, beta_rep_out, phi_rep_out, sig_rep_out)))
      colnames(sr_rep_out) <- c("alpha","beta","rho","sigma")

      ## average years & size-at-age
      nref <- 5
      nyrs_obs <- dim(dataObs)[1]
      yrs_averaged <- seq((nyrs_obs - nref + 1), nyrs_obs, 1)
      av_mean_SaAs <- colMeans(meanSaA[yrs_averaged, ])

      ret_by_age_M <- ret_by_age_sex[,, which(names(ret_by_age_sex[1,1,]) == "M")]
      ret_by_age_M[is.na(ret_by_age_M)] <- 0
      ret_by_age_F <- ret_by_age_sex[,, which(names(ret_by_age_sex[1,1,]) == "F")]
      ret_by_age_F[is.na(ret_by_age_F)] <- 0
      ret_propFy_age <- round(ret_by_age_F / (ret_by_age_F + ret_by_age_M), 4)

      av_propF_by_age <- round(colMeans(ret_propFy_age[yrs_averaged, ], na.rm = TRUE), 4)
      av_propF_by_age[is.na(av_propF_by_age)] <- 0

      ret_by_age <- ret_by_age_F + ret_by_age_M
      ret_by_age_prop <- prop.table(ret_by_age, margin = 1)
      ret_by_age_prop[is.na(ret_by_age_prop)] <- 0
      av_p_ret_by_age <- round(colMeans(ret_by_age_prop[yrs_averaged, ], na.rm = TRUE), 4)
      av_p_ret_by_age[is.na(av_p_ret_by_age)] <- 0

      if (sim_recruits == "fecundity") av_rep_out_by_age <- round(exp(.calc_reprod_output(av_mean_SaAs, allometry)[[1]]) * av_propF_by_age)
      if (sim_recruits == "eggmass")   av_rep_out_by_age <- round(exp(.calc_reprod_output(av_mean_SaAs, allometry)[[2]]) * av_propF_by_age)

      ## selectivity by age (may be stochastic; seed for alignment)
      age_selectivity <- withr::with_seed(seednum + 7400L,.calc_selectivity(av_mean_SaAs, maxsel, sdsel))
      age_selectivity <- age_selectivity / max(age_selectivity)

      ## ypr functions
      ypr_func <- function(log_F_max, data) {
        z_as <- av_rep_out_by_age
        pi_as <- av_p_ret_by_age
        zPR0_as <- z_as * pi_as
        zPR0 <- sum(zPR0_as)
        R0 <- log(alpha_rep_out * zPR0) / (beta_rep_out * zPR0)
        F_max <- exp(log_F_max)
        F_as <- F_max * age_selectivity
        U_as <- 1 - exp(-F_as)
        zPR_F_as <- (1 - U_as) * z_as * pi_as
        zPR_F <- sum(zPR_F_as)
        RF <- log(alpha_rep_out * zPR_F) / (beta_rep_out * zPR_F)
        N_as <- RF * pi_as
        S_as <- N_as * (1 - U_as)
        S <- sum(S_as)
        Z_as <- S_as * z_as
        Z <- sum(Z_as)
        H_as <- N_as * U_as
        H <- sum(H_as)
        output <- c(H = H, S = S, R = RF, Z_million = Z / 1e6)
        return(output)
      }
      ypr_max <- function(log_F_max, q, data) { ypr_func(log_F_max, data)["H"] * -1 }

      optim_data <- list(age_selectivity, av_rep_out_by_age, av_p_ret_by_age, alpha_rep_out, beta_rep_out)
      fit <- optim(par = log(1), fn = ypr_max, method = "Brent", lower = -10, upper = 2, data = optim_data)
      F_eq <- fit$par
      eq_out <- ypr_func(F_eq, data = optim_data)
      H_eq <- H_eq <- round(eq_out[1])
      S_eq <- S_eqs <- round(eq_out[2])
      U_eq <- U_eqs <- round(H_eq / (H_eq + S_eq), 4)
    } ## end YPR branch

    ######################################################### set new goal
    if (harvmgmt == "smsy_goal") {
      if (S_msy_obs > msygoal * 0.5 & S_msy_obs < msygoal * 2) { msygoal <- round(S_msy_obs) } else { if (i == 1) { msygoal <- msygoalini } }
    }
    if (harvmgmt == "umsy_goal") {
      if (U_msy_obs > 0 & U_msy_obs < 0.85) { msygoal <- round(U_msy_obs, 2) } else { if (i == 1) { msygoal <- msygoalini } }
    }
    if (harvmgmt == "s_eq_goal") {
      if (S_eq > msygoal * 0.5 & S_eq < msygoal * 2) { msygoal <- round(S_eq) } else { if (i == 1) { msygoal <- msygoalini } }
    }
    if (harvmgmt == "u_eq_goal") {
      if (U_eq > 0 & U_eq < 0.85) { msygoal <- round(U_eq, 2) } else { if (i == 1) { msygoal <- msygoalini } }
    }
    if (harvmgmt == "smsy_dlm_goal") {
      if (Smsy_dlm > msygoal * 0.5 & Smsy_dlm < msygoal * 2) { msygoal <- round(Smsy_dlm) } else { if (i == 1) { msygoal <- msygoalini } }
    }
    if (harvmgmt == "umsy_dlm_goal") {
      if (Umsy_dlm > 0 & Umsy_dlm < 0.85) { msygoal <- round(Umsy_dlm, 2) } else { if (i == 1) { msygoal <- msygoalini } }
    }

    ## save S_msy estimate for this review (to match legacy storage)
    if (harvmgmt == "smsy_goal")      S_msy_estimate[i] <- round(S_msy_obs)
    if (harvmgmt == "s_eq_goal")      S_msy_estimate[i] <- round(S_eq)
    if (harvmgmt == "smsy_dlm_goal")  S_msy_estimate[i] <- round(Smsy_dlm)

    ##=============================================## mgmt strategy factor
    if (harvmgmt %in% c("smsy_goal","s_eq_goal","smsy_dlm_goal")) {
      msygoal <- msygoal * factorMSY
    }

    ##==================================================================##
  } ## end of loop over escapement goal reviews

  ##==================================================================##
  year_index <- (nyi + 1):(nyi + ny - nage - 1) ## reconstructed data period

  ##=============## overall mean size of escapement, harvest, and return
  esc_by_age <- rowSums(esc_by_age_sex, dims = 2); esc_by_age[is.na(esc_by_age)] <- 0; esc_by_age <- esc_by_age[year_index, ]
  ret_by_age <- rowSums(ret_by_age_sex, dims = 2); ret_by_age[is.na(ret_by_age)] <- 0; ret_by_age <- ret_by_age[year_index, ]
  harv_by_age <- rowSums(harv_by_age_sex, dims = 2); harv_by_age[is.na(harv_by_age)] <- 0; harv_by_age <- harv_by_age[year_index, ]
  sizeData <- meanSaA[year_index, ]
  meansize_esc <- meansize_ret <- meansize_harv <- NA
  for (y in 1:nyrec) {
    sizes_esc_y <- sizes_ret_y <- sizes_harv_y <- NA
    for (a in 1:nage) {
      if (esc_by_age[y, a] > 0) { sizes_esc_new <- rep(sizeData[y, a], esc_by_age[y, a]) } else { sizes_esc_new <- NA }
      if (ret_by_age[y, a] > 0) { sizes_ret_new <- rep(sizeData[y, a], ret_by_age[y, a]) } else { sizes_ret_new <- NA }
      if (harv_by_age[y, a] > 0){ sizes_harv_new <- rep(sizeData[y, a], harv_by_age[y, a]) } else { sizes_harv_new <- NA }
      sizes_esc_y <- c(sizes_esc_y, sizes_esc_new)
      sizes_ret_y <- c(sizes_ret_y, sizes_ret_new)
      sizes_harv_y <- c(sizes_harv_y, sizes_harv_new)
    }
    sizes_esc_y <- sizes_esc_y[!is.na(sizes_esc_y)]
    meansize_esc[y] <- mean(as.numeric(sizes_esc_y), na.rm = TRUE)
    sizes_ret_y <- sizes_ret_y[!is.na(sizes_ret_y)]
    meansize_ret[y] <- mean(as.numeric(sizes_ret_y), na.rm = TRUE)
    sizes_harv_y <- sizes_harv_y[!is.na(sizes_harv_y)]
    meansize_harv[y] <- mean(as.numeric(sizes_harv_y), na.rm = TRUE)
  }
  data$meanSesc <- meansize_esc
  data$meanSret <- meansize_ret
  data$meanSharv <- meansize_harv

  ##=======================## change in reproductive potential over time
  use_ny <- nyh
  ny_obs <- dim(dataObs)[1]
  esc_by_age_M <- esc_by_age_sex[,, which(names(esc_by_age_sex[1,1,]) == "M")]
  esc_by_age_F <- esc_by_age_sex[,, which(names(esc_by_age_sex[1,1,]) == "F")]
  esc_by_age_F[is.na(esc_by_age_F)] <- 0
  esc_by_age_F <- esc_by_age_F[year_index, ]
  mean_SaAs <- meanSaA[year_index,]
  fec_by_age <- eggm_by_age <- data.frame(array(NA, dim = c(ny_obs, nage)))
  potential <- data.frame(array(NA, dim = c(ny_obs, 6)))
  names(potential) <- c("Fecundity","Eggmass","AvPerFemFecund","AvPerFemEggmass","AvPerSpawnerFecund","AvPerSpawnerEggmass")
  for (y in 1:ny_obs) {
    sizes_y <- NA
    for (a in 1:nage) {
      if (esc_by_age_F[y, a] > 0) { size_new <- mean_SaAs[y, a] } else { size_new <- NA }
      rep_out_per_cap <- .calc_reprod_output(size_new, allometry)
      fec_by_age[y, a]  <- sum(rep(exp(rep_out_per_cap[[1]]), esc_by_age_F[y, a]))
      eggm_by_age[y, a] <- sum(rep(exp(rep_out_per_cap[[2]]), esc_by_age_F[y, a]))
      sizes_y <- c(sizes_y, rep(size_new, esc_by_age_F[y, a]))
    }
    sizes_y <- sizes_y[!is.na(sizes_y)]
    potential[y, 1] <- sum(exp(.calc_reprod_output(sizes_y, allometry)[[1]]))
    potential[y, 2] <- sum(exp(.calc_reprod_output(sizes_y, allometry)[[2]]))
    potential[y, 3] <- potential[y, 1] / sum(esc_by_age_F[y, ]) ## females
    potential[y, 4] <- potential[y, 2] / sum(esc_by_age_F[y, ]) ## females
    potential[y, 5] <- potential[y, 1] / sum(esc_by_age[y, ])   ## all spawners
    potential[y, 6] <- potential[y, 2] / sum(esc_by_age[y, ])   ## all spawners
  }
  fec_by_age[is.na(fec_by_age)] <- 0
  eggm_by_age[is.na(eggm_by_age)] <- 0
  potential[is.na(potential)] <- 0
  dataObs <- data.frame(cbind(dataObs, round(potential, 2)))

  ##----------------------------------------## percent change calculations
  change <- data.frame(array(NA, dim = c(4, 3)))
  colnames(change) <- c("3yrs","5yrs","10yrs")
  rownames(change) <- c("AvPerFemFecund","AvPerFemEggmass","AvPerSpawnerFecund","AvPerSpawnerEggmass")
  nyrefs <- c(3, 5, 10)
  for (i in 1:4) {
    if (i == 1) trend <- dataObs$AvPerFemFecund[1:use_ny]
    if (i == 2) trend <- dataObs$AvPerFemEggmass[1:use_ny]
    if (i == 3) trend <- dataObs$AvPerSpawnerFecund[1:use_ny]
    if (i == 4) trend <- dataObs$AvPerSpawnerEggmass[1:use_ny]
    for (h in 1:3) {
      nyref <- nyrefs[h]
      earlymean <- mean(trend[1:nyref])
      latemean  <- mean(trend[(use_ny - nyref + 1):use_ny])
      change[i, h] <- round((latemean - earlymean) * 100 / earlymean, 2)
    }
  }

  ##======================================================## save output
  output.list <- list(
    "para"                 = srparms,
    "sr_sim"               = sr_sim,
    "fec"                  = change[3, ],
    "egg"                  = change[4, ],
    "S_msy"                = S_msy_estimate,
    "data"                 = data,
    "obs"                  = dataObs,
    "ret_by_age"           = ret_by_age,
    "meanSaA"              = meanSaA,
    "propfemale"           = propfemale,
    "selectivities_by_age" = selectivities_by_age,
    "MSY_Goals"            = MSY_Goals,
    "impl_errors"          = impl_errors
  )
  return(output.list)
} ## end function

