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

run_model<-function(config) {

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

      set.seed(seednum)

      #################################### generate stock-recruit parameters
      ##====================================## Ricker parameters and regimes
      ##------------------------------------------------------------## alpha
      alpha.low<-rlnorm(1,meanlog=log(alpha_mean),sdlog=sr_parms_sd)
      alpha.high<-alpha.low*regstr ## high productivity regime
      aset<-1 ## initial regime is assumed to be low productivity
      if(reglength==0) aset<-1 ## also without regime shifts
      alpha<-c(alpha.low,alpha.high)[aset] ## select alpha parameter
      ##--------------------------------------------------------## beta/rmax
      if(sr_corr!=0) { ## in case of correlated parameters
        beta_corr<-beta_mean*(1+(alpha-alpha_mean)*sr_corr)
        beta<-rlnorm(1,meanlog=log(beta_corr),sdlog=sr_parms_sd)
      } else {
        beta<-rlnorm(1,meanlog=log(beta_mean),sdlog=sr_parms_sd)
      }
      maxr<-(alpha/beta)*exp(-1) ## Hilborn & Walters 1992 (max recruitment)
      srparms<-c(alpha.low=alpha.low,alpha.high=alpha.high,beta=beta,maxr=maxr)
      ##------------------------------------------------------## S_msy/U_msy
      ## true value in absence of observation error (Scheuerell 2016)
      S_msy_true<-round((1-gsl::lambert_W0(exp(1-log(alpha))))/beta)
      U_msy_true<-round((1-gsl::lambert_W0(exp(1-log(alpha)))),3)

      ##===========================================## alternative parameters
      ## parameters for fecundity/eggmass models (Staton et al. 2021)
      for(i in 1:6){ assign(names(alt_sr_param)[i],as.numeric(alt_sr_param[i])) }
      ##------------------------------------------------------------## Rmax
      rmax_N0<-(alpha_N0/beta_N0)*exp(-1)
      rmax_EASL<-(alpha_EASL/beta_EASL)*exp(-1)
      rmax_EMASL<-(alpha_EMASL/beta_EMASL)*exp(-1)
      ##-------------------------------------------------------## fecundity
      if(sim_recruits=="fecundity"){
      alpha_fec<-(alpha_EASL/alpha_N0)*alpha
      beta_fec<-(beta_EASL/beta_N0)*beta
      rmax_fec<-(rmax_EASL/rmax_N0)*maxr
      }
      ##--------------------------------------------------------## egg mass
      if(sim_recruits=="eggmass"){
      alpha_egg<-(alpha_EMASL/alpha_N0)*alpha
      beta_egg<-(beta_EMASL/beta_N0)*beta
      rmax_egg<-(rmax_EMASL/rmax_N0)*maxr
      }

      #################################### arrays for storing simulated data
      nage<-length(ages)
      sexes<-c("M","F")
      nsex<-length(sexes)
      nyr<-nyi+ny+nage+1 ## additional years dropped during reconstruction
      Yvec<-c(seq(1,nyr))
      ##---------------------------------------------## true population data
      pnames<-c("Year","Ret","Rec","eps","reg","Esc","Harv","RepOut")
      PopDat<-data.frame(array(NA,dim=c(nyr,length(pnames))))
      dimnames(PopDat)<-list(Yvec,c(pnames))
      PopDat$Year<-seq(nyr)
      ##-------------------------------## calendar year age-sex compositions
      ret_by_age_sex<-array(NA,dim=c(nyr,nage,nsex))
      dimnames(ret_by_age_sex)<-list(Yvec,paste0("Ret",ages),sexes)
      harv_by_age_sex<-array(NA,dim=c(nyr,nage,nsex))
      dimnames(harv_by_age_sex)<-list(Yvec,paste0("Harv",ages),sexes)
      esc_by_age_sex<-array(NA,dim=c(nyr,nage,nsex))
      dimnames(esc_by_age_sex)<-list(Yvec,paste0("Esc",ages),sexes)
      ##-----------------------------------------------------## observations
      onames<-c("obsRet","obsEsc","obsHarv","RepOut")
      ObsDat<-data.frame(array(NA,dim=c(nyr,length(onames))))
      dimnames(ObsDat)<-list(c(seq(1,nyr)),onames)
      ##----------------------------------------------## harvest selectivity
      selectivities_by_age<-array(0,dim=c(nyr,nage)) ## selectivity by age

      ################################################### demographic trends
      max_yh<-nyi+nyh ## last year of historical period
      yrs_h<-(nyi+1):max_yh ## historical years index
      yrs_f<-(nyi+nyh+1):nyr ## future years index

      ##=========================================================## mean age
      ## deterministic trend but age proportions drawn based on mean by year
      meanage<-NA
      meanage[1:nyi]<-meanageini
      for(y in yrs_h) { meanage[y]<-meanage[y-1]+agetrend/(nyh) }
      if(futureT=="yes") { for(y in yrs_f) meanage[y]<-meanage[y-1]+agetrend/(nyh) } else { meanage[yrs_f]<-meanage[max_yh] }

      ##===============================================## proportion female
      ## propF trend given as change over entire period in logit space
      propF_y<-NA;propF_y[1:nyi]<-stats::qlogis(propF)
      n_ref<-5
      # for(y in yrs_h){ propF_y[y]<-propF_y[y-1]+propFtrend/(ny) }
      for(y in yrs_h){ propF_y[y]<-propF_y[y-1]+propFtrend/(nyh) }
      if(futureT=="yes"){ for(y in yrs_f) propF_y[y]<-propF_y[y-1]+propFtrend/(nyh) } else { propF_y[yrs_f]<-mean(propF_y[(max_yh-n_ref+1):max_yh]) }
      propF_y<-propF_y+rnorm(nyr,0,0.05) ## with normal error on logit scale
      propfemale<-stats::plogis(propF_y)

      ##=================================================## mean size-at-age
      meanSaA<-meanSaAdet<-SaA_anoms<-array(0,dim=c(nyr,nage))
      dimnames(meanSaA)<-list(Yvec,paste0("age",ages))
      ##----------------------------------------## initial mean sizes at age
      oceanage<-ages[1:(length(ages)-2)] ## ocean ages
      growth<-function(age,Linf,k){ Linf*(1-exp(-k*age)) } ## growth function
      meansini<-round((ocean0s+growth(oceanage,Linf=(vonB_Linf-ocean0s),k=vonB_k))/10)*10 ## rounded to cm
      meanSaA_ini<-c(0,ocean0s,meansini) ## fish enter the ocean at age 2 (BY+2)
      ##--------------------------------------## size-at-age trend (mm/year)
      trendsSaA<-sizetrends/(nyh) ## for reconstructed time period
      ##--------------------------------------## deterministic size at age
      for(y in 1:nyi){ meanSaAdet[y,]<-meanSaA_ini } ## same for initial years
      for(y in yrs_h){ for(a in 1:nage){ meanSaAdet[y,a]<-meanSaAdet[y-1,a]+ trendsSaA[a] } }
      if(futureT=="yes"){ for(y in yrs_f){ for(a in 1:nage){ meanSaAdet[y,a]<-meanSaAdet[y-1,a]+ trendsSaA[a] } } } else { for(a in 1:nage){ meanSaAdet[yrs_f,a]<-mean(meanSaAdet[(max_yh-n_ref+1):max_yh,a]) } }
      ##-------------------------------------------## size-at-age anomalies
      for(a in 1:nage){ SaA_anoms[,a]<-rnorm(nyr,mean=0,sd=meanSaA_ini[a]*sdSaA) }
      ##-------------------------------------------## stochastic size at age
      meanSaA<-meanSaAdet+SaA_anoms

      ##==================================================================##
      ##=========================================## generate population data
      ##==================================================================##
      age_comp<-array(0,dim=c(nyr,nage,nsex))

      ##################################### generate first few years of data
      PopDat$Esc[1:nyi]<-round(exp(rnorm(nyi,log(S_msy_true)-0.5*0.25^2,0.25)))
      for(y in 1:nyi) {
      PopDat$reg[y]<-aset ## regime (no regime shifts yet)
      ##----------------------------------## generate recruits by brood year
      if(sim_recruits=="spawners"){
      sr<-SizeShiftsImplicationsV2::ricker(spawn=PopDat$Esc[y],sigma=procerr,alpha=alpha,beta=beta,rho=rho,last.eps=ifelse(y==1,0,PopDat$eps[y-1]))
      PopDat$RepOut[y]<-PopDat$Esc[y]
      }
      ##----------------------------------## alternative reproductive units
      if(sim_recruits!="spawners"){
      age_comp_y<-SizeShiftsImplicationsV2::agecomp(ages,round(PopDat$Esc[y]*(1-propfemale[y])),meanage[y]-agediff/2,sdage) ## draw age structure to simulate initial sizes
      sizes_y<-rep(meanSaA[y,],prop.table(age_comp_y)*PopDat$Esc[y])
      ##-----------------------------## generate recruits based on fecundity
      if(sim_recruits=="fecundity"){
      fecundity_total<-round(sum(exp(SizeShiftsImplicationsV2::reprod_output(sizes_y,allometry)[[1]])))
      if(ricker_type=="const_beta")sr<-SizeShiftsImplicationsV2::ricker(spawn=fecundity_total,sigma=procerr,alpha=alpha_fec,beta=beta_fec,rho=rho,last.eps=ifelse(y==1,0,PopDat$eps[y-1]))
      PopDat$RepOut[y]<-fecundity_total
      }
      ##------------------------------## generate recruits based on egg mass
      if(sim_recruits=="eggmass"){
      eggmass_total<-round(sum(exp(SizeShiftsImplicationsV2::reprod_output(sizes_y,allometry)[[2]])))
      if(ricker_type=="const_beta")sr<-SizeShiftsImplicationsV2::ricker(spawn=eggmass_total,sigma=procerr,alpha=alpha_egg,beta=beta_egg,rho=rho,last.eps=ifelse(y==1,0,PopDat$eps[y-1]))
      PopDat$RepOut[y]<-eggmass_total
      }
      } ## end if 'sim_recruits' statement
      ##---------------------------------## store recruits from brood year y
      PopDat$Rec[y]<-round(sr$rec)
      ##-------------------## store anomaly (epsilon for AR1 process year y)
      PopDat$eps[y]<-round(sr$eps,5)
      ##-------------------## generate sex-specific recruit age compositions
      age_comp[y,,1]<-SizeShiftsImplicationsV2::agecomp(ages,round(PopDat$Rec[y]*(1-propfemale[y])),meanage[y]-agediff/2,sdage) ## males
      age_comp[y,,2]<-SizeShiftsImplicationsV2::agecomp(ages,round(PopDat$Rec[y]*propfemale[y]),meanage[y]+agediff/2,sdage) ## females
      for(k in 1:nage) {
        ret_by_age_sex[y+k,k,1]<-age_comp[y,k,1]
        ret_by_age_sex[y+k,k,2]<-age_comp[y,k,2]
      }
      ##-----------------------## generate escapement observations from data
      ObsDat$obsEsc[y]<-PopDat$Esc[y]*exp(rnorm(1,-0.5*obserr^2,obserr))
      } ## end loop over initial years

      ##################################### years for re-assessing MSY goals
      if(harvmgmt=="fix_harv_rate") goalfreq<-F
      firstrev<-20
      if(goalfreq) { goalrev<-seq(nyi+firstrev,ny,goalfreq) } else { goalrev<-ny }
      nrev<-length(goalrev)
      ##-------## set initial escapement goal as mean of initial escapements
      if(harvmgmt %in% c("smsy_goal","s_eq_goal","smsy_dlm_goal")) {
      msygoal<-S_msy_true;msygoalini<-msygoal
      }
      if(harvmgmt %in% c("umsy_goal","u_eq_goal","umsy_dlm_goal")) {
      msygoal<-U_msy_true;msygoalini<-msygoal
      }
      ##-------------------------------------## storing errors, bias, goals
      impl_errors<-MSY_Goals<-HarvRates<-S_msy_estimate<-NA

      ######################################## loop through esc goal reviews
      for(i in 1:nrev) {
      ##=========================================## calculate years of data
      y_rev<-goalrev[i]
      if(i==1) nycnt<-nyi else nycnt<-nyi+goalrev[i-1]
      yindex<-(nycnt+1):(nyi+y_rev)
      ##===============================## generate subsequent years of data
      MSY_Goals[yindex]<-msygoal
      for(y in yindex) {
      ## if regime shift switch productivity and change regime indicator
      if(reglength!=0){
      if(sample(c(rep(0,reglength-1),1),1)){
      if(aset==1) alpha<-alpha.low else alpha<-alpha.high
      aset<-abs(aset-3) ## change regime indicator
      }
      }
      PopDat$reg[y]<-aset ## regime
      ##===============================## total return across ages and sexes
      PopDat$Ret[y]<-sum(ret_by_age_sex[y,,],na.rm=T)

      ##==============================## management using fixed harvest rate
      if(harvmgmt=="fix_harv_rate"){
      ## fishery implementation error on logit scale
      impl_err<-impl_errors[y]<-rnorm(1,mean=0,sd=harverr)
      harv_thresh<-round(maxr*0.25) ## 25% max recruits as conservation threshold
      if(PopDat$Ret[y]>=harv_thresh) { harvrate_det<-harvrate } else {
      harvrate_det<-harvrate*(PopDat$Ret[y]/(PopDat$Ret[y]+round(maxr*0.1))) }
      harvrate_real<-stats::plogis(stats::qlogis(harvrate_det)+impl_err)
      if(is.na(harvrate_real)) harvrate_real<-0
      HarvRates[y]<-harvrate_real
      PopDat$Harv[y]<-round(PopDat$Ret[y]*harvrate_real) ## total harvest
      PopDat$Esc[y]<-PopDat$Ret[y]-PopDat$Harv[y] ## escapement=return-harvest
      }  ## end if statement

      ##==================================## management by harvest rate goal
      if(harvmgmt %in% c("umsy_goal","u_eq_goal","umsy_dlm_goal")){
      ## fishery implementation error on logit scale
      impl_err<-impl_errors[y]<-rnorm(1,mean=0,sd=harverr)
      harvrate_real<-stats::plogis(stats::qlogis(msygoal)+impl_err) ## realized harvest rate
      PopDat$Harv[y]<-round(PopDat$Ret[y]*harvrate_real) ## total harvest
      PopDat$Esc[y]<-PopDat$Ret[y]-PopDat$Harv[y] ## escapement=return-harvest
      }  ## end if statement

      ##===================================## management by escapement goal
      if(harvmgmt %in% c("smsy_goal","s_eq_goal","smsy_dlm_goal")){
      ## fishery implementation normal error and bias on log-scale
      impl_err<-impl_errors[y]<-rnorm(1,mean=-0.5*harverr^2,sd=harverr)
      bias<-0
      escape<-round(msygoal*exp(impl_err+bias))
      PopDat$Esc[y]<-min(escape,PopDat$Ret[y]) ## no larger than return
      PopDat$Harv[y]<-PopDat$Ret[y]-PopDat$Esc[y]
      } ## end if statement

      ##===========================================## harvest by age and sex
      ## harvest by age based on size-at-age (unequal probability sampling)
      ##-----------------------------------## age-sex proportions in return
      ret_prop_by_sex_y<-prop.table(ret_by_age_sex[y,,],1) ## this year
      ##----------------------------------## sum return by age across sexes
      ret_by_age_all_y<-rowSums(ret_by_age_sex[y,,]) ## this year
      ##------------------------------------------## approximate 'sampling'
      fac<-10 ## approximation to speed up the computation
      nsample<-round(PopDat$Harv[y]/fac)
      return_by_age<-round(ret_by_age_all_y/fac)
      return_by_age[is.na(return_by_age)]<-0
      age_vector<-rep(ages,return_by_age)
      ##------------------------------------------------## apply selectivity
      selectivities_by_age[y,]<-SizeShiftsImplicationsV2::selectivity(meanSaA[y,],maxsel,sdsel)
      select<-rep(selectivities_by_age[y,],return_by_age)
      ##-------------------------------## sample harvest without replacement
      harvested<-sample(age_vector,size=nsample,prob=select,replace=F)
      tab<-table(harvested)
      ##------------------------------------## redistribute harvest to sexes
      ## using sex-by-age proportions prior to harvest
      harv_prop_M<-ret_prop_by_sex_y[as.numeric(names(tab)),1]
      harv_by_age_sex[y,as.numeric(names(tab)),1]<-round(tab*fac*harv_prop_M)
      harv_prop_F<-ret_prop_by_sex_y[as.numeric(names(tab)),2]
      harv_by_age_sex[y,as.numeric(names(tab)),2]<-round(tab*fac*harv_prop_F)
      harv_by_age_sex[is.na(harv_by_age_sex)]<-0
      ##-----------------------## avoid over-harvesting of age-by-sex groups
      ## only needed with approximation
      index<-which(harv_by_age_sex[y,,1]>ret_by_age_sex[y,,1])
      if(length(index)>=1){harv_by_age_sex[y,index,1]<-ret_by_age_sex[y,index,1]}
      index<-which(harv_by_age_sex[y,,2]>ret_by_age_sex[y,,2])
      if(length(index)>=1){harv_by_age_sex[y,index,2]<-ret_by_age_sex[y,index,2]}

      ##=======================================## escapement by age and size
      esc_by_age_sex[y,,]<-ret_by_age_sex[y,,]-harv_by_age_sex[y,,]

      ##================================================## generate recruits
      if(sim_recruits=="spawners"){
      sr<-SizeShiftsImplicationsV2::ricker(spawn=PopDat$Esc[y],sigma=procerr,alpha=alpha,beta=beta,rho=rho,last.eps=ifelse(y==1,0,PopDat$eps[y-1]))
      PopDat$RepOut[y]<-PopDat$Esc[y]
      }
      ##===================================## alternative reproductive units
      if(sim_recruits!="spawners"){
      ##---------------------------------------------## all individual sizes
      sizes_y<-NA
      for(a in 1:nage) {
      esc_F_a_y<-esc_by_age_sex[y,a,which(names(esc_by_age_sex[1,1,])=="F")]
      if(is.na(esc_F_a_y)) { size_new<-NA  } else { if(esc_F_a_y>0) { size_new<-meanSaA[y,a] } else { size_new<-NA } }
      sizes_y<-c(sizes_y,rep(size_new,esc_F_a_y))
      }
      sizes_y<-sizes_y[!is.na(sizes_y)]
      ##-----------------------------## generate recruits based on fecundity
      if(sim_recruits=="fecundity"){
      fecundity_total<-round(sum(exp(SizeShiftsImplicationsV2::reprod_output(sizes_y,allometry)[[1]])))
      if(ricker_type=="const_beta")sr<-SizeShiftsImplicationsV2::ricker(spawn=fecundity_total,sigma=procerr,alpha=alpha_fec,beta=beta_fec,rho=rho,last.eps=ifelse(y==1,0,PopDat$eps[y-1]))
      PopDat$RepOut[y]<-fecundity_total
      }
      ##------------------------------## generate recruits based on egg mass
      if(sim_recruits=="eggmass"){
      eggmass_total<-round(sum(exp(SizeShiftsImplicationsV2::reprod_output(sizes_y,allometry)[[2]])))
      if(ricker_type=="const_beta")sr<-SizeShiftsImplicationsV2::ricker(spawn=eggmass_total,sigma=procerr,alpha=alpha_egg,beta=beta_egg, rho=rho,last.eps=ifelse(y==1,0,PopDat$eps[y-1]))
      PopDat$RepOut[y]<-eggmass_total
      }

      } ## end if 'sim_recruits' statement

      ##---------------------------------## store recruits from brood year y
      PopDat$Rec[y]<-round(sr$rec)
      ##-------------------## store anomaly (epsilon for AR1 process year y)
      PopDat$eps[y]<-round(sr$eps,5)
      ##--------------------## generate sex-specific recruit age composition
      age_comp[y,,1]<-SizeShiftsImplicationsV2::agecomp(ages,round(PopDat$Rec[y]*(1-propfemale[y])),meanage[y]-agediff/2,sdage) ## males
      age_comp[y,,2]<-SizeShiftsImplicationsV2::agecomp(ages,round(PopDat$Rec[y]*propfemale[y]),meanage[y]+agediff/2,sdage) ## females
      for(k in 1:nage) {
        ret_by_age_sex[y+k,k,1]<-age_comp[y,k,1]
        ret_by_age_sex[y+k,k,2]<-age_comp[y,k,2]
      }
      ##=================================## generate observations from data
      ## multiply true escapement by escapement observation error
      ObsDat$obsEsc[y]<-PopDat$Esc[y]*exp(rnorm(1,-0.5*obserr^2,obserr))
      ## multiply true harvest by harvest observation error
      ObsDat$obsHarv[y]<-PopDat$Harv[y]*exp(rnorm(1,-0.5*hobserr^2,hobserr))
      ## store observed return (harvest + escapement)
      ObsDat$obsRet[y]<-ObsDat$obsEsc[y]+ObsDat$obsHarv[y]
      ##----------------------------## what if fecundity/eggmass are known?
      ObsDat$RepOut[y]<-PopDat$RepOut[y]
      } ## end year loop

      ##==========================## round population numbers to individuals
      roundit<-c("Ret","Rec","Esc","Harv","RepOut")
      PopDat[,names(PopDat)%in%roundit]<-round(PopDat[,names(PopDat)%in%roundit])
      ObsDat<-round(ObsDat)

      ##========================================## simulated population data
      data_all<-PopDat
      data_all[is.na(data_all)]<-0
      year_index<-(nyi+1):(nyi+y_rev-nage-1) ## match reconstructed data
      nyrec<-length(year_index) ## number of reconstructed years
      data<-data.frame(data_all[year_index,])

      ##################################################### assess MSY goals
      ##=====================================## Ricker fit to simulated data
      data$lnRS<-log(data$Rec/data$Esc)
      ##--------------------------------------------------------------## GLS
      mod_sim<-nlme::gls(lnRS~Esc,data=data,correlation=nlme::corAR1())
      sig_sim<-summary(mod_sim)$sigma ## sigma to correct for log-normal error
      phi_sim<-coef(mod_sim$modelStruct$corStruct,unconstrained=FALSE)
      log_a_sim<-summary(mod_sim)$coefficients[1] ## log-productivity
      a_sim<-exp(log_a_sim) ## productivity (alpha)
      b_sim<- -summary(mod_sim)$coefficients[2] ## density-dependence (beta)
      alpha_sim<-exp(log_a_sim+((0.5*sig_sim^2)/(1-phi_sim^2))) ## corrected
      # alpha_sim<-exp(log_a_sim+(sig_sim^2)/(2*(1-phi_sim^2))) ## same
      beta_sim<-b_sim*(alpha_sim/a_sim) ## corrected beta value
      S_msy_sim<-(1-gsl::lambert_W0(exp(1-log(alpha_sim))))/beta_sim
      U_msy_sim<-(1-gsl::lambert_W0(exp(1-log(alpha_sim))))
      S_max_sim<-1/beta_sim
      ##----------------------------------------------------------## results
      sr_sim<-data.frame(rbind(c(alpha_sim,beta_sim,phi_sim,sig_sim)))
      colnames(sr_sim)<-c("alpha","beta","rho","sigma")
      ##===================================================## reconstruction
      ## reconstruct recruitment from year y using observed escapement and harvest
      dataObs<-ObsDat ## all years simulated
      dataObs[is.na(dataObs)]<-0
      ##--------------------------------------## escapement age proportions
      age_comp_esc<-rowSums(esc_by_age_sex,dims=2)
      age_comp_esc[is.na(age_comp_esc)]<-0
      age_prop_esc<-prop.table(as.matrix(age_comp_esc),margin=1)
      age_comp_esc_obs<-array(0,dim=dim(age_prop_esc))
      for(y in 1:nyr){age_comp_esc_obs[y,]<-age_prop_esc[y,]*dataObs$obsEsc[y]}
      age_comp_esc_obs[is.na(age_comp_esc_obs)]<-0
      age_comp_esc_obs<-round(age_comp_esc_obs)
      ##------------------------------------------## harvest age proportions
      age_comp_harv<-rowSums(harv_by_age_sex,dims=2)
      age_comp_harv[is.na(age_comp_harv)]<-0
      age_prop_harv<-prop.table(as.matrix(age_comp_harv),margin=1)
      age_comp_harv_obs<-array(0,dim=dim(age_prop_harv))
      for(y in 1:nyr){age_comp_harv_obs[y,]<-age_prop_harv[y,]*dataObs$obsHarv[y]}
      age_comp_harv_obs[is.na(age_comp_harv_obs)]<-0
      age_comp_harv_obs<-round(age_comp_harv_obs)
      ##----------------------------------------## recruitment by brood year
      ## drop last nage years for full reconstruction
      nage<-length(ages)
      dataObs$recRec<-NA
      age_comp_ret_obs_byBY<-array(NA,dim=dim(age_comp_harv_obs))
      for(yy in (nyi+1):(nyi+ny-nage)) {
      sum_year_ages<-0
      for(a in 1:nage){
      sum_broodyear<-age_comp_esc_obs[yy+a,a]+age_comp_harv_obs[yy+a,a]
      age_comp_ret_obs_byBY[yy,a]<-sum_broodyear
      sum_year_ages<-sum_year_ages+sum_broodyear
      } ## end age loop
      dataObs$recRec[yy]<-sum_year_ages
      } ## end year loop
      dataObs<-data.frame(dataObs[year_index,])
      dataObs<-dataObs[!is.na(dataObs$recRec),]
      dataObs<-dataObs[dataObs$recRec>0,]

      ##================================## Ricker fits to reconstructed data
      ## fewer years of data due to reconstruction
      dataObs$obsRpS<-dataObs$recRec/dataObs$obsEsc
      dataObs$obslnRS<-log(dataObs$recRec/dataObs$obsEsc)
      ##---------------------------------------------------------## MSY gols
      if(harvmgmt %in% c("smsy_goal","umsy_goal")){
      ##--------------------------------------------------------------## GLS
      mod_obs<-nlme::gls(obslnRS~obsEsc,data=dataObs,correlation=nlme::corAR1())
      sig_obs<-summary(mod_obs)$sigma ## sigma to correct for log-normal error
      phi_obs<-coef(mod_obs$modelStruct$corStruct,unconstrained=FALSE)
      dataObs$mod_obs_resid<-residuals(mod_obs)
      log_a_obs<-summary(mod_obs)$coefficients[1] ## log-productivity
      a_obs<-exp(log_a_obs) ## productivity (alpha)
      b_obs<- -summary(mod_obs)$coefficients[2] ## density-dependence (beta)
      alpha_obs<-exp(log_a_obs+((0.5*sig_obs^2)/(1-phi_obs^2))) ## corrected alpha
      beta_obs<-b_obs*(alpha_obs/a_obs) ## corrected beta value
      S_msy_obs<-(1-gsl::lambert_W0(exp(1-log(alpha_obs))))/beta_obs
      U_msy_obs<-(1-gsl::lambert_W0(exp(1-log(alpha_obs))))
      S_max_obs<-1/beta_obs
      sr_obs<-data.frame(rbind(c(alpha_obs,beta_obs,phi_obs,sig_obs)))
      colnames(sr_obs)<-c("alpha","beta","rho","sigma")
      } ## end if statement 'harvmgmt'

      ##==============================================================## DLM
      ## dynamic linear model with time-varying alpha and/or beta parameters
      if(harvmgmt %in% c("smsy_dlm_goal","umsy_dlm_goal")){
      dataDLM<-dplyr::select(dataObs,Rec=recRec,Esc=obsEsc)
      nrd<-dim(dataDLM)[1]
      if(var=="alpha") mod_dlm<-suppressWarnings(SizeShiftsImplicationsV2::DLMfit(data=dataDLM,var_alpha=TRUE, var_beta=FALSE))
      if(var=="beta") mod_dlm<-suppressWarnings(SizeShiftsImplicationsV2::DLMfit(data=dataDLM,var_alpha=FALSE, var_beta=TRUE))
      if(var=="both") mod_dlm<-suppressWarnings(SizeShiftsImplicationsV2::DLMfit(data=dataDLM,var_alpha=TRUE, var_beta=TRUE))
      sig_dlm<-mod_dlm$sigma
      log_a_dlm<-mod_dlm$results$alpha_y ## log-productivity
      a_dlm<-exp(log_a_dlm) ## productivity (alpha)
      b_dlm<- -mod_dlm$results$beta_y ## density-dependence (beta)
      alpha_dlm<-signif(exp(log_a_dlm+((0.5*sig_dlm^2))),4) ## corrected
      beta_dlm<-signif(b_dlm*(alpha_dlm/a_dlm),4) ## corrected
      S_msy_dlm<-round((1-gsl::lambert_W0(exp(1-log(alpha_dlm))))/beta_dlm)
      U_msy_dlm<-round((1-gsl::lambert_W0(exp(1-log(alpha_dlm)))),3)
      dlm_out<-data.frame(cbind(alpha_dlm,beta_dlm,S_msy_dlm,U_msy_dlm))
      ##------------------------------------------## recent years Smsy/Umsy
      ndlm<-5
      alpha.dlm<-mean(alpha_dlm[(nrd-ndlm+1):nrd])
      beta.dlm<-mean(beta_dlm[(nrd-ndlm+1):nrd])
      Smsy_dlm<-round((1-gsl::lambert_W0(exp(1-log(alpha.dlm))))/beta.dlm)
      Umsy_dlm<-round((1-gsl::lambert_W0(exp(1-log(alpha.dlm)))),3)
      Smax_dlm<-1/beta.dlm
      }

      ########################################## yield-per-recruit analysis
      if(harvmgmt %in% c("s_eq_goal","u_eq_goal")){
      ##----------------## Ricker fit to reproductive output unit using GLS
      dataObs$lnRperRU<-log(dataObs$recRec/dataObs$RepOut)
      mod_rep_out<-nlme::gls(lnRperRU~RepOut,data=dataObs,correlation=nlme::corAR1())
      sig_rep_out<-summary(mod_rep_out)$sigma
      phi_rep_out<-coef(mod_rep_out$modelStruct$corStruct,unconstrained=FALSE)
      log_a_rep_out<-summary(mod_rep_out)$coefficients[1] ## log-productivity
      a_rep_out<-exp(log_a_rep_out) ## productivity
      b_rep_out<--summary(mod_rep_out)$coefficients[2] ## density-dependence
      alpha_rep_out<-exp(log_a_rep_out+((0.5*sig_rep_out^2)/(1-phi_rep_out^2)))
      beta_rep_out<-b_rep_out*(alpha_rep_out/a_rep_out) ## correct beta value
      sr_rep_out<-data.frame(rbind(c(alpha_rep_out,beta_rep_out,phi_rep_out,sig_rep_out)))
      colnames(sr_rep_out)<-c("alpha","beta","rho","sigma")
      ##--------------------------------------## average reproductive output
      nref<-5
      nyrs_obs<-dim(dataObs)[1]
      yrs_averaged<-seq((nyrs_obs-nref+1),nyrs_obs,1) ## last 10 years
      ##----------------------------------------------## average size-at-age
      av_mean_SaAs<-colMeans(meanSaA[yrs_averaged,])
      ##-----------------------------------------## return by sex and by age
      ret_by_age_M<-ret_by_age_sex[,,which(names(ret_by_age_sex[1,1,])=="M")]
      ret_by_age_M[is.na(ret_by_age_M)]<-0
      ret_by_age_F<-ret_by_age_sex[,,which(names(ret_by_age_sex[1,1,])=="F")]
      ret_by_age_F[is.na(ret_by_age_F)]<-0
      ret_propFy_age<-round(ret_by_age_F/(ret_by_age_F+ret_by_age_M),4)
      ##------------## average proportion female by age for reference period
      av_propF_by_age<-round(colMeans(ret_propFy_age[yrs_averaged,],na.rm=T),4)
      av_propF_by_age[is.na(av_propF_by_age)]<-0
      ##--------------------------## average probability of returning by age
      ret_by_age<-ret_by_age_F+ret_by_age_M
      ret_by_age_prop<-prop.table(ret_by_age,margin=1)
      ret_by_age_prop[is.na(ret_by_age_prop)]<-0
      av_p_ret_by_age<-round(colMeans(ret_by_age_prop[yrs_averaged,],na.rm=T),4)
      av_p_ret_by_age[is.na(av_p_ret_by_age)]<-0
      ##-------------------## average per spawner reproductive output by age
      if(sim_recruits=="fecundity") av_rep_out_by_age<-round(exp(SizeShiftsImplicationsV2::reprod_output(av_mean_SaAs,allometry)[[1]])*av_propF_by_age)
      if(sim_recruits=="eggmass") av_rep_out_by_age<-round(exp(SizeShiftsImplicationsV2::reprod_output(av_mean_SaAs,allometry)[[2]])*av_propF_by_age)
      ##-----------------------------------------------## selectivity by age
      age_selectivity<-SizeShiftsImplicationsV2::selectivity(av_mean_SaAs,maxsel,sdsel)
      age_selectivity<-age_selectivity/max(age_selectivity)
      ##-------------------------------------------------## maximize harvest
      ##--------------------## function for 'yield-per-recruit' analyses
      ypr_func<-function(log_F_max,data) {
      z_as<-av_rep_out_by_age ## average reproductive output by age
      pi_as<-av_p_ret_by_age ## average probability of returning by age
      zPR0_as<-z_as*pi_as
      zPR0<-sum(zPR0_as) ## unfished z per recruit
      R0<-log(alpha_rep_out*zPR0)/(beta_rep_out*zPR0) ## unfished eq rec
      F_max<-exp(log_F_max) ## instantaneous F of fully-selected age
      F_as<-F_max*age_selectivity ## instantaneous F of each age
      U_as<-1-exp(-F_as) ## proportion of each age harvested
      zPR_F_as<-(1-U_as)*z_as*pi_as ## reprod. output per recruit by age
      zPR_F<-sum(zPR_F_as) ## reprod. output per recruit after fishing
      RF<-log(alpha_rep_out*zPR_F)/(beta_rep_out*zPR_F) ## fished eq recruitment
      N_as<-RF*pi_as ## equilibrium return by age
      S_as<-N_as*(1-U_as) ## equilibrium spawners by age
      S<-sum(S_as) ## total equilibrium spawners: Smsy if maximizing harvest
      Z_as<-S_as*z_as ## total equilibrium reproductive output by age
      Z<-sum(Z_as) ## total equilibrium reproductive output
      H_as<-N_as*U_as ## equilibrium harvest by age
      H<-sum(H_as) ## total equilibrium harvest
      output=c(H=H,S=S,R=RF,Z_million=Z/1e6) ## bundle the output
      return(output) ## return the output
      }
      ##--------------## function to optimize ("H" for Smsy | "R" for Smax)
      ypr_max<-function(log_F_max,q,data) {
      ypr_func(log_F_max,data)["H"] * -1
      } ## multiply by -1 to maximize a quantity (instead of minimize)
      optim_data<-list(age_selectivity,av_rep_out_by_age,av_p_ret_by_age, 		alpha_rep_out,beta_rep_out)
      fit<-optim(par=log(1),fn=ypr_max,method="Brent",lower=-10,upper=2,data=optim_data)
      # fit<-nlminb(start=log(1),objective=ypr_max,lower=-10,upper=2,data=optim_data)
      F_eq<-fit$par ## F that maximizes equilibrium harvest or recruitment
      # if(fit$convergence==0) print("optim function converged")
      eq_out<-ypr_func(F_eq,data) ## use F_eq to get equilibrium values
      H_eq<-H_eq<-round(eq_out[1])
      S_eq<-S_eqs<-round(eq_out[2])
      U_eq<-U_eqs<-round(H_eq/(H_eq+S_eq),4)
      } ## end if statement for 'harvmgmt'

      ######################################################### set new goal
      ## for management strategies that rely on re-evaluating MSY goal
      ##------------------------------------------------------------## S_msy
      if(harvmgmt=="smsy_goal"){
      if(S_msy_obs>msygoal*0.5 & S_msy_obs<msygoal*2){ msygoal<-round(S_msy_obs) } else { if(i==1) { msygoal<-msygoalini } else { msygoal<-msygoal } }
      } ## set S_msy as new goal if within 0.5-2 times old goal
      ##------------------------------------------------------------## U_msy
      if(harvmgmt=="umsy_goal"){
      if(U_msy_obs>0 & U_msy_obs<0.85){ msygoal<-round(U_msy_obs,2) } else { if(i==1) { msygoal<-msygoalini } else { msygoal<-msygoal } }
      } ## set U_msy as new goal if smaller than 85% harvest rate
      ##-------------------------------------------------------------## S_eq
      if(harvmgmt=="s_eq_goal"){
      if(S_eq>msygoal*0.5 & S_eq<msygoal*2){ msygoal<-round(S_eq) } else { if(i==1) { msygoal<-msygoalini } else { msygoal<-msygoal } }
      } ## set S_eq as new goal if within 0.5-2 times old goal
      ##-------------------------------------------------------------## U_eq
      if(harvmgmt=="u_eq_goal"){
      if(U_eq>0 & U_eq<0.85){ msygoal<-round(U_eq,2) } else { if(i==1) { msygoal<-msygoalini } else { msygoal<-msygoal } }
      } ## set U_eq as new goal if smaller than 85% harvest rate
      ##---------------------------------------------------------## Smsy_dlm
      if(harvmgmt=="smsy_dlm_goal"){
        if(Smsy_dlm>msygoal*0.5 & Smsy_dlm<msygoal*2){ msygoal<-round(Smsy_dlm) } else { if(i==1) { msygoal<-msygoalini } else { msygoal<-msygoal } }
      } ## set S_eq as new goal if within 0.5-2 times old goal
      ##---------------------------------------------------------## Umsy_dlm
      if(harvmgmt=="umsy_dlm_goal"){
        if(Umsy_dlm>0 & Umsy_dlm<0.85){ msygoal<-round(Umsy_dlm,2) } else { if(i==1) { msygoal<-msygoalini } else { msygoal<-msygoal } }
      } ## set U_eq as new goal if smaller than 85% harvest rate
      ##-------------------------------------------## S_msy estimate to save
      if(harvmgmt=="smsy_goal") S_msy_estimate[i]<-round(S_msy_obs)
      if(harvmgmt=="s_eq_goal") S_msy_estimate[i]<-round(S_eq)
      if(harvmgmt=="smsy_dlm_goal") S_msy_estimate[i]<-round(Smsy_dlm)

      ##=============================================## management strategy
      ## conservative (0.8*MSY), neutral, aggressive (1.2*MSY)
      if(harvmgmt %in% c("smsy_goal","s_eq_goal","smsy_dlm_goal")){
      msygoal<-msygoal*factorMSY
      }

      ##==================================================================##
      } ## end of loop over escapement goal reviews > end of simulation part

      ##==================================================================##
      year_index<-(nyi+1):(nyi+ny-nage-1) ## reconstructed data period

      ##==================================================================##
      ##=============## overall mean size of escapement, harvest, and return
      ##==================================================================##
      ## means across sexes, i.e. not sex-specific
      ##----------------------------------------------## age comp escapement
      esc_by_age<-rowSums(esc_by_age_sex,dims=2)
      esc_by_age[is.na(esc_by_age)]<-0
      esc_by_age<-esc_by_age[year_index,]
      ##--------------------------------------------------## age comp return
      ret_by_age<-rowSums(ret_by_age_sex,dims=2)
      ret_by_age[is.na(ret_by_age)]<-0
      ret_by_age<-ret_by_age[year_index,]
      ##--------------------------------------------------## age comp return
      harv_by_age<-rowSums(harv_by_age_sex,dims=2)
      harv_by_age[is.na(harv_by_age)]<-0
      harv_by_age<-harv_by_age[year_index,]
      ##--------------------------------------------------------## size data
      sizeData<-meanSaA[year_index,]
      ##----------------------------------------------## loop years and ages
      meansize_esc<-meansize_ret<-meansize_harv<-NA
      for(y in 1:nyrec) {
      sizes_esc_y<-sizes_ret_y<-sizes_harv_y<-NA
      for(a in 1:nage) {
      if(esc_by_age[y,a]>0) { sizes_esc_new<-rep(sizeData[y,a],esc_by_age[y,a]) } else { sizes_esc_new<-NA }
      if(ret_by_age[y,a]>0) { sizes_ret_new<-rep(sizeData[y,a],ret_by_age[y,a]) } else { sizes_ret_new<-NA }
      if(harv_by_age[y,a]>0) {sizes_harv_new<-rep(sizeData[y,a],harv_by_age[y,a]) } else { sizes_harv_new<-NA }
      sizes_esc_y<-c(sizes_esc_y,sizes_esc_new)
      sizes_ret_y<-c(sizes_ret_y,sizes_ret_new)
      sizes_harv_y<-c(sizes_harv_y,sizes_harv_new)
      }
      sizes_esc_y<-sizes_esc_y[!is.na(sizes_esc_y)]
      meansize_esc[y]<-mean(as.numeric(sizes_esc_y),na.rm=T)
      sizes_ret_y<-sizes_ret_y[!is.na(sizes_ret_y)]
      meansize_ret[y]<-mean(as.numeric(sizes_ret_y),na.rm=T)
      sizes_harv_y<-sizes_harv_y[!is.na(sizes_harv_y)]
      meansize_harv[y]<-mean(as.numeric(sizes_harv_y),na.rm=T)
      }
      ##------------------------------------------------------## add to data
      data$meanSesc<-meansize_esc
      data$meanSret<-meansize_ret
      data$meanSharv<-meansize_harv

      ##==================================================================##
      ##=======================## change in reproductive potential over time
      ##==================================================================##
      use_ny<-nyh ## all years of historical trends
      ny_obs<-dim(dataObs)[1]
      ##------------------------------------## simulated escapement age comp
      esc_by_age_M<-esc_by_age_sex[,,which(names(esc_by_age_sex[1,1,])=="M")]
      esc_by_age_F<-esc_by_age_sex[,,which(names(esc_by_age_sex[1,1,])=="F")]
      esc_by_age_F[is.na(esc_by_age_F)]<-0
      esc_by_age_F<-esc_by_age_F[year_index,]
      mean_SaAs<-meanSaA[year_index,] ## no observation error
      ##--------------------------------## calculate reproductive potential
      fec_by_age<-eggm_by_age<-data.frame(array(NA,dim=c(ny_obs,nage)))
      potential<-data.frame(array(NA,dim=c(ny_obs,6)))
      names(potential)<-c("Fecundity","Eggmass","AvPerFemFecund","AvPerFemEggmass", "AvPerSpawnerFecund","AvPerSpawnerEggmass")
      for(y in 1:ny_obs) {
      sizes_y<-NA
      for(a in 1:nage) {
      if(esc_by_age_F[y,a]>0) { size_new<-mean_SaAs[y,a] } else { size_new<-NA }
      rep_out_per_cap<-SizeShiftsImplicationsV2::reprod_output(size_new,allometry) ## fecundity / egg mass
      fec_by_age[y,a]<-sum(rep(exp(rep_out_per_cap[[1]]),esc_by_age_F[y,a]))
      eggm_by_age[y,a]<-sum(rep(exp(rep_out_per_cap[[2]]),esc_by_age_F[y,a]))
      sizes_y<-c(sizes_y,rep(size_new,esc_by_age_F[y,a]))
      }
      sizes_y<-sizes_y[!is.na(sizes_y)]
      potential[y,1]<-sum(exp(SizeShiftsImplicationsV2::reprod_output(sizes_y,allometry)[[1]]))
      potential[y,2]<-sum(exp(SizeShiftsImplicationsV2::reprod_output(sizes_y,allometry)[[2]]))
      potential[y,3]<-potential[y,1]/sum(esc_by_age_F[y,]) ## females
      potential[y,4]<-potential[y,2]/sum(esc_by_age_F[y,]) ## females
      potential[y,5]<-potential[y,1]/sum(esc_by_age[y,]) ## all spawners
      potential[y,6]<-potential[y,2]/sum(esc_by_age[y,]) ## all spawners
      }
      fec_by_age[is.na(fec_by_age)]<-0
      eggm_by_age[is.na(eggm_by_age)]<-0
      potential[is.na(potential)]<-0
      ##------------------------------------------------------## add to data
      dataObs<-data.frame(cbind(dataObs,round(potential,2)))
      ##----------------------------------------## calculate percent change
      change<-data.frame(array(NA,dim=c(4,3)))
      colnames(change)<-c("3yrs","5yrs","10yrs")
      rownames(change)<-c("AvPerFemFecund","AvPerFemEggmass", "AvPerSpawnerFecund","AvPerSpawnerEggmass")
      nyrefs<-c(3,5,10)
      for(i in 1:4) {
      if(i==1) trend<-dataObs$AvPerFemFecund[1:use_ny]
      if(i==2) trend<-dataObs$AvPerFemEggmass[1:use_ny]
      if(i==3) trend<-dataObs$AvPerSpawnerFecund[1:use_ny]
      if(i==4) trend<-dataObs$AvPerSpawnerEggmass[1:use_ny]
      for(h in 1:3){
      nyref<-nyrefs[h]
      earlymean<-mean(trend[1:nyref])
      latemean<-mean(trend[(use_ny-nyref+1):use_ny])
      change[i,h]<-round((latemean-earlymean)*100/earlymean,2)
      }
      }

      ##==================================================================##
      ##======================================================## save output
      ##==================================================================##
      output.list<-list("para"=srparms,"sr_sim"=sr_sim,"fec"=change[3,], "egg"=change[4,],"S_msy"=S_msy_estimate,"data"=data,"obs"=dataObs, "ret_by_age"=ret_by_age,"meanSaA"=meanSaA,"propfemale"=propfemale, "selectivities_by_age"=selectivities_by_age,"MSY_Goals"=MSY_Goals, "impl_errors"=impl_errors)
      ##----------------------------------------------------## return lists
      return(output.list)
} ## end function

##==================================================================##
##==================================================================##
##==================================================================##
