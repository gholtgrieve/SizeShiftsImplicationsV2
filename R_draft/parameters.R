## List of input parameters for simulation model
parms.list<-list(    				
seednum=seednum, ## seed
j=j, ## scenario
k=k, ## iteration
print=F, ## print progress for diagnostics 
plot=F, ## plot data during model run (only first few iterations)
sim_recruits="eggmass", ##  from 'spawners' 'fecundity' 'eggmass'
est_method="GLS", ## SR fit method ('GLS' - previously 'JAGS','STAN')
ricker_type="const_beta", ## 'const_beta' or 'const_rmax'
var="both", ## vary 'alpha' or 'beta' or 'both' in DLM?
optimize="logFmax", ## optimize F as 'Fmax' or logFmax'
nyi=10, ## initial years not analyzed (to produce return)
nyh=50, ## number of 'historical' years with trends
ny=110, ## number of total years with trends (ny-10 reconstructed)
futureT=scen$futureT[j], ## continue trends beyond historical period?
obserr=0.2, ## observation error of escapement (default=0.2)
hobserr=0.1, ## observation error of harvest (default=0.1)
harverr=0.15, ## harvest implementation error (default=0.3) 
reglength=20, ## inverse probability of regime shift
regstr=1, ## low/high productivity factor (default=1)
harvmgmt=paste0(scen$mgmt)[j], ## 'smsy', 'umsy' or 'fix_harv_rate'
factorMSY=scen$factorMSY[j], ## aggressive or precautionary
harvrate=0.5, ## harvest rate (only used for 'fix_harv_rate')
goalfreq=10, ## review freq in yrs ('F'=no review)
procerr=0.35, ## recruitment process error (default=0.35)
rho=0.4, ## recruitment autocorrelation (default=0.4)
alpha_mean=5, ## mean productivity at low abundance (default=5)
beta_mean=5e-5, ## mean density effect (default=5e-5)
sr_corr=0.0, ## strength of alpha-beta correlation (default=0)
sr_parms_sd=0.2, ## stochastic variation SR parameters
maxsel=scen$maxsel[j], ## size of maximum selectivity (default=0.2)
sdsel=scen$sdsel[j], ## selectivity (see Bromaghin 2005)
ages=seq(1,9,1), ## all ages in model (total/BY age)
meanageini=5.5, ## initial mean age (default=5.5)
agetrend=-0.4*scen$ageT[j], ## mean age trend (default=-0.4)
sdage=0.6, ## variation in cohort age composition (default=0.6)
vonB_Linf=1200, ## von B average asymptotic length (default=1200)
vonB_k=0.325, ## von B growth rate coefficient (default=0.325)
ocean0s=150, ## size of ocean-0 fish (default=150)
sdSaA=0.01, ## sd of size-at-age anomalies (default=0.01)
sizetrends=c(0,0,30,10,-30,-60,rep(-90,3))*scen$sizeT[j], ## trends
propF=0.45, ## proportion female in return (default=0.45)
propFtrend=-0.4*scen$sexT[j], ## trend on logit scale (default=0.4)
agediff=1, ## age difference between sexes (default= 1)
allometry=c(size_ref=800,
            fec_ref=6600,
            b_fec=2.4,
            egg_ref=916,
            b_eggs=4.8),
alt_sr_param=c(alpha_N0=5.07,
               alpha_EASL=0.00177,
               alpha_EMASL=0.01036,
               beta_N0=8.6e-6,
               beta_EASL=2.973e-09,
               beta_EMASL=1.6967e-08) ## see Staton et al. 2021 
)
