###########################################################################
## AGE-DEPENDANT SEIRS MODEL WITH 5-YEAR AGE CLASSES USING UN DEMOG DATA ##
## Modified for Thailand by sompob@tropmedres.ac (8 May 2020)
###########################################################################
#rm(list=ls()) 
#setwd("~/Documents/MORU 2017/COVID19_modelling/CoMO/Model/TH_model")
setwd("~/OneDrive - tropmedres.ac/projects/2020/COVID-Vac/MM_VacCO19-master")
library("deSolve")
library("dplyr")


#########  INCIDENCE DATA
incdata_X<-read.csv("Thcovidcases.csv")
incdata_X[,1]<-as.Date(incdata_X[,1],"%d-%m-%y")


########## POP Structure
load('THpopstruct.RData')
load('mort_sever_default.Rda')

popstruc <- THpop %>% 
  select(age_category, pop) %>% 
  rename(agefloor = age_category) %>% 
  as.data.frame()

popbirth <- THpop %>% 
  select(age_category, birth) %>% 
  as.data.frame() # unit should be per person per day

mort <- THpop %>% 
  pull(death) # unit should be per person per day

ihr <- mort_sever_default %>% 
  select(age_category, ihr) %>% 
  as.data.frame()

ifr <- mort_sever_default %>% 
  select(age_category, ifr) %>% 
  as.data.frame()


##########   CONTACT DATA
load('THcontacts.RData')


#########    POP AGEING
# per year ageing matrix
A <- length(age_categories)
dd<-seq(1:A)/seq(1:A)
ageing <- t(diff(diag(dd),lag=1)/(5*365.25))
ageing<-cbind(ageing,0*seq(1:A)) # no ageing from last compartment


#########   INITIALISE SIMULATION/INTERVENTION START TIMES
startdate<-as.Date("2020-02-15") 
# stopdate<-Sys.Date() # today
stopdate<-as.Date("2022-12-31") #defult ("2020-12-31)
# stopdate<-as.Date("2020-03-18")
day_start <- as.numeric(startdate-startdate)
day_stop <- as.numeric(stopdate-startdate)
times <- seq(day_start, day_stop)

tin<-as.numeric(startdate-as.Date("2020-01-01"))/365.25
initP<-sum(popstruc[,2])       # population size 
ageindcase<-20                 # age of index case (years)
aci <- floor((ageindcase/5)+1) # age class of index case


######    THESE ARE JUST VARIABLE DEFINITIONS - PLEASE DO NOT CHANGE   #################################
######     GO TO LINE 652 TO Choose Interventions' start dates         #################################

# date to start the self isolation intervention
date_selfis_on<-as.Date("2121-12-15")
# date to start the social distancing intervention
#date_dist_on<-as.Date("2020-03-17")
date_dist_on<-as.Date("2120-6-1")
# date to start the handwashing intervention
# date_hand_on<-as.Date("2020-02-01")
date_hand_on<-as.Date("2121-12-31")
# date to start the working from home
#date_work_on<-as.Date("2020-03-19")
date_work_on<-as.Date("2121-12-15")
# date to start the school closure
#date_school_on<-as.Date("2020-03-23")
date_school_on<-as.Date("2121-12-15")
# date to start cocooning the elderly
date_cocoon_on<-as.Date("2121-12-14")
# date to start international travel ban
date_travelban_on<-as.Date("2121-12-31")
# date to start screening
date_screen_on<-as.Date("2121-12-21")
# date to start voluntary quarantine
#date_quarantine_on<-as.Date("2020-03-15")
date_quarantine_on<-as.Date("2121-12-19")
# date to start lockdown low 

date_lockdown_low_on<-as.Date("2121-12-20")
# date to start lockdown high 
date_lockdown_high_on<-as.Date("2121-12-20")

# date to start lockdown mid 
date_lockdown_mid_on<-as.Date("2020-03-23")
#date_lockdown_mid_on<-as.Date("2021-12-18")

# date to start vaccinating ("2021-12-31")
date_vaccine_on<-as.Date("2020-03-1")

#Define vaccination strategy based on 21 age group
AgeAll<-rep(1,21)
AgeHighI<-rep(0,21)
AgeHighI[5:8]<-1    #Vac age 20-39
AgeAdult<-rep(1,21)
AgeAdult[1:3]<-0    #Vac age >15
AgeElder<-rep(0,21)  
AgeElder[13:21]<-1  #Vac age >60

#choose Vaccination strategy
AgeVac<-rep(0,21)
AgeVac<-AgeAll
AgeVac<-AgeHighI

################   DEFINE PARAMETERS
parameters <- c(
  # Transmission instrinsic
  p=0.042/2,           # probabilty of infection given a contact 
  rho = 25,          # relative infectiousness of incubation phase (%) min 0 max 100 step 0.5 
  omega=5,         # average duration of immunity (years) min 0.5 max 100 step 0.5  (default 200)
  omegav=5,        # average duration of vaccinated immunity (years) min 0.5 max 100 step 0.5
  gamma=3.5,         # average incubation period (days) min 1 max 7 step 0.5 
  nui=4.5,             # average duration of symptomatic infection period (days) min 1 max 7 step 0.5
  report=0,          # percentage of all asymptomatic infections that are reported (%) min 0 max 100 step 1
  reportc=100,         # percentage of all symptomatic infections that are reported (%) min 0 max 100 step 1
  reporth=100,        # percentage of all infections requiring hospitalisation that are actually admitted to hospital (%) min 0 max 100 step 1
  beds_available = sum(popstruc[,2])*2.54/1000,#80000, # maximum number of hospital beds - numeric 
  icu_beds_available = sum(popstruc[,2])*6.6/10000,#8000, # maximum number of hospital beds - numeric 
  ventilators_available = 10000, # maximum number of ventilators - numeric
  give = 85 ,        # system capacity stressor  
  pdeath_h = 30,     # probability of dying when hospitalised 
  pdeath_hc = 35,    # probability of dying when denied hospitalisation 
  pdeath_icu = 50,   # probability of dying when admitted to ICU 
  pdeath_icuc = 75,  # probability of dying when admission to ICU denied 
  pdeath_vent = 75,  # probability of dying when ventilated 
  pdeath_ventc = 85, # probability of dying when ventilator denied 
  ihr_scaling = 1,   # scaling factor for infection hospitalisation rate
  nus = 10,          # duration of non-fatal hospitalised infection (days) min 1 max 20 step 0.5
  nusc = 10,         # duration of non-fatal denied hospitalisation infection (days) min 1 max 20 step 0.5
  nu_icu = 10,       # duration of non-fatal icu infection (days) min 1 max 20 step 0.5
  nu_icuc = 10,      # duration of non-fatal denied icu infection (days) min 1 max 20 step 0.5
  nu_vent = 10,      # duration of non-fatal ventilated infection (days) min 1 max 20 step 0.5
  nu_ventc = 10,     # duration of non-fatal denied ventilation infection (days) min 1 max 20 step 0.5
  rhos= 5,           # relative level of contacts from severely ill patients (%) min 0 max 100 step 1
  amp=0,            # relative amplitude of seasonal forcing (%) min 0 max 100 step 1
  phi=12,            # month of peak in seasonal forcing
  pclin=15,          # probability upon infection of developing clinical symptoms
  prob_icu = 70,     # probability upon hospitalisation of requiring icu admission   
  prob_vent = 80,    # probability upon admission to the UCI of requiring a ventilator
  # INTERVENTIONS
  # self isolation
  selfis_on=as.numeric(date_selfis_on-startdate),
  selfis_dur=12,    # duration of self-isolation protocol (weeks) min 1 max 52 step 1
  selfis_cov=50,    # coverage of self isolation (%) min 0 max 100 step 1
  selfis_eff=50,    # adherence to self isolation (%) min 0 max 100 step 1
  # social distancing
  dist_on=as.numeric(date_dist_on-startdate),
  dist_dur=520,      # duration of social distancing protocol (weeks) min 1 max 52 step 1
  dist_cov=80,      # coverage of social distancing (%) min 0 max 100 step 1
  dist_eff=100,     # adherence to social distancing (%) min 0 max 100 step 1
  # hand washing
  hand_on=as.numeric(date_hand_on-startdate),
  hand_dur=8,      # duration of increased hand hygiene protocol (weeks) min 1 max 52 step 1
  hand_eff=80,       # efficacy of hand hygiene  (%) min 0 max 100 step 1 (default 5)
  # working at home
  work_on=as.numeric(date_work_on-startdate),
  work_dur=8,      # duration of working from home protocol (weeks) min 1 max 52 step 1
  work_cov=50,      # coverage of working from home (%) min 0 max 100 step 1
  work_eff=85,      # efficacy of working from home (%) min 0 max 100 step 1
  w2h = 10,         # work contacts that get attibuted to home when working from home (%) min 0 max 100 step 1
  # school closures
  school_on=as.numeric(date_school_on-startdate),
  school_dur=12,    # duration of school closure (weeks) min 1 max 52 step 1
  school_eff=90,    # efficacy of school closure (%) min 0 max 100 step 1
  s2h = 20,         # school contacts that get attibuted to home when school closes (%) min 0 max 100 step 1
  # cocooning the elderly
  cocoon_on = as.numeric(date_cocoon_on-startdate), 
  cocoon_dur=16,    # duration of elderly cocoon protocol (weeks) min 1 max 52 step 1
  cocoon_eff=35,    # efficacy of elderly cocoon (%) min 0 max 100 step 1
  cocoon_cov=75,    # coverage of elderly cocoon (%) min 0 max 100 step 1
  age_cocoon=70,    # minimum age for elderly cocoon min 0 max 100 step 5
  # vaccination
  vaccine_on= as.numeric(date_vaccine_on-startdate),
  vaccine_eff1=10,   # vaccine efficacy (%)- min 0 max 100 step 1; reduce infection
  vaccine_eff2=0,   # vaccine efficacy (%)- min 0 max 100 step 1; reduce transmission
  vaccine_eff3=0,   # vaccine efficacy (%)- min 0 max 100 step 1; reduce severity
  vaccine_cov=30,    # vaccine coverage (%)- min 0 max 100 step 1
  vac_campaign=2, # Number of weeks it takes to reach maximum coverage - min 1 max 8 step 1
  
  
  # imported cases 
  mean_imports = 1,           # user defined - mean number of infectious migrants per day (number) - min 0 max 500 step 1
  travelban_on= as.numeric(date_travelban_on-startdate),
  travelban_dur = 16,         # duration of internation travel restrictions (weeks) - min 1 max 52 step 1
  travelban_eff=50,           # travel restriction efficacy (%) - min 0 max 100 step 1
  # screening - increases the rate of isolation of infectious people in the model
  screen_on = as.numeric(date_screen_on-startdate), 
  screen_dur = 12,            # duration of intensified screening (week) - min 1 max 52 step 1
  screen_cov = 90,            # sensitivity of screening test min 25 max 100 step 1
  screen_overdispersion = 4,  # overdispersion of cases around index case. If  1 likelihood same as general population min 1 max 5 step 0.2 
  screen_contacts = 4,        # number of contacts screened per index case min 1 max 10 step 1
  # quarantine - This is the bi-product of increasing testing of suspected cases with a certain false positivity rate and voluntary home quarantining of people sharing a house with an infectious case
  quarantine_on = as.numeric(date_quarantine_on-startdate),
  quarantine_cov = 70,        # coverage of quarantine (%)- min 0 max 100 step 1
  quarantine_dur = 24,        # duration of quarantine (weeks) - min 1 max 52 step 1
  quarantine_days = 14,       # days in isolation for average person (days)  - min 5 max 21 step 1
  quarantine_effort = 2,      # days to implement maximum quarantine coverage - min 1 max 5
  quarantine_eff_home = 50,   # increase in the number of contacts at home when quarantined (%) - min 0 max 100 step 52
  quarantine_eff_other = 90,  # reduction in the number of other contacts when quarantined (%) - min 0 max 100 step 52
  # lockdown
  lockdown_low_on=as.numeric(date_lockdown_low_on-startdate),
  lockdown_low_dur = 16,
  lockdown_mid_on=as.numeric(date_lockdown_mid_on-startdate),
  lockdown_mid_dur = 16,
  lockdown_high_on=as.numeric(date_lockdown_high_on-startdate),
  lockdown_high_dur = 16,
  # mean household size
  household_size = 2          # mean household size (number) - min 1 max 10 step 1 
)
######

####################
# Scale parameters to percentages/ rates
parameters["rho"]<-parameters["rho"]/100
parameters["omega"]<-(1/(parameters["omega"]*365))
parameters["omegav"]<-(1/(parameters["omegav"]*365))
parameters["gamma"]<-1/parameters["gamma"]
parameters["nui"]<-1/parameters["nui"]
parameters["report"]<-parameters["report"]/100
parameters["reportc"]<-parameters["reportc"]/100
parameters["reporth"]<-parameters["reporth"]/100
parameters["nus"]<-1/parameters["nus"]
parameters["rhos"]<-parameters["rhos"]/100
parameters["amp"]<-parameters["amp"]/100
parameters["selfis_dur"]<-parameters["selfis_dur"]*7
parameters["selfis_cov"]<-parameters["selfis_cov"]/100
parameters["selfis_eff"]<-parameters["selfis_eff"]/100
parameters["dist_dur"]<-parameters["dist_dur"]*7
parameters["dist_cov"]<-parameters["dist_cov"]/100
parameters["dist_eff"]<-parameters["dist_eff"]/100
parameters["hand_dur"]<-parameters["hand_dur"]*7
parameters["hand_eff"]<-parameters["hand_eff"]/100
parameters["work_dur"]<-parameters["work_dur"]*7
parameters["work_cov"]<-parameters["work_cov"]/100
parameters["work_eff"]<-parameters["work_eff"]/100
parameters["w2h"]<-parameters["w2h"]/100
parameters["school_dur"]<-parameters["school_dur"]*7
parameters["schoolcov"]<-parameters["schoolcov"]/100
parameters["school_eff"]<-parameters["school_eff"]/100
parameters["s2h"]<-parameters["s2h"]/100
parameters["cocoon_dur"]<-parameters["cocoon_dur"]*7
parameters["cocoon_cov"]<-parameters["cocoon_cov"]/100
parameters["cocoon_eff"]<-parameters["cocoon_eff"]/100
parameters["age_cocoon"]<-floor((parameters["age_cocoon"]/5)+1)
parameters["travelban_eff"]<-parameters["travelban_eff"]/100
parameters["vaccine_eff1"]<-parameters["vaccine_eff1"]/100
parameters["vaccine_eff2"]<-parameters["vaccine_eff2"]/100
parameters["vaccine_eff3"]<-parameters["vaccine_eff3"]/100
parameters["vaccine_cov"]<-parameters["vaccine_cov"]/100
parameters["vac_campaign"]<-parameters["vac_campaign"]*7
parameters["travelban_dur"]<-parameters["travelban_dur"]*7
parameters["screen_dur"]<-parameters["screen_dur"]*7
parameters["screen_cov"]<-parameters["screen_cov"]/100
parameters["quarantine_cov"]<-parameters["quarantine_cov"]/100
parameters["quarantine_dur"]<-parameters["quarantine_dur"]*7
parameters["quarantine_days"]<-parameters["quarantine_days"]
parameters["quarantine_effort"]<-1/parameters["quarantine_effort"]
parameters["quarantine_eff_home"]<-parameters["quarantine_eff_home"]/-100
parameters["quarantine_eff_other"]<-parameters["quarantine_eff_other"]/100
parameters["give"]<-parameters["give"]/100
parameters["pdeath_h"]<-parameters["pdeath_h"]/100
parameters["pdeath_hc"]<-parameters["pdeath_hc"]/100
parameters["pdeath_icu"]<-parameters["pdeath_icu"]/100
parameters["pdeath_icuc"]<-parameters["pdeath_icuc"]/100
parameters["pdeath_vent"]<-parameters["pdeath_vent"]/100
parameters["pdeath_ventc"]<-parameters["pdeath_ventc"]/100
parameters["nusc"]<-1/parameters["nusc"]
parameters["nu_icu"]<-1/parameters["nu_icu"]
parameters["nu_icuc"]<-1/parameters["nu_icuc"]
parameters["nu_vent"]<-1/parameters["nu_vent"]
parameters["nu_ventc"]<-1/parameters["nu_ventc"]
parameters["pclin"]<-parameters["pclin"]/100
parameters["prob_icu"]<-parameters["prob_icu"]/100
parameters["prob_vent"]<-parameters["prob_vent"]/100
parameters["lockdown_low_dur"]<-parameters["lockdown_low_dur"]*7
parameters["lockdown_mid_dur"]<-parameters["lockdown_mid_dur"]*7
parameters["lockdown_high_dur"]<-parameters["lockdown_high_dur"]*7
#########
# parameters2<-parameters
# 
# #########    SEVERITY AND MORTALITY
# # age dependent hosp and mort - correction parameters allow to change the shape of the IFR curve
# ifr_correction_young<-2
# ifr_correction_old<-1.75
# ihr <- read.csv("covidagehosp_X.csv",header=TRUE) # hospitalisation rate given infection
# ifr <- read.csv("covidagefrpi_X.csv",header=TRUE) # fatality rate given infection
# ihr<- parameters["ihr_scaling"]*ihr/100   # csv data is in percentages
# ifr_original<-ifr/100   # csv data is in percentages
# for (i in 1:A){
#   ifr[i,2]=ifr[i,2]/max(ifr[,2])    # transform ifr into a normalised age profile (highest value turns into 1)
# }
# ifr[1:14,2]<-ifr[1:14,2]/ifr_correction_young
# ihr$severe[15:21]<-ihr$severe[15:21]*ifr_correction_old

###########################################################################
# Define the indices for each variable
Sindex<-1:A
Eindex<-(A+1):(2*A)
Iindex<-(2*A+1):(3*A)
Rindex<-(3*A+1):(4*A)
Xindex<-(4*A+1):(5*A)
Hindex<-(5*A+1):(6*A)
HCindex<-(6*A+1):(7*A)
Cindex<-(7*A+1):(8*A)
CMindex<-(8*A+1):(9*A)

SVindex<-(9*A+1):(10*A)
EVindex<-(10*A+1):(11*A)
IVindex<-(11*A+1):(12*A)
CLVindex<-(12*A+1):(13*A)
HVindex<-(13*A+1):(14*A)
ICUVindex<-(14*A+1):(15*A)
VentVindex<-(15*A+1):(16*A)
RVindex<-(16*A+1):(17*A)
#dSVdt,dEVdt,dIVdt,dCLVdt,dHVdt,dICUVdt,dVentVdt,dRVdt

QSindex<-(17*A+1):(18*A)
QEindex<-(18*A+1):(19*A)
QIindex<-(19*A+1):(20*A)
QRindex<-(20*A+1):(21*A)
CLindex<-(21*A+1):(22*A)
QCindex<-(22*A+1):(23*A)
ICUindex<-(23*A+1):(24*A)
ICUCindex<-(24*A+1):(25*A)
Ventindex<-(25*A+1):(26*A)
VentCindex<-(26*A+1):(27*A)
CMCindex<-(27*A+1):(28*A)

###########################################################################
# MODEL INITIAL CONDITIONS
initI<-0*popstruc[,2]  # Infected and symptomatic
initE<-0*popstruc[,2]  # Incubating
initE[aci]<-1          # place random index case in E compartment
initR<-0*popstruc[,2]  # Immune
initX<-0*popstruc[,2]  # Isolated 
#initV<-0*popstruc[,2]  # Vaccinated 
initSV<-0*popstruc[,2]  # Vaccinated and susceptible
initEV<-0*popstruc[,2]  # Vaccinated cases
initIV<-0*popstruc[,2]  # Vaccinated asymptomatic cases
initCLV<-0*popstruc[,2]  # Vaccinated symptomatic cases
initHV<-0*popstruc[,2]  # Vaccinated hospitalised cases
initICUV<-0*popstruc[,2]  # Vaccinated hospi ICU cases
initVentV<-0*popstruc[,2]  # Vaccinated hosp ICU with ventilator
initRV<-0*popstruc[,2]  # Vaccinated infected and recovered

initQS<-0*popstruc[,2] # quarantined S 
initQE<-0*popstruc[,2] # quarantined E  
initQI<-0*popstruc[,2] # quarantined I  
initQR<-0*popstruc[,2] # quarantined R  
initH<-0*popstruc[,2]  # hospitalised 
initHC<-0*popstruc[,2] # hospital critical 
initC<-0*popstruc[,2]  # Cumulative cases (true)
initCM<-0*popstruc[,2] # Cumulative deaths (true)
initCL<-0*popstruc[,2] # symptomatic cases
initQC<-0*popstruc[,2] # quarantined C 
initICU<-0*popstruc[,2]   # icu
initICUC<-0*popstruc[,2]  # icu critical
initVent<-0*popstruc[,2]  # icu vent
initVentC<-0*popstruc[,2] # icu vent crit
initCMC<-0*popstruc[,2]   # Cumulative deaths (true)
initS<-popstruc[,2]-initE-initI-initR-initX-initSV-initEV-initIV-initCLV-initHV-initICUV-initVentV-initRV-initH-initHC-initQS-initQE-initQI-initQR-initCL-initQC-initICU-initICUC-initVent-initVentC  # Susceptible (non-immune)

###############
# set progress bar
pb <- txtProgressBar(min = 0, max = length(times), style = 3)

# set up a function to solve the equations
covid<-function(t, Y, parameters) 
{
  
  with(as.list(c(Y, parameters)),
       {
         S <- Y[Sindex]
         E <- Y[Eindex]
         I <- Y[Iindex]
         R <- Y[Rindex]
         X <- Y[Xindex]
         H <- Y[Hindex]
         HC <- Y[HCindex]
         C <- Y[Cindex]
         CM <- Y[CMindex]
         #V <- Y[Vindex]
         SV <- Y[SVindex]
         EV <- Y[EVindex]
         IV <- Y[IVindex]
         CLV <- Y[CLVindex]
         HV <- Y[HVindex]
         ICUV <- Y[ICUVindex]
         VentV <- Y[VentVindex]
         RV <- Y[RVindex]
         QS <- Y[QSindex]
         QE <- Y[QEindex]
         QI <- Y[QIindex]
         QR <- Y[QRindex]
         CL <- Y[CLindex]
         QC <- Y[QCindex]
         ICU <- Y[ICUindex]
         ICUC <- Y[ICUCindex]
         Vent <- Y[Ventindex]
         VentC <- Y[VentCindex]
         CMC <- Y[CMCindex]
         P <- (S+E+I+R+X+SV+EV+IV+CLV+HV+ICUV+VentV+RV+H+HC+QS+QE+QI+QR+CL+QC+ICU+ICUC+Vent+VentC)
         # print(sum(P))
         
         # health system performance
         f <- c(1,(1+give)/2,(1-give)/2,0)
         KH<-beds_available
         KICU<- icu_beds_available
         Kvent<- ventilators_available
         x.H <- c(0,(1+give)*KH/2,(3-give)*KH/2,2*KH)
         x.ICU <- c(0,(1+give)*KICU/2,(3-give)*KICU/2,2*KICU)
         x.Vent <- c(0,(1+give)*Kvent/2,(3-give)*Kvent/2,2*Kvent)
         fH <- splinefun(x.H, f, method = "hyman")
         fICU <- splinefun(x.ICU, f, method = "hyman")
         fVent<- splinefun(x.Vent, f, method = "hyman")
         critH<-min(1-fH(sum(H)+sum(ICUC))+(1-reporth),1)
         crit<-min(1-fICU(sum(ICU)+sum(Vent)+sum(VentC)),1)
         critV<-min(1-fVent(sum(Vent)),1)
         # print(fH(sum(H)))
         
         # interventions
         isolation<-(t>=selfis_on)*(t<=selfis_on+selfis_dur)
         distancing<-(t>=dist_on)*(t<=(dist_on+dist_dur))
         handwash<-(t>=hand_on)*(t<=(hand_on+hand_dur))
         workhome<-(t>=work_on)*(t<=(work_on+work_dur))
         schoolclose<-(t>=school_on)*(t<=(school_on+school_dur))
         cocoon<-(t>=cocoon_on)*(t<=(cocoon_on+cocoon_dur))*cocoon_cov
         vaccine<-(t>=(vaccine_on))*(t<=vaccine_on+vac_campaign)
         
         
         travelban<-(t>=travelban_on)*(t<=(travelban_on+travelban_dur))
         screen<-(t>=screen_on)*(t<=(screen_on+screen_dur))
         quarantine<-(t>=quarantine_on)*(t<=(quarantine_on+quarantine_dur))
         lockdown_low<-(t>=lockdown_low_on)*(t<=(lockdown_low_on+lockdown_low_dur))
         lockdown_mid<-(t>=lockdown_mid_on)*(t<=(lockdown_mid_on+lockdown_mid_dur))
         lockdown_high<-(t>=lockdown_high_on)*(t<=(lockdown_high_on+lockdown_high_dur))

         screen_eff<-0
         selfis<-0
         school<-1
         dist<-1
         hand<-0
         vaccinate<-0
         trvban_eff<-0
         quarantine_rate<-0
         
         if (lockdown_low || lockdown_mid || lockdown_high){
           if(lockdown_low){
             selfis<-0.5
             dist<-0.25
             school<-0
             trvban_eff<-0
             quarantine_rate<-0
             work<-0
             cocoon<-0.95
             hand<-0.05
             vaccinate<-0
           }
           if(lockdown_mid){
             selfis<-0.5
             dist<-0.35
             school<-0.85
             trvban_eff<-0
             quarantine_rate<-0.05
             work<-0.5
             cocoon<-0.95
             hand<-0.05
             vaccinate<-0
           }
           if(lockdown_high){
             selfis<-0.95
             dist<-0.95
             school<-0.85
             trvban_eff<-0.95
             quarantine_rate<-0.9
             work<-0.75
             cocoon<-0.95
             hand<-0.075
             vaccinate<-0
           }
         }
         else{
           if (workhome){
             work<-work_cov*work_eff
           }else{work<-1}
           if (isolation){
             selfis<-selfis_cov
             if(screen){
               screen_eff<-min((report*I+reportc*(CL)+H+ICU+Vent+reporth*HC+ICUC+VentC)*screen_contacts*(screen_overdispersion*I/P)*screen_cov/P,1) 
             }
           }
           if (schoolclose){
             school<-school_eff
           }
           if(distancing){
             dist<-dist_cov*dist_eff
           }
           if(handwash){
             hand<-hand_eff
           }
           if(vaccine){
             vac_rate <- (-log(1-vaccine_cov)/vac_campaign)
             vaccinate <- vac_rate
           }

           if(travelban){
             trvban_eff<-travelban_eff
           }
           if(quarantine){
             quarantine_rate<-min(((I+CL+H+ICU+Vent+HC+ICUC+VentC)*(household_size-1)/P),1)*quarantine_cov*quarantine_effort
           }
         }
         
         
         # cocooning the elderly
         cocoon_mat<-matrix((1-cocoon_eff),nrow = length(popstruc$pop),ncol = length(popstruc$pop))
         cocoon_mat[1:(age_cocoon-1),1:(age_cocoon-1)]<-1
         
         # contact matrices
         cts<-(contact_home+distancing*(1-dist)*contact_other+(1-distancing)*contact_other
               +(1-schoolclose)*contact_school # school on
               +schoolclose*(1-school)*contact_school # school close
               +schoolclose*contact_home*school*s2h # inflating contacts at home when school closes
               +(1-workhome)*contact_work  # normal work
               +workhome*(1-work)*contact_work # people not working from home when homework is active
               +contact_home*workhome*work*w2h # inflating contacts at home when working from home
         )
         
         # Final transmission related parameters
         contacts <- (1-cocoon)*cts+cocoon*cts*cocoon_mat+cocoon*(1+school*(1-school_eff)+work*(1-work_eff))*contact_home*(1-cocoon_mat)
         seas <- 1+amp*cos(2*3.14*(t-(phi*365.25/12))/365.25)
         importation <- mean_imports*(1-trvban_eff)
         HH<-H+ICU+Vent
         HHC<-HC+ICUC+VentC
         lam <- (1-hand)*p*seas*(contacts%*%((rho*E+(I+CL+importation)+(1-selfis_eff)*(X+HHC)+rhos*(HH))/P))
         # contacts under home quarantine
         lamq<-(1-hand)*p*seas*((1-quarantine_eff_home)*contact_home%*%(((1-selfis_eff)*(X+HHC))/P))+(1-hand)*p*seas*(1-quarantine_eff_other)*(contact_other%*%((rho*E+(I+CL+importation)+(1-selfis_eff)*(X+HHC)+rhos*(HH))/P))
         
         # birth/death
         b1<-sum(popbirth[,2]*popstruc[,2])
         birth<-0*popbirth[,2]
         birth[1]<-b1
         
         # ODE system
         dSdt <- -S*lam-S*vaccinate+omega*R+ageing%*%S-mort*S+birth-quarantine_rate*S +(1/quarantine_days)*QS+omegav*SV+omega*RV
         dEdt <- S*lam-gamma*E+ageing%*%E-mort*E-quarantine_rate*E+(1/quarantine_days)*QE 
         dIdt <- gamma*(1-pclin)*(1-screen_eff)*(1-ihr[,2])*E-nui*I+ageing%*%I-mort*I + (1/quarantine_days)*QI - quarantine_rate*I
         dCLdt<- gamma*pclin*(1-selfis)*(1-ihr[,2])*E-nui*CL+ageing%*%CL-mort*CL + (1/quarantine_days)*QC
         dRdt <- nui*I-omega*R+nui*X+nui*CL+ageing%*%R-mort*R + (1/quarantine_days)*QR + nus*(1-pdeath_h*ifr[,2])*H + (1-pdeath_icu*ifr[,2])*nu_icu*ICU + (1-pdeath_icuc*ifr[,2])*nu_icuc*ICUC + (1-pdeath_hc*ifr[,2])*nusc*HC + (1-pdeath_vent*ifr[,2])*nu_vent*Vent+ (1-pdeath_ventc*ifr[,2])*nu_ventc*VentC - vaccinate*R + omegav*RV
         dXdt <- gamma*selfis*pclin*(1-ihr[,2])*E+gamma*(1-pclin)*screen_eff*(1-ihr[,2])*E-nui*X+ageing%*%X-mort*X 
         ############
         #dVdt <- vaccinate*S -(1-vaccine_eff)*lam*V +ageing%*%V - mort*V
         
         #to add terms (vaccine_eff1=%reduction in infection, vaccine_eff2=%reduction in duration of infection, vaccine_eff3=%reduction in risk of sevrerity, hospitalisation in this case, but we have hospi, ICU and ICUVent [need to clarify the effect of vac])
         #SV, EV, IV, CL,HV,ICUV, VentV, RV set initial condition
         #
         #if (vaccine){cf
         #Add vaccine compartment - add (AgeVac*)S*vaccinate to dSdt and dSVdt
         dSVdt <- S*vaccinate - (1-vaccine_eff1)*SV*lam + ageing%*%SV-mort*SV-omegav*SV #Assuming the lam is the same as general population
         
         dEVdt <- (1-vaccine_eff1)*SV*lam - gamma*EV +ageing%*%EV-mort*EV
         dIVdt <- gamma*(1-pclin)*(1-ihr[,2]*(1+vaccine_eff3))*EV-nui*(1+vaccine_eff2)*IV+ageing%*%IV - mort*IV
         dCLVdt<- gamma*(pclin)*(1-ihr[,2]*(1+vaccine_eff3))*EV-nui*(1+vaccine_eff2)*CLV+ageing%*%CLV - mort*CLV
         dHVdt <- gamma*ihr[,2]*(1+vaccine_eff3)*(1-prob_icu)*EV-nus*(1+vaccine_eff2)*HV+ ageing%*%HV - mort*HV
         dICUVdt <- gamma*ihr[,2]*(1-vaccine_eff3)*prob_icu*(1-prob_vent)*EV-nu_icu*(1+vaccine_eff2)*ICUV +ageing%*%ICUV - mort*ICUV
         dVentVdt <- gamma*ihr[,2]*(1-vaccine_eff3)*prob_icu*prob_vent*EV-nu_vent*(1+vaccine_eff2)*VentV +ageing%*%VentV - mort*VentV
         dRVdt <- nui*(1+vaccine_eff2)*IV + nui*(1+vaccine_eff2)*CLV +nus*(1+vaccine_eff2)*HV+nu_icu*(1+vaccine_eff2)*ICUV+nu_vent*(1+vaccine_eff2)*VentV + ageing%*%RV-mort*RV
         + nus*(1-pdeath_h*ifr[,2])*HV + (1-pdeath_icu*ifr[,2])*nu_icu*ICU + (1-pdeath_icuc*ifr[,2])*nu_icuc*ICUC + (1-pdeath_hc*ifr[,2])*nusc*HC + (1-pdeath_vent*ifr[,2])*nu_vent*Vent+ (1-pdeath_ventc*ifr[,2])*nu_ventc*VentC-omega*RV + vaccinate*R - omegav*RV
         #nui*I-omega*R+nui*X+nui*CL+ageing%*%R-mort*R + (1/quarantine_days)*QR + nus*(1-pdeath_h*ifr[,2])*H + (1-pdeath_icu*ifr[,2])*nu_icu*ICU + (1-pdeath_icuc*ifr[,2])*nu_icuc*ICUC + (1-pdeath_hc*ifr[,2])*nusc*HC + (1-pdeath_vent*ifr[,2])*nu_vent*Vent+ (1-pdeath_ventc*ifr[,2])*nu_ventc*VentC
        # from RV back to S   
         #} else {
        #   dSVdt<-0*SV
         #   dEVdt<-0*EV
         #   dIVdt<-0*IV
         #   dCLVdt<-0*CLV
         # dHVdt<-0*HV
         # dICUVdt<-0*ICUV
         # dVentVdt<-0*VentV
         # dRVdt<-0*RV
         #}
         #modified to include Cumulative cases and deaths for those with vaccination
         
         
         ############ 
         dQSdt <- quarantine_rate*S+ ageing%*%QS-mort*QS - (1/quarantine_days)*QS - lamq*QS
         dQEdt <- quarantine_rate*E - gamma*QE + ageing%*%QE-mort*QE - (1/quarantine_days)*QE + lamq*QS 
         dQIdt <- quarantine_rate*I + gamma*(1-ihr[,2])*(1-pclin)*QE-nui*QI+ageing%*%QI-mort*QI - (1/quarantine_days)*QI
         dQCdt <- gamma*(1-ihr[,2])*pclin*QE-nui*QC+ageing%*%QC-mort*QC - (1/quarantine_days)*QC
         dQRdt <- nui*QI+nui*QC+ageing%*%QR-mort*QR - (1/quarantine_days)*QR
         
         dHdt <- gamma*ihr[,2]*(1-prob_icu)*(1-critH)*E + gamma*ihr[,2]*(1-prob_icu)*(1-critH)*QE - nus*H + ageing%*%H-mort*H  # all pdeath have to be lower than
         dHCdt <- gamma*ihr[,2]*(1-prob_icu)*critH*E + gamma*ihr[,2]*(1-prob_icu)*critH*QE - nusc*HC + ageing%*%HC-mort*HC 
         dICUdt <- gamma*ihr[,2]*prob_icu*(1-crit)*(1-prob_vent)*E + gamma*ihr[,2]*prob_icu*(1-crit)*(1-prob_vent)*QE - nu_icu*ICU +ageing%*%ICU - mort*ICU 
         dICUCdt <- gamma*ihr[,2]*prob_icu*crit*(1-prob_vent)*E + gamma*ihr[,2]*prob_icu*crit*(1-prob_vent)*QE - nu_icuc*ICUC +ageing%*%ICUC - mort*ICUC 
         dVentdt <- gamma*ihr[,2]*prob_icu*(1-crit)*(1-critV)*prob_vent*E + gamma*ihr[,2]*prob_icu*(1-crit)*(1-critV)*prob_vent*QE + (1-critV)*VentC*1/2 - nu_vent*Vent +ageing%*%Vent - mort*Vent 
         dVentCdt <- gamma*ihr[,2]*prob_icu*prob_vent*(1-crit)*critV*E +gamma*ihr[,2]*prob_icu*prob_vent*crit*E+
           gamma*ihr[,2]*prob_icu*prob_vent*(1-crit)*critV*QE + gamma*ihr[,2]*prob_icu*prob_vent*crit*QE - 
           (1-critV)*VentC*1/2-nu_ventc*VentC +ageing%*%VentC - mort*VentC
         
         #Add Terms on accumulate case and death from vaccine compartments       
         dCdt <- report*gamma*(1-pclin)*(1-ihr[,2])*(E+QE)+reportc*gamma*pclin*(1-ihr[,2])*(E+QE)+ 
           gamma*ihr[,2]*(1-critH)*(1-prob_icu)*(E+QE)+gamma*ihr[,2]*critH*reporth*(1-prob_icu)*(E+QE)+gamma*ihr[,2]*prob_icu*(E+QE)+
           report*gamma*(1-pclin)*(1-ihr[,2])*(1-vaccine_eff3)*EV+reportc*gamma*(pclin)*(1-ihr[,2]*(1-vaccine_eff3))*EV+
           gamma*ihr[,2]*(1-vaccine_eff3)*(1-prob_icu)*EV+gamma*ihr[,2]*(1-vaccine_eff3)*prob_icu*(1-prob_vent)*EV+gamma*ihr[,2]*(1-vaccine_eff3)*prob_icu*prob_vent*EV
         
         dCMdt<- nus*pdeath_h*ifr[,2]*H + nusc*pdeath_hc*ifr[,2]*HC + nu_icu*pdeath_icu*ifr[,2]*ICU + nu_icuc*pdeath_icuc*ifr[,2]*ICUC +  nu_vent*pdeath_vent*ifr[,2]*Vent + nu_ventc*pdeath_ventc*ifr[,2]*VentC + 
           mort*H + mort*HC + mort*ICU + mort*ICUC + mort*Vent + mort*VentC+
           nus*(1-vaccine_eff2)*pdeath_h*ifr[,2]*HV + nu_icu*(1-vaccine_eff2)*pdeath_icu*ifr[,2]*ICUV + nu_vent*(1-vaccine_eff2)*pdeath_vent*ifr[,2]*VentV + mort*HV + mort*ICUV + mort*VentV
         
         dCMCdt <- nusc*pdeath_hc*ifr[,2]*HC+nu_icuc*pdeath_icuc*ifr[,2]*ICUC + nu_ventc*pdeath_ventc*ifr[,2]*VentC + 
           mort*HC + mort*ICUC + mort*VentC

         setTxtProgressBar(pb, t)
         
         # return the rate of change
         list(c(dSdt,dEdt,dIdt,dRdt,dXdt,dHdt,dHCdt,dCdt,dCMdt,dSVdt,dEVdt,dIVdt,dCLVdt,dHVdt,dICUVdt,dVentVdt,dRVdt,dQSdt,dQEdt,dQIdt,dQRdt,dCLdt,dQCdt,dICUdt,dICUCdt,dVentdt,dVentCdt,dCMCdt))
       }
  ) 
}

#modify to get the target outcomes

#continue working on model outputs including V compartment
process_ode_outcome <- function(out){
  # Start Bridge ----
  critH<-c()
  crit<-c()
  critV<-c()
  # End Bridge ----
  
  # START Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  f <- c(1,(1+parameters["give"])/2,(1-parameters["give"])/2,0) 
  KH<-parameters["beds_available"]
  KICU<- parameters["icu_beds_available"]
  Kvent<- parameters["ventilators_available"]
  x.H <- c(0,(1+parameters["give"])*KH/2,(3-parameters["give"])*KH/2,2*KH) 
  x.ICU <- c(0,(1+parameters["give"])*KICU/2,(3-parameters["give"])*KICU/2,2*KICU) 
  x.Vent <- c(0,(1+parameters["give"])*Kvent/2,(3-parameters["give"])*Kvent/2,2*Kvent) 
  fH <- splinefun(x.H, f, method = "hyman") 
  fICU <- splinefun(x.ICU, f, method = "hyman") 
  fVent<- splinefun(x.Vent, f, method = "hyman") 
  for (i in 1:length(times)){
    critH[i]<-min(1-fH(sum(out[i,(Hindex+1)]))+(1-parameters["reporth"]),1)
    crit[i]<-min(1-fICU((sum(out[i,(ICUindex+1)]))+(sum(out[i,(Ventindex+1)]))+(sum(out[i,(VentCindex+1)]))))
    critV[i]<-min(1-fVent((sum(out[i,(Ventindex+1)]))),1)
  }
  
  # total population
  pop1<-out[,(Sindex+1)]+out[,(Eindex+1)]+out[,(Iindex+1)]+out[,(CLindex+1)]+out[,(Rindex+1)]+out[,(Xindex+1)]+
    out[,(SVindex+1)]+out[,(EVindex+1)]+out[,(IVindex+1)]+out[,(CLVindex+1)]+out[,(HVindex+1)]+out[,(ICUVindex+1)]+out[,(VentVindex+1)]+out[,(RVindex+1)]+
    out[,(QSindex+1)]+out[,(QEindex+1)]+out[,(QIindex+1)]+out[,(QCindex+1)]+out[,(QRindex+1)]+
    out[,(Hindex+1)]+out[,(HCindex+1)]+out[,(ICUindex+1)]+out[,(ICUCindex+1)]+out[,(Ventindex+1)]+out[,(VentCindex+1)] 
  tpop1<-rowSums(pop1)
  time<-as.Date(out[,1]+startdate)
  # daily incidence
  inc1 <- parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*out[,(Eindex+1)]%*%(1-ihr[,2])+
    parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*out[,(Eindex+1)]%*%(1-ihr[,2])+
    parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*out[,(QEindex+1)]%*%(1-ihr[,2])+
    parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*out[,(QEindex+1)]%*%(1-ihr[,2])+
    parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*out[,(EVindex+1)]%*%(1-ihr[,2])+
    parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*out[,(EVindex+1)]%*%(1-ihr[,2])
  
  inc1h<- parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*(1-critH)*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*critH*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*critH*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*parameters["prob_icu"]+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*parameters["prob_icu"]+
    parameters["gamma"]*out[,(EVindex+1)]%*%ihr[,2]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(EVindex+1)]%*%ihr[,2]*parameters["prob_icu"]

  
  dailyinc1<-rowSums(inc1)+rowSums(inc1h)      # daily incidence
  cuminc1<-colSums(inc1)+colSums(inc1h)        # cumulative incidence
  
  previcureq0<-rowSums(out[,(Iindex+1)])+ rowSums(out[,(IVindex+1)])  #asymptomatic cases
  previcureq01<-rowSums(out[,(CLindex+1)])+ rowSums(out[,(CLVindex+1)]) #mild symp cases

  
  previcureq1<-rowSums(out[,(Hindex+1)])+ rowSums(out[,(ICUCindex+1)])+ rowSums(out[,(HVindex+1)])+ rowSums(out[,(ICUVindex+1)])    # requirement for beds
  previcureq21<-rowSums(out[,(ICUindex+1)])+rowSums(out[,(VentCindex+1)])+rowSums(out[,(ICUVindex+1)])+rowSums(out[,(VentVindex+1)])   # requirement for icu
  previcureq31<-rowSums(out[,(Ventindex+1)])+rowSums(out[,(VentVindex+1)])   # requirement for icu
  cmortality1<-rowSums(out[,(CMindex+1)])      # cumulative mortality
  overloadH1<-rowSums(out[,(HCindex+1)])       # requirement for beds
  overloadICU1<-rowSums(out[,(ICUCindex+1)])   # requirement for beds
  overloadVent1<-rowSums(out[,(VentCindex+1)]) # requirement for beds
  ccases1<-rowSums(out[,(Cindex+1)])           # cumulative cases
  
  inc_overloadH1<-((parameters["gamma"]*(1-parameters["prob_icu"])*out[,(Eindex+1)])) #not taken V compartment into account
  inc_overloadICU1<-((parameters["gamma"]*parameters["prob_icu"]*(1-parameters["prob_vent"])*out[,(Eindex+1)]))
  for (i in 1:length(times)) {
    inc_overloadH1[i,]<-inc_overloadH1[i,]*critH[i]*ihr[,2]
    inc_overloadICU1[i,]<-inc_overloadICU1[i,]*crit[i]*ihr[,2]
  }
  inc_overloadH1<-cumsum(rowSums(inc_overloadH1))
  inc_overloadICU1<-cumsum(rowSums(inc_overloadICU1))
  # END Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  
  # START Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  ##########################    CALCULATE MORTALITY 
  pdeath_hc<-parameters["pdeath_hc"]
  prob_icu<-parameters["prob_icu"]
  prob_vent<-parameters["prob_vent"]
  pdeath_icuc<-parameters["pdeath_icuc"]
  pdeath_ventc<-parameters["pdeath_ventc"]
  
  
  cinc_mort_H1 <- cumsum(rowSums(parameters["nus"]*parameters["pdeath_h"]*(out[,(Hindex+1)]%*%ifr[,2])+ out[,(Hindex+1)]%*%mort))+
    cumsum(rowSums(parameters["nus"]*parameters["pdeath_h"]*(out[,(HVindex+1)]%*%ifr[,2])+ out[,(HVindex+1)]%*%mort))
  cinc_mort_HC1 <- cumsum(rowSums(parameters["nusc"]*parameters["pdeath_hc"]*(out[,(HCindex+1)]%*%ifr[,2]) + out[,(HCindex+1)]%*%mort))
  cinc_mort_ICU1 <- cumsum(rowSums(parameters["nu_icu"]*parameters["pdeath_icu"]*out[,(ICUindex+1)]%*%ifr[,2] + out[,(ICUindex+1)]%*%mort))+
    cumsum(rowSums(parameters["nu_icu"]*parameters["pdeath_icu"]*out[,(ICUVindex+1)]%*%ifr[,2] + out[,(ICUVindex+1)]%*%mort))
  cinc_mort_ICUC1 <- cumsum(rowSums(parameters["nu_icuc"]*parameters["pdeath_icuc"]*out[,(ICUCindex+1)]%*%ifr[,2] + out[,(ICUCindex+1)]%*%mort))
  cinc_mort_Vent1 <- cumsum(rowSums(parameters["nu_vent"]*parameters["pdeath_vent"]*out[,(Ventindex+1)]%*%ifr[,2] + out[,(Ventindex+1)]%*%mort))+
    cumsum(rowSums(parameters["nu_vent"]*parameters["pdeath_vent"]*out[,(VentVindex+1)]%*%ifr[,2] + out[,(VentVindex+1)]%*%mort))
  cinc_mort_VentC1 <- cumsum(rowSums(parameters["nu_ventc"]*parameters["pdeath_ventc"]*out[,(VentCindex+1)]%*%ifr[,2] + out[,(VentCindex+1)]%*%mort))
  
  base_mort_H1 <- cumsum(rowSums(out[,(Hindex+1)]%*%mort)) + cumsum(rowSums(out[,(HVindex+1)]%*%mort))
  base_mort_HC1 <- cumsum(rowSums(out[,(HCindex+1)]%*%mort))
  base_mort_ICU1 <- cumsum(rowSums(out[,(ICUindex+1)]%*%mort)) +cumsum(rowSums(out[,(ICUVindex+1)]%*%mort))
  base_mort_ICUC1 <- cumsum(rowSums(out[,(ICUCindex+1)]%*%mort))
  base_mort_Vent1 <- cumsum(rowSums(out[,(Ventindex+1)]%*%mort))+cumsum(rowSums(out[,(VentVindex+1)]%*%mort))
  base_mort_VentC1 <- cumsum(rowSums(out[,(VentCindex+1)]%*%mort))
  
  base_mort_S1 <- cumsum(rowSums(out[,(Sindex+1)]%*%mort))+cumsum(rowSums(out[,(SVindex+1)]%*%mort))
  base_mort_E1 <-  cumsum(rowSums(out[,(Eindex+1)]%*%mort))+cumsum(rowSums(out[,(EVindex+1)]%*%mort))
  base_mort_I1 <-  cumsum(rowSums(out[,(Iindex+1)]%*%mort))+cumsum(rowSums(out[,(IVindex+1)]%*%mort))
  base_mort_CL1 <- cumsum(rowSums(out[,(CLindex+1)]%*%mort))+cumsum(rowSums(out[,(CLVindex+1)]%*%mort))
  base_mort_X1 <-  cumsum(rowSums(out[,(Xindex+1)]%*%mort))    
  base_mort_QS1 <- cumsum(rowSums(out[,(QSindex+1)]%*%mort)) 
  base_mort_QE1 <- cumsum(rowSums(out[,(QEindex+1)]%*%mort)) 
  base_mort_QI1 <- cumsum(rowSums(out[,(QIindex+1)]%*%mort)) 
  base_mort_QC1 <- cumsum(rowSums(out[,(QCindex+1)]%*%mort)) 
  base_mort_QR1 <- cumsum(rowSums(out[,(QRindex+1)]%*%mort)) 
  base_mort_R1 <-  cumsum(rowSums(out[,(Rindex+1)]%*%mort))+cumsum(rowSums(out[,(RVindex+1)]%*%mort))
  
  # END Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----

  
  # Export in a cohesive format ----
  results <- list()
  results$time <- startdate + times  # dates
  #results$Rt <- Rt
  results$cum_mortality <- round(cmortality1)  # cumulative mortality
  results$pct_total_pop_infected <- round(100 * tail(cumsum(rowSums(parameters["gamma"]*out[,(Eindex+1)])),1)/sum(popstruc[,2]), 1) +round(100 * tail(cumsum(rowSums(parameters["gamma"]*out[,(EVindex+1)])),1)/sum(popstruc[,2]), 1) # proportion of the  population that has been infected at the end of the simulation
  results$doubling_time <- round(log(2)*7 / (log(dailyinc1[2+7] / dailyinc1[2])), 2)  # (Baseline only) to double the number of infections at inception
  results$required_beds <- round(previcureq1)  # required beds
  results$saturation <- parameters["beds_available"]  # saturation
  results$daily_incidence <- round(dailyinc1)  # daily incidence (Reported)
  results$daily_total_cases <- round(rowSums(parameters["gamma"]*out[,(Eindex+1)]+parameters["gamma"]*out[,(QEindex+1)]+parameters["gamma"]*out[,(EVindex+1)])) # daily incidence (Reported + Unreported)  # daily incidence (Reported + Unreported)
  results$cum_cases<- ccases1
  results$asymp_cases <- round(previcureq0)
  results$symp_cases <-  round(previcureq01) 
  results$hospital_surge_beds <- round(previcureq1)
  results$icu_beds <- round(previcureq21)
  results$ventilators <- round(previcureq31)

  results$death_natural_non_exposed <- round(base_mort_S1)
  results$death_natural_exposed <- round(base_mort_E1 + base_mort_I1 + base_mort_CL1 + base_mort_X1 + base_mort_QS1 + 
                                           base_mort_QE1 + base_mort_QI1 + base_mort_QC1 + base_mort_QR1 + base_mort_R1)
  results$death_treated_hospital <- round(cinc_mort_H1)
  results$death_treated_icu <- round(cinc_mort_ICU1)
  results$death_treated_ventilator <- round(cinc_mort_Vent1)
  results$death_untreated_hospital <- round(cinc_mort_HC1)
  results$death_untreated_icu <- round(cinc_mort_ICUC1)
  results$death_untreated_ventilator <- round(cinc_mort_VentC1)
  results$total_deaths <- results$death_treated_hospital + results$death_treated_icu + results$death_treated_ventilator +
    results$death_untreated_hospital + results$death_untreated_icu + results$death_untreated_ventilator
  results$total_deaths_end <- last(results$total_deaths)
  
  
  # !!!! code re-using of variable names but with different str() !!! - request some cleaning
  # START Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  ## AGE DEPENDENT MORTALITY
  cinc_mort_H1 <- parameters["nus"]*parameters["pdeath_h"]*(out[,(Hindex+1)]) +parameters["nus"]*parameters["pdeath_h"]*(out[,(HVindex+1)])
  cinc_mort_HC1 <- parameters["nusc"]*parameters["pdeath_hc"]*(out[,(HCindex+1)])
  cinc_mort_ICU1 <- parameters["nu_icu"]*parameters["pdeath_icu"]*out[,(ICUindex+1)] + parameters["nu_icu"]*parameters["pdeath_icu"]*out[,(ICUVindex+1)]
  cinc_mort_ICUC1 <- parameters["nu_icuc"]*parameters["pdeath_icuc"]*out[,(ICUCindex+1)] 
  cinc_mort_Vent1 <- parameters["nu_vent"]*parameters["pdeath_vent"]*out[,(Ventindex+1)] + parameters["nu_vent"]*parameters["pdeath_vent"]*out[,(VentVindex+1)]
  cinc_mort_VentC1 <- parameters["nu_ventc"]*parameters["pdeath_ventc"]*out[,(VentCindex+1)] 
  totage1<-as.data.frame(cinc_mort_H1+cinc_mort_HC1+cinc_mort_ICU1+cinc_mort_ICUC1+cinc_mort_Vent1+cinc_mort_VentC1)
  basemort_H1<-(out[,(Hindex+1)])
  basemort_HC1<-(out[,(HCindex+1)])
  basemort_ICU1<-(out[,(ICUindex+1)])
  basemort_ICUC1<-(out[,(ICUCindex+1)])
  basemort_Vent1<-(out[,(Ventindex+1)])
  basemort_VentC1<-(out[,(VentCindex+1)])
  basemort_HV<-(out[,(HVindex+1)])
  basemort_ICUV<-(out[,(ICUVindex+1)])
  basemort_VentV<-(out[,(VentVindex+1)])
  
  totbase1<-as.data.frame(basemort_H1+basemort_HC1+basemort_ICU1+basemort_ICUC1+basemort_Vent1+basemort_VentC1+basemort_HV+basemort_ICUV+basemort_VentV)
  tc<-c()
  # END Placeholder for Ricardo/Lisa code (DO NOT EDIT) ----
  for (i in 1:dim(cinc_mort_H1)[1]) {
    for (j in 1:dim(cinc_mort_H1)[2]) {
      tc<-rbind(tc,c(i, j, totage1[i,j]*ifr[j,2]+totbase1[i,j]*mort[j])) 
    }
  }
  tc<-as.data.frame(tc)
  colnames(tc)<-c("Day","Age","value")
  #to add result asymp and symp with age group in tc
  results$tc <- tc %>%
    mutate(Date = startdate + Day,
           age_cat = case_when(
             Age >=  1 & Age <= 6   ~ "≤ 30 y.o.",
             Age >  6 & Age <= 8    ~ "30-40 y.o.",
             Age >  8 & Age <= 10    ~ "40-50 y.o.",
             Age >  10 & Age <= 12    ~ "50-60 y.o.",
             Age >  12 & Age <= 14    ~ "60-70 y.o.",
             Age >=  15  ~ "≥ 70 y.o.")) %>%
    mutate(age_cat = factor(age_cat, levels = rev(c("≤ 30 y.o.", "30-40 y.o.",
                                                    "40-50 y.o.", "50-60 y.o.", "60-70 y.o.", "≥ 70 y.o."))))
  
  mortality_lag <- data.frame(Age = popstruc$agefloor)
  if(nrow(out) >= 30)  mortality_lag <- bind_cols(mortality_lag, 
                                                  data.frame(day30 = out[30,CMindex+1]/out[30,Cindex+1]) %>%
                                                    mutate(day30 = ifelse(is.infinite(day30), 0, day30)) %>%
                                                    rename(`Day 30` = day30))
  if(nrow(out) >= 60)  mortality_lag <- bind_cols(mortality_lag, 
                                                  data.frame(day60 = out[60,CMindex+1]/out[60,Cindex+1]) %>%
                                                    mutate(day60 = ifelse(is.infinite(day60), 0, day60)) %>%
                                                    rename(`Day 60` = day60))
  if(nrow(out) >= 90)  mortality_lag <- bind_cols(mortality_lag, 
                                                  data.frame(day90 = out[90,CMindex+1]/out[90,Cindex+1]) %>%
                                                    mutate(day90 = ifelse(is.infinite(day90), 0, day90)) %>%
                                                    rename(`Day 90` = day90))
  if(nrow(out) >= 120)  mortality_lag <- bind_cols(mortality_lag, 
                                                   data.frame(day120 = out[120,CMindex+1]/out[120,Cindex+1]) %>%
                                                     mutate(day120 = ifelse(is.infinite(day120), 0, day120)) %>%
                                                     rename(`Day 120` = day120))
  
  results$mortality_lag <- mortality_lag
  
  return(results)
}

########### 
###########    RUN BASELINE MODEL - start time for interventions is set to day 1e5, i.e. interventions are always off

Y<-c(initS,initE,initI,initR,initX,initH,initHC,initC,initCM,initSV,initEV,initIV,initCLV,initHV,initICUV,initVentV,initRV, initQS, initQE, initQI,initQR, initCL, initQC, initICU, initICUC, initVent, initVentC, initCMC) # initial conditions for the main solution vector
# out <- ode(y = Y, times = times, func = covid, parms = parameters)
# results <- process_ode_outcome(out)
# plot(x=results$time,y=results$daily_incidence,col="red",type='l',xlab="date",ylab="daily incidence cases",ylim=c(0, 10000))
# lines(incdata_X,type='l')



#### Sompob
# set progress bar
pb <- txtProgressBar(min = 0, max = length(times), style = 3);
parameters["vaccine_eff1"] <- 0.9
parameters["vaccine_cov"] <- 0.5
out.vaceff190 <- ode(y = Y, times = times, func = covid, parms = parameters, method = 'euler',hini=0.05)
results.vaceff90 <- process_ode_outcome(out.vaceff190)

pb <- txtProgressBar(min = 0, max = length(times), style = 3)
parameters["vaccine_eff1"] <- 0.3
parameters["vaccine_cov"] <- 0.5
out.vaceff110 <- ode(y = Y, times = times, func = covid, parms = parameters, method = 'euler',hini=0.05)
results.vaceff10 <- process_ode_outcome(out.vaceff110)



## compare incidence
plot(x=times, y=results.vaceff90$daily_incidence,type = 'l',col='red',ylim = c(0,200000))
lines(x=times,y=results.vaceff10$daily_incidence,type = 'l',col='blue')

## compare incidence
plot(x=times, y=cumsum(results.vaceff90$daily_incidence),type = 'l',col='red')
lines(x=times,y=cumsum(results.vaceff10$daily_incidence),type = 'l',col='blue')


svall.vaceff90 <- out.vaceff190[,SVindex+1]
svall.vaceff10 <- out.vaceff110[,SVindex+1]
rvall.vaceff90 <- out.vaceff190[,RVindex+1]
rvall.vaceff10 <- out.vaceff110[,RVindex+1]

svall.sum.vaceff90 <- rowSums(svall.vaceff90)
svall.sum.vaceff10 <- rowSums(svall.vaceff10)
rvall.sum.vaceff90 <- rowSums(rvall.vaceff90)
rvall.sum.vaceff10 <- rowSums(rvall.vaceff10)

plot(x=times,y=svall.sum.vaceff90,type='l',col='red')
lines(x=times,y=svall.sum.vaceff10,type = 'l',col='blue')

plot(x=times,y=rvall.sum.vaceff90,type='l',col='red',ylim = c(0,10^7))
lines(x=times,y=rvall.sum.vaceff10,type = 'l',col='blue')


sall <- out.vaceff190[,Sindex+1]

sall.sum <- rowSums(sall)

plot(x=times,y=log10(sall.sum),type='l',col='red')


#######



sum(results$asymp_cases)/(11)
sum(results$symp_cases)/11
sum(results$hospital_surge_beds)/11 
sum(results$icu_beds)/11
sum(results$ventilators)/11
results

par(mfrow=c(1,1))
plot(out,select = c(2,3))

plot(out)
out[200:300,SVindex]

plot(times,rowSums(out[,Sindex]))
plot(times,rowSums(out[,SVindex]))

plot(times,rowSums(out[,Eindex]))
plot(times,rowSums(out[,EVindex]))

plot(times,rowSums(out[,Rindex]))
plot(times,rowSums(out[,RVindex]))



colSums()

plot(out,select = c(43,65))
plot(out,select = )
plot(out,select = c(199,211))

plot(out,select = c(199,200))
plot(out,select = c(201,221)) #plot SV and EV

plot(out,select = c(341,351))

sum(rowSums(out[,(Cindex+1)]))/11 #cumulative cases
sum(rowSums(out[,(Iindex+1)]))/11 + sum(rowSums(out[,(IVindex+1)]))/11  #asymptomatic cases
sum(rowSums(out[,(CLindex+1)]))/11 + sum(rowSums(out[,(CLVindex+1)]))/11 #mild symp cases
plot(results$asymp_cases)
sum(results$asymp_cases)/(11)
sum(results$symp_cases)/11
sum(results$hospital_surge_beds)/11 
sum(results$icu_beds)/11
sum(results$ventilators)/11
results








 Y<-c(initS,initE,initI,initR,initX,initH,initHC,initC,initCM,initSV,initEV,initIV,initCLV,initHV,initICUV,initVentV,initRV, initQS, initQE, initQI,initQR, initCL, initQC, initICU, initICUC, initVent, initVentC, initCMC) # initial conditions for the main solution vector
out1 <- ode(y = Y, times = times, method = "euler", hini = 0.05, func = covid, parms = parameters)

Y<-c(initS,initE,initI,initR,initX,initH,initHC,initC,initCM,initSV,initEV,initIV,initCLV,initHV,initICUV,initVentV,initRV, initQS, initQE, initQI,initQR, initCL, initQC, initICU, initICUC, initVent, initVentC, initCMC) # initial conditions for the main solution vector
out2 <- ode(y = Y, times = times, method = "euler", hini = 0.05, func = covid, parms = parameters)

#out0 <- ode(y = Y, times = times, func = covid, parms = parameters)


...
#for vaccination; add date and vaccine intervention par then re run with out1/2/3 and parameters1/2/3, as well as plot with line

plot(x=results$time,y=results$asymp_cases,col="red",type='l',xlab="date",ylab="daily asymp cases")
lines(incdata_X,type='l')





results <- process_ode_outcome(out)
results1 <- process_ode_outcome(out1)
results1
plot(x=results1$time,y=results1$daily_incidence,col="red",type='l',xlab="date",ylab="daily incidence cases",ylim=c(0, 40000))
lines(incdata_X,type='l')


results2 <- process_ode_outcome(out1)
results2
plot(x=results2$time,y=results2$daily_incidence,col="red",type='l',xlab="date",ylab="daily incidence cases",ylim=c(0, 15000))
lines(incdata_X,type='l')



##### list of variables that can be visualised
#names(results)
results$tc
results1$tc
results2$tc

results$pct_total_pop_infected
results1$pct_total_pop_infected
results2$pct_total_pop_infected

results$mortality_lag
results1$mortality_lag
results2$mortality_lag

sum(results$daily_total_cases)
sum(results1$daily_total_cases)

results$total_deaths_end/sum(results$daily_total_cases)

results$total_deaths_end
results1$total_deaths_end
results2$total_deaths_end


### daily incidence cases
plot(x=results$time,y=results$daily_incidence,col="red",type='l',xlab="date",ylab="daily incidence cases")
lines(incdata_X,type='l')


plot(x=results1$time,y=results$daily_incidence,col="red",type='l',xlab="date",ylab="daily incidence cases")
lines(incdata_X,type='l')

### R(t)
plot(x=results$time, y=results$Rt, type='l',xlab = "date",ylab = "predicted R(t)")




######################
par1<-list(vaccine_eff1,output_model_vac1_age,output_model_vac4_age,output_model_vac6_age,output_model_vac7_age,output_model_vac8_age)

par2<-list(output_model_vac2_age,output_model_vac1_age,output_model_vac4_age,output_model_vac6_age,output_model_vac7_age,output_model_vac8_age)

Vac_TPP <- function (par1, par2, par3, par4)

#Step2:Create function
##Function###
# to evaluate Rotasiil, do it quickly just change par2 to be cost of ROTASIIL instead of ROTAVAC
#line416 from actual cost to par2 and line415 (Rotavac) to actual cost 
# also change output to get results from SIIL line 1073-74 -> CQICER<-c(mean_Incost_vac4_base,mean_InQALY_vac4_base,mean_ICER_vac4_base)
###total case by age and year IPD
par2<-list(output_model_vac2_age,output_model_vac1_age,output_model_vac4_age,output_model_vac6_age,output_model_vac7_age,output_model_vac8_age)

TWSA <- function(par1, par2){
  Total_case_IPD <- array(NA, dim=c(4,5,1000)) 
  ###total case by age and year OPD
  Total_case_OPD <- array(NA, dim=c(4,5,1000)) 
  ###total case by age and year self
  Total_case_self <- array(NA, dim=c(4,5,1000)) 
  
  for (k in 1:1000) 
  {
    Total_case_IPD[1,1:5,k]<-prob_IPD1[k]*output_model_age[1,1:5,k]
    Total_case_IPD[2,1:5,k]<-prob_IPD2[k]*output_model_age[2,1:5,k]
    Total_case_IPD[3,1:5,k]<-prob_IPD3[k]*output_model_age[3,1:5,k]
    Total_case_IPD[4,1:5,k]<-prob_IPD4[k]*output_model_age[4,1:5,k]
    
    Total_case_OPD[1,1:5,k]<-prob_OPD1[k]*output_model_age[1,1:5,k]
    Total_case_OPD[2,1:5,k]<-prob_OPD2[k]*output_model_age[2,1:5,k]
    Total_case_OPD[3,1:5,k]<-prob_OPD3[k]*output_model_age[3,1:5,k]
    Total_case_OPD[4,1:5,k]<-prob_OPD4[k]*output_model_age[4,1:5,k]
    
    Total_case_self[1,1:5,k]<-prob_self1[k]*output_model_age[1,1:5,k]
    Total_case_self[2,1:5,k]<-prob_self2[k]*output_model_age[2,1:5,k]
    Total_case_self[3,1:5,k]<-prob_self3[k]*output_model_age[3,1:5,k]
    Total_case_self[4,1:5,k]<-prob_self4[k]*output_model_age[4,1:5,k]
  }
  Total_case<-Total_case_self+Total_case_OPD+Total_case_IPD
  Total_case
  
  ###case with vaccine 1
  ###total case by age and year IPD
  Total_case_IPD_vac1 <- array(NA, dim=c(4,5,1000)) 
  ###total case by age and year OPD
  Total_case_OPD_vac1 <- array(NA, dim=c(4,5,1000)) 
  ###total case by age and year self
  Total_case_self_vac1 <- array(NA, dim=c(4,5,1000)) 
  
  for (k in 1:1000) 
  {
    Total_case_IPD_vac1[1,1:5,k]<-prob_IPD1[k]*par2[1,1:5,k]
    Total_case_IPD_vac1[2,1:5,k]<-prob_IPD2[k]*par2[2,1:5,k]
    Total_case_IPD_vac1[3,1:5,k]<-prob_IPD3[k]*par2[3,1:5,k]
    Total_case_IPD_vac1[4,1:5,k]<-prob_IPD4[k]*par2[4,1:5,k]
    
    Total_case_OPD_vac1[1,1:5,k]<-prob_OPD1[k]*par2[1,1:5,k]
    Total_case_OPD_vac1[2,1:5,k]<-prob_OPD2[k]*par2[2,1:5,k]
    Total_case_OPD_vac1[3,1:5,k]<-prob_OPD3[k]*par2[3,1:5,k]
    Total_case_OPD_vac1[4,1:5,k]<-prob_OPD4[k]*par2[4,1:5,k]
    
    Total_case_self_vac1[1,1:5,k]<-prob_self1[k]*par2[1,1:5,k]
    Total_case_self_vac1[2,1:5,k]<-prob_self2[k]*par2[2,1:5,k]
    Total_case_self_vac1[3,1:5,k]<-prob_self3[k]*par2[3,1:5,k]
    Total_case_self_vac1[4,1:5,k]<-prob_self4[k]*par2[4,1:5,k]
  }
  Total_case_vac1<-Total_case_self_vac1+Total_case_OPD_vac1+Total_case_IPD_vac1
  Total_case_vac1
  ########################
  #EPI part finished
  
  #PAR4: Risk of each condition (Mild and Severe)
  #rc_all<-array(data=NA, c(no.of.dtp,4,no.of.agegr)) #number is age group
  
  #COST ESTIMATION: Cost of Vaccination and Cost of treatment due to RV (Direct medical cost, Direct non medical cost, and Indirect cost (Gov=1+2, Soc=1+2+3) by 4 age groups
  ##Cost of Vaccination
  wastage1 <- 1
  wastage2 <- 1
  wastage3 <- 30
  wastage4 <- 5
  Rotarix <-	274.98/(1-wastage1/100)   #2 doses
  RotaTeq <-	183.32/(1-wastage2/100)   #3 doses
  RotaVAC <-	31.31/(1-wastage3/100) #3 doses 31.31
  RotaSIIL <-	par1/(1-wastage4/100) #3 doses 31.31
  
  
  #Cost_logis <-	4.85 #additioal cost per dose for logistics [Cost of logistics (per additional dose)	4.85	1.73	6.59] This include expandsion of the fridge
  #Cost_adv <-	11.94 #Cost_adverse_event form onwipa (adverse=diarrhea) may assume +/-10%
  #Cost_intus <- 61.06 #per vaccinated population or 6,106 baht per intusseception case
  
  
  ################# 
  #Vaccination Cost
  #################
  # cost_vac1 <- Rotarix + Cost_logistics +Cost_adverse_event
  # cost_vac2 <- RotaTeq + Cost_logistics +Cost_adverse_event
  # cost_vac3 <- RotaVAC + Cost_logistics +Cost_adverse_event
  # cost_vac4 <- RotaSIIL + Cost_logistics +Cost_adverse_event
  
  Cost_intus<-rep(0, 1000) # for base case
  
  #cost vac1 among population age 0-2 m, 3-4m, 5-6m in 2020
  pop_age1_y1 <- 700000
  pop_age2_y1 <- 700000
  pop_age3_y1 <- 700000
  #cost vac1 among population age 0-2 m, 3-4m, 5-6m in 2021
  pop_age1_y2 <- pop_age1_y1-2000
  pop_age2_y2 <- pop_age1_y1-2000
  pop_age3_y2 <- pop_age1_y1-2000
  #cost vac1 among population age 0-2 m, 3-4m, 5-6m in 2022
  pop_age1_y3 <- pop_age1_y2-2000
  pop_age2_y3 <- pop_age1_y2-2000
  pop_age3_y3 <- pop_age1_y2-2000
  #cost vac1 among population age 0-2 m, 3-4m, 5-6m in 2023
  pop_age1_y4 <- pop_age1_y3-2000
  pop_age2_y4 <- pop_age1_y3-2000
  pop_age3_y4 <- pop_age1_y3-2000
  #cost vac1 among population age 0-2 m, 3-4m, 5-6m in 2024
  pop_age1_y5 <- pop_age1_y4-2000
  pop_age2_y5 <- pop_age1_y4-2000
  pop_age3_y5 <- pop_age1_y4-2000
  
  #Vaccince Costs
  #Vac1 (2 doses)
  #Add loop from here to 593 to add PSA effect (1,000 run)
  #use it as vector
  Cost_per_dose_vac1d1_y1 <- rep(NA, 1000) #vector 
  Cost_per_dose_vac1d2_y1 <- rep(NA, 1000)
  Cost_per_dose_vac1d1_y2 <- rep(NA, 1000)
  Cost_per_dose_vac1d2_y2 <- rep(NA, 1000)
  Cost_per_dose_vac1d1_y3 <- rep(NA, 1000)
  Cost_per_dose_vac1d2_y3 <- rep(NA, 1000)
  Cost_per_dose_vac1d1_y4 <- rep(NA, 1000)
  Cost_per_dose_vac1d2_y4 <- rep(NA, 1000)
  Cost_per_dose_vac1d1_y5 <- rep(NA, 1000)
  Cost_per_dose_vac1d2_y5 <- rep(NA, 1000)
  
  
  Total_cost_vac1_y1<-rep(NA, 1000)
  Total_cost_vac1_y2<-rep(NA, 1000)
  Total_cost_vac1_y3<-rep(NA, 1000)
  Total_cost_vac1_y4<-rep(NA, 1000)
  Total_cost_vac1_y5<-rep(NA, 1000)
  #Vac2
  Cost_per_dose_vac2d1_y1<-rep(NA, 1000)
  Cost_per_dose_vac2d2_y1<-rep(NA, 1000)
  Cost_per_dose_vac2d3_y1<-rep(NA, 1000)
  Cost_per_dose_vac2d1_y2<-rep(NA, 1000)
  Cost_per_dose_vac2d2_y2<-rep(NA, 1000)
  Cost_per_dose_vac2d3_y2<-rep(NA, 1000)
  Cost_per_dose_vac2d1_y3<-rep(NA, 1000)
  Cost_per_dose_vac2d2_y3<-rep(NA, 1000)
  Cost_per_dose_vac2d3_y3<-rep(NA, 1000)
  Cost_per_dose_vac2d1_y4<-rep(NA, 1000)
  Cost_per_dose_vac2d2_y4<-rep(NA, 1000)
  Cost_per_dose_vac2d3_y4<-rep(NA, 1000)
  Cost_per_dose_vac2d1_y5<-rep(NA, 1000)
  Cost_per_dose_vac2d2_y5<-rep(NA, 1000)
  Cost_per_dose_vac2d3_y5<-rep(NA, 1000)
  
  Total_cost_vac2_y1<-rep(NA, 1000)
  Total_cost_vac2_y2<-rep(NA, 1000)
  Total_cost_vac2_y3<-rep(NA, 1000)
  Total_cost_vac2_y4<-rep(NA, 1000)
  Total_cost_vac2_y5<-rep(NA, 1000)
  #Vac3
  Cost_per_dose_vac3d1_y1<-rep(NA, 1000)
  Cost_per_dose_vac3d2_y1<-rep(NA, 1000)
  Cost_per_dose_vac3d3_y1<-rep(NA, 1000)
  Cost_per_dose_vac3d1_y2<-rep(NA, 1000)
  Cost_per_dose_vac3d2_y2<-rep(NA, 1000)
  Cost_per_dose_vac3d3_y2<-rep(NA, 1000)
  Cost_per_dose_vac3d1_y3<-rep(NA, 1000)
  Cost_per_dose_vac3d2_y3<-rep(NA, 1000)
  Cost_per_dose_vac3d3_y3<-rep(NA, 1000)
  Cost_per_dose_vac3d1_y4<-rep(NA, 1000)
  Cost_per_dose_vac3d2_y4<-rep(NA, 1000)
  Cost_per_dose_vac3d3_y4<-rep(NA, 1000)
  Cost_per_dose_vac3d1_y5<-rep(NA, 1000)
  Cost_per_dose_vac3d2_y5<-rep(NA, 1000)
  Cost_per_dose_vac3d3_y5<-rep(NA, 1000)
  
  Total_cost_vac3_y1<-rep(NA, 1000)
  Total_cost_vac3_y2<-rep(NA, 1000)
  Total_cost_vac3_y3<-rep(NA, 1000)
  Total_cost_vac3_y4<-rep(NA, 1000)
  Total_cost_vac3_y5<-rep(NA, 1000)
  
  #Vac4
  Cost_per_dose_vac4d1_y1<-rep(NA, 1000)
  Cost_per_dose_vac4d2_y1<-rep(NA, 1000)
  Cost_per_dose_vac4d3_y1<-rep(NA, 1000)
  Cost_per_dose_vac4d1_y2<-rep(NA, 1000)
  Cost_per_dose_vac4d2_y2<-rep(NA, 1000)
  Cost_per_dose_vac4d3_y2<-rep(NA, 1000)
  Cost_per_dose_vac4d1_y3<-rep(NA, 1000)
  Cost_per_dose_vac4d2_y3<-rep(NA, 1000)
  Cost_per_dose_vac4d3_y3<-rep(NA, 1000)
  Cost_per_dose_vac4d1_y4<-rep(NA, 1000)
  Cost_per_dose_vac4d2_y4<-rep(NA, 1000)
  Cost_per_dose_vac4d3_y4<-rep(NA, 1000)
  Cost_per_dose_vac4d1_y5<-rep(NA, 1000)
  Cost_per_dose_vac4d2_y5<-rep(NA, 1000)
  Cost_per_dose_vac4d3_y5<-rep(NA, 1000)
  
  Total_cost_vac4_y1<-rep(NA, 1000)
  Total_cost_vac4_y2<-rep(NA, 1000)
  Total_cost_vac4_y3<-rep(NA, 1000)
  Total_cost_vac4_y4<-rep(NA, 1000)
  Total_cost_vac4_y5<-rep(NA, 1000)
  #then, vac2 and 3 year 1-5 and wrap with for loop
  
  for (k in 1:1000) 
  {
    #Rotarix 2 doses  
    Cost_per_dose_vac1d1_y1[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age1_y1*cov1)
    Cost_per_dose_vac1d2_y1[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age2_y1*cov2)
    Cost_per_dose_vac1d1_y2[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age1_y2*cov1)
    Cost_per_dose_vac1d2_y2[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age2_y2*cov2)
    Cost_per_dose_vac1d1_y3[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age1_y3*cov1)
    Cost_per_dose_vac1d2_y3[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age2_y3*cov2)
    Cost_per_dose_vac1d1_y4[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age1_y4*cov1)
    Cost_per_dose_vac1d2_y4[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age2_y4*cov2)
    Cost_per_dose_vac1d1_y5[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age1_y5*cov1)
    Cost_per_dose_vac1d2_y5[k] <- (Rotarix+Cost_logis[k]+Cost_adv[k])*(pop_age2_y5*cov2)
    
    Total_cost_vac1_y1[k]<- Cost_per_dose_vac1d1_y1[k] +Cost_per_dose_vac1d2_y1[k]
    Total_cost_vac1_y2[k]<- Cost_per_dose_vac1d1_y2[k] +Cost_per_dose_vac1d2_y2[k] 
    Total_cost_vac1_y3[k]<- Cost_per_dose_vac1d1_y3[k] +Cost_per_dose_vac1d2_y3[k]
    Total_cost_vac1_y4[k]<- Cost_per_dose_vac1d1_y4[k] +Cost_per_dose_vac1d2_y4[k] 
    Total_cost_vac1_y5[k]<- Cost_per_dose_vac1d1_y5[k] +Cost_per_dose_vac1d2_y5[k] 
  }
  
  for (k in 1:1000) 
  {
    #Rotateq 3 doses
    Cost_per_dose_vac2d1_y1[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age1_y1*cov1)
    Cost_per_dose_vac2d2_y1[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age2_y1*cov2)
    Cost_per_dose_vac2d3_y1[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age3_y1*cov3)
    
    Cost_per_dose_vac2d1_y2[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age1_y2*cov1)
    Cost_per_dose_vac2d2_y2[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age2_y2*cov2)
    Cost_per_dose_vac2d3_y2[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age3_y2*cov3)
    
    Cost_per_dose_vac2d1_y3[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age1_y3*cov1)
    Cost_per_dose_vac2d2_y3[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age2_y3*cov2)
    Cost_per_dose_vac2d3_y3[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age3_y3*cov3)
    
    Cost_per_dose_vac2d1_y4[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age1_y4*cov1)
    Cost_per_dose_vac2d2_y4[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age2_y4*cov2)
    Cost_per_dose_vac2d3_y4[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age3_y4*cov3)
    
    Cost_per_dose_vac2d1_y5[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age1_y5*cov1)
    Cost_per_dose_vac2d2_y5[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age2_y5*cov2)
    Cost_per_dose_vac2d3_y5[k] <- (RotaTeq+Cost_logis[k]+Cost_adv[k])*(pop_age3_y5*cov3)
    
    Total_cost_vac2_y1[k]<- Cost_per_dose_vac2d1_y1[k]+Cost_per_dose_vac2d2_y1[k]+Cost_per_dose_vac2d3_y1[k]
    Total_cost_vac2_y2[k]<- Cost_per_dose_vac2d1_y2[k]+Cost_per_dose_vac2d2_y2[k]+Cost_per_dose_vac2d3_y2[k]
    Total_cost_vac2_y3[k]<- Cost_per_dose_vac2d1_y3[k]+Cost_per_dose_vac2d2_y3[k]+Cost_per_dose_vac2d3_y3[k]
    Total_cost_vac2_y4[k]<- Cost_per_dose_vac2d1_y4[k]+Cost_per_dose_vac2d2_y4[k]+Cost_per_dose_vac2d3_y4[k]
    Total_cost_vac2_y5[k]<- Cost_per_dose_vac2d1_y5[k]+Cost_per_dose_vac2d2_y5[k]+Cost_per_dose_vac2d3_y5[k]
  }
  for (k in 1:1000) 
  {
    #cost vac3 among population age 0-2 m, 3-4m, 5-6m in 2020
    #Vac3
    Cost_per_dose_vac3d1_y1[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y1*cov1)
    Cost_per_dose_vac3d2_y1[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y1*cov2)
    Cost_per_dose_vac3d3_y1[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y1*cov3)
    Cost_per_dose_vac3d1_y2[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y2*cov1)
    Cost_per_dose_vac3d2_y2[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y2*cov2)
    Cost_per_dose_vac3d3_y2[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y2*cov3)
    Cost_per_dose_vac3d1_y3[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y3*cov1)
    Cost_per_dose_vac3d2_y3[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y3*cov2)
    Cost_per_dose_vac3d3_y3[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y3*cov3)
    Cost_per_dose_vac3d1_y4[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y4*cov1)
    Cost_per_dose_vac3d2_y4[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y4*cov2)
    Cost_per_dose_vac3d3_y4[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y4*cov3)
    Cost_per_dose_vac3d1_y5[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y5*cov1)
    Cost_per_dose_vac3d2_y5[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y5*cov2)
    Cost_per_dose_vac3d3_y5[k] <- (RotaVAC+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y5*cov3)
    
    
    Total_cost_vac3_y1[k]<- Cost_per_dose_vac3d1_y1[k] +Cost_per_dose_vac3d2_y1[k] +Cost_per_dose_vac3d3_y1[k]
    Total_cost_vac3_y2[k]<- Cost_per_dose_vac3d1_y2[k] +Cost_per_dose_vac3d2_y2[k] +Cost_per_dose_vac3d3_y2[k]
    Total_cost_vac3_y3[k]<- Cost_per_dose_vac3d1_y3[k] +Cost_per_dose_vac3d2_y3[k] +Cost_per_dose_vac3d3_y3[k]
    Total_cost_vac3_y4[k]<- Cost_per_dose_vac3d1_y4[k] +Cost_per_dose_vac3d2_y4[k] +Cost_per_dose_vac3d3_y4[k]
    Total_cost_vac3_y5[k]<- Cost_per_dose_vac3d1_y5[k] +Cost_per_dose_vac3d2_y5[k] +Cost_per_dose_vac3d3_y5[k]
    
    #cost vac4 among population age 0-2 m, 3-4m, 5-6m in 2020
    #Vac4
    Cost_per_dose_vac4d1_y1[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y1*cov1)
    Cost_per_dose_vac4d2_y1[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y1*cov2)
    Cost_per_dose_vac4d3_y1[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y1*cov3)
    Cost_per_dose_vac4d1_y2[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y2*cov1)
    Cost_per_dose_vac4d2_y2[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y2*cov2)
    Cost_per_dose_vac4d3_y2[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y2*cov3)
    Cost_per_dose_vac4d1_y3[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y3*cov1)
    Cost_per_dose_vac4d2_y3[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y3*cov2)
    Cost_per_dose_vac4d3_y3[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y3*cov3)
    Cost_per_dose_vac4d1_y4[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y4*cov1)
    Cost_per_dose_vac4d2_y4[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y4*cov2)
    Cost_per_dose_vac4d3_y4[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y4*cov3)
    Cost_per_dose_vac4d1_y5[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age1_y5*cov1)
    Cost_per_dose_vac4d2_y5[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age2_y5*cov2)
    Cost_per_dose_vac4d3_y5[k] <- (RotaSIIL+Cost_logis[k]+Cost_adv[k]+Cost_intus[k])*(pop_age3_y5*cov3)
    
    Total_cost_vac4_y1[k]<-Cost_per_dose_vac4d1_y1[k]+Cost_per_dose_vac4d2_y1[k]+Cost_per_dose_vac4d3_y1[k]
    Total_cost_vac4_y2[k]<-Cost_per_dose_vac4d1_y2[k]+Cost_per_dose_vac4d2_y2[k]+Cost_per_dose_vac4d3_y2[k]
    Total_cost_vac4_y3[k]<-Cost_per_dose_vac4d1_y3[k]+Cost_per_dose_vac4d2_y3[k]+Cost_per_dose_vac4d3_y3[k]
    Total_cost_vac4_y4[k]<-Cost_per_dose_vac4d1_y4[k]+Cost_per_dose_vac4d2_y4[k]+Cost_per_dose_vac4d3_y4[k]
    Total_cost_vac4_y5[k]<-Cost_per_dose_vac4d1_y5[k]+Cost_per_dose_vac4d2_y5[k]+Cost_per_dose_vac4d3_y5[k]
    
  }
  #problem with Total_cost_vac2-4_yx
  #test
  Cost_per_dose_vac1d1_y4
  Cost_per_dose_vac2d3_y4
  Cost_per_dose_vac3d1_y4
  Cost_per_dose_vac4d2_y4
  Total_cost_vac1_y3
  Total_cost_vac1_y2-Total_cost_vac2_y2
  Total_cost_vac1_y3-Total_cost_vac3_y3
  Total_cost_vac3_y3-Total_cost_vac2_y3
  Total_cost_vac4_y1-Total_cost_vac2_y1
  Total_cost_vac3_y2-Total_cost_vac4_y2
  Total_cost_vac3_y2-Total_cost_vac3_y1
  mean(Total_cost_vac1_y1)
  mean(Total_cost_vac2_y1)
  mean(Total_cost_vac3_y1)
  mean(Total_cost_vac1_y2)
  mean(Total_cost_vac2_y2)
  mean(Total_cost_vac3_y2)
  mean(Total_cost_vac1_y3)
  mean(Total_cost_vac2_y3)
  mean(Total_cost_vac3_y3)
  mean(Total_cost_vac1_y4)
  mean(Total_cost_vac2_y4)
  mean(Total_cost_vac3_y4)
  
  ###################################################
  #Treatment Cost and Utility by each condition by age group (4 groups) [Loop with 1000 iterations]
  ###################################################
  ###Direct medical cost, Direct non medical cost, and Indirect cost (Gov=1, Soc=1+2+3) by 4 age groups
  
  
  ## without vaccine aged 0-5 in 2019
  #societal cost 1(direact med) cost 2 (direct non-med) cost 3 (indirect)
  #prob_IPD1 (risk IPD) prob_OPD1 (risk OPD) prob_self1 (risk self)
  ##cost without vaccince
  #Total cost outcome
  Total_Cost_soc <- array(NA, dim=c(4,5,1000)) 
  Total_Cost_gov <- array(NA, dim=c(4,5,1000)) 
  Total_Cost_vac1_soc <- array(NA, dim=c(4,5,1000)) 
  Total_Cost_vac1_gov <- array(NA, dim=c(4,5,1000)) 
  
  for (k in 1:1000) 
  {
    Total_Cost_soc[1,1:5,k]<-((cost_IPD1[k]*prob_IPD1[k])+(cost_IPD2[k]*prob_IPD1[k])+(cost_IPD3[k]*prob_IPD1[k])+(cost_OPD1[k]*prob_OPD1[k])+(cost_OPD2[k]*prob_OPD1[k])+(cost_OPD3[k]*prob_OPD1[k])+(cost_self1[k]*prob_self1[k])+(cost_self2[k]*prob_self1[k])+(cost_self3[k]*prob_self1[k]))*output_model_age[1,1:5,k]
    Total_Cost_soc[2,1:5,k]<-((cost_IPD4[k]*prob_IPD2[k])+(cost_IPD5[k]*prob_IPD2[k])+(cost_IPD6[k]*prob_IPD2[k])+(cost_OPD4[k]*prob_OPD2[k])+(cost_OPD5[k]*prob_OPD2[k])+(cost_OPD6[k]*prob_OPD2[k])+(cost_self4[k]*prob_self2[k])+(cost_self5[k]*prob_self2[k])+(cost_self6[k]*prob_self2[k]))*output_model_age[2,1:5,k]
    Total_Cost_soc[3,1:5,k]<-((cost_IPD7[k]*prob_IPD3[k])+(cost_IPD8[k]*prob_IPD3[k])+(cost_IPD9[k]*prob_IPD3[k])+(cost_OPD7[k]*prob_OPD3[k])+(cost_OPD8[k]*prob_OPD3[k])+(cost_OPD9[k]*prob_OPD3[k])+(cost_self7[k]*prob_self3[k])+(cost_self8[k]*prob_self3[k])+(cost_self9[k]*prob_self3[k]))*output_model_age[3,1:5,k]
    Total_Cost_soc[4,1:5,k]<-((cost_IPD10[k]*prob_IPD4[k])+(cost_IPD11[k]*prob_IPD4[k])+(cost_IPD12[k]*prob_IPD4[k])+(cost_OPD10[k]*prob_OPD4[k])+(cost_OPD11[k]*prob_OPD4[k])+(cost_OPD12[k]*prob_OPD4[k])+(cost_self10[k]*prob_self4[k])+(cost_self11[k]*prob_self4[k])+(cost_self12[k]*prob_self4[k]))*output_model_age[4,1:5,k]
    
    Total_Cost_gov[1,1:5,k]<-((cost_IPD1[k]*prob_IPD1[k])+(cost_IPD2[k]*prob_IPD1[k])+(cost_IPD3[k]*prob_IPD1[k]))*output_model_age[1,1:5,k]
    Total_Cost_gov[2,1:5,k]<-((cost_IPD4[k]*prob_IPD2[k])+(cost_IPD5[k]*prob_IPD2[k])+(cost_IPD6[k]*prob_IPD2[k]))*output_model_age[2,1:5,k]
    Total_Cost_gov[3,1:5,k]<-((cost_IPD7[k]*prob_IPD3[k])+(cost_IPD8[k]*prob_IPD3[k])+(cost_IPD9[k]*prob_IPD3[k]))*output_model_age[3,1:5,k]
    Total_Cost_gov[4,1:5,k]<-((cost_IPD10[k]*prob_IPD4[k])+(cost_IPD11[k]*prob_IPD4[k])+(cost_IPD12[k]*prob_IPD4[k]))*output_model_age[4,1:5,k]
    
    Total_Cost_vac1_soc[1,1:5,k]<-((cost_IPD1[k]*prob_IPD1[k])+(cost_IPD2[k]*prob_IPD1[k])+(cost_IPD3[k]*prob_IPD1[k])+(cost_OPD1[k]*prob_OPD1[k])+(cost_OPD2[k]*prob_OPD1[k])+(cost_OPD3[k]*prob_OPD1[k])+(cost_self1[k]*prob_self1[k])+(cost_self2[k]*prob_self1[k])+(cost_self3[k]*prob_self1[k]))*par2[1,1:5,k]
    Total_Cost_vac1_soc[2,1:5,k]<-((cost_IPD4[k]*prob_IPD2[k])+(cost_IPD5[k]*prob_IPD2[k])+(cost_IPD6[k]*prob_IPD2[k])+(cost_OPD4[k]*prob_OPD2[k])+(cost_OPD5[k]*prob_OPD2[k])+(cost_OPD6[k]*prob_OPD2[k])+(cost_self4[k]*prob_self2[k])+(cost_self5[k]*prob_self2[k])+(cost_self6[k]*prob_self2[k]))*par2[2,1:5,k]
    Total_Cost_vac1_soc[3,1:5,k]<-((cost_IPD7[k]*prob_IPD3[k])+(cost_IPD8[k]*prob_IPD3[k])+(cost_IPD9[k]*prob_IPD3[k])+(cost_OPD7[k]*prob_OPD3[k])+(cost_OPD8[k]*prob_OPD3[k])+(cost_OPD9[k]*prob_OPD3[k])+(cost_self7[k]*prob_self3[k])+(cost_self8[k]*prob_self3[k])+(cost_self9[k]*prob_self3[k]))*par2[3,1:5,k]
    Total_Cost_vac1_soc[4,1:5,k]<-((cost_IPD10[k]*prob_IPD4[k])+(cost_IPD11[k]*prob_IPD4[k])+(cost_IPD12[k]*prob_IPD4[k])+(cost_OPD10[k]*prob_OPD4[k])+(cost_OPD11[k]*prob_OPD4[k])+(cost_OPD12[k]*prob_OPD4[k])+(cost_self10[k]*prob_self4[k])+(cost_self11[k]*prob_self4[k])+(cost_self12[k]*prob_self4[k]))*par2[4,1:5,k]
    
    Total_Cost_vac1_gov[1,1:5,k]<-((cost_IPD1[k]*prob_IPD1[k])+(cost_IPD2[k]*prob_IPD1[k])+(cost_IPD3[k]*prob_IPD1[k]))*par2[1,1:5,k]
    Total_Cost_vac1_gov[2,1:5,k]<-((cost_IPD4[k]*prob_IPD2[k])+(cost_IPD5[k]*prob_IPD2[k])+(cost_IPD6[k]*prob_IPD2[k]))*par2[2,1:5,k]
    Total_Cost_vac1_gov[3,1:5,k]<-((cost_IPD7[k]*prob_IPD3[k])+(cost_IPD8[k]*prob_IPD3[k])+(cost_IPD9[k]*prob_IPD3[k]))*par2[3,1:5,k]
    Total_Cost_vac1_gov[4,1:5,k]<-((cost_IPD10[k]*prob_IPD4[k])+(cost_IPD11[k]*prob_IPD4[k])+(cost_IPD12[k]*prob_IPD4[k]))*par2[4,1:5,k]
  }
  
  Total_Cost_gov
  Total_Cost_soc
  
  Total_Cost_vac1_soc
  Total_Cost_vac1_gov
  
  ####################################################
  # LE and (QA)LY (by age groups; 0-5 and >5)
  ####################################################
  # estimate the QALY loss per rotavirus infection case by age group using LE (WHO or others) (early death) and Qol during infection(Utility loss*Length of infection) 
  #Total QALY loss outcome without vac and with vaccine
  #Total treatment cost outcome vac1
  
  Total_QALY_loss_novac <- array(NA, dim=c(4,5,1000)) 
  Total_QALY_loss_vac1 <- array(NA, dim=c(4,5,1000)) 
  
  
  for (k in 1:1000) 
  {
    Total_QALY_loss_novac[1,1:5,k]<-((IPD1[k]/365*(1-uti_IPD1[k]))+(OPD[k]/365*(1-uti_OPD1[k]))+(self[k]/365*(1-uti_self1[k])))*output_model_age[1,1:5,k]
    Total_QALY_loss_novac[2,1:5,k]<-((IPD2[k]/365*(1-uti_IPD2[k]))+(OPD[k]/365*(1-uti_OPD2[k]))+(self[k]/365*(1-uti_self2[k])))*output_model_age[2,1:5,k]
    Total_QALY_loss_novac[3,1:5,k]<-((IPD3[k]/365*(1-uti_IPD2[k]))+(OPD[k]/365*(1-uti_OPD2[k]))+(self[k]/365*(1-uti_self2[k])))*output_model_age[3,1:5,k]
    Total_QALY_loss_novac[4,1:5,k]<-((IPD4[k]/365*(1-uti_IPD2[k]))+(OPD[k]/365*(1-uti_OPD2[k]))+(self[k]/365*(1-uti_self2[k])))*output_model_age[4,1:5,k]
    
    Total_QALY_loss_vac1[1,1:5,k]<-((IPD1[k]/365*(1-uti_IPD1[k]))+(OPD[k]/365*(1-uti_OPD1[k]))+(self[k]/365*(1-uti_self1[k])))*par2[1,1:5,k]
    Total_QALY_loss_vac1[2,1:5,k]<-((IPD2[k]/365*(1-uti_IPD2[k]))+(OPD[k]/365*(1-uti_OPD2[k]))+(self[k]/365*(1-uti_self2[k])))*par2[2,1:5,k]
    Total_QALY_loss_vac1[3,1:5,k]<-((IPD3[k]/365*(1-uti_IPD2[k]))+(OPD[k]/365*(1-uti_OPD2[k]))+(self[k]/365*(1-uti_self2[k])))*par2[3,1:5,k]
    Total_QALY_loss_vac1[4,1:5,k]<-((IPD4[k]/365*(1-uti_IPD2[k]))+(OPD[k]/365*(1-uti_OPD2[k]))+(self[k]/365*(1-uti_self2[k])))*par2[4,1:5,k]
    
    
  }
  #Parameter lists
  #DISCOUNTING RATE
  disc.rate<- 0.03 #disocunting rate0.03 
  disc.2020<- 1
  disc.2021<- disc.2020-disc.rate
  disc.2022<- disc.2021-disc.rate
  disc.2023<- disc.2022-disc.rate
  disc.2024<- disc.2023-disc.rate
  
  #Vaccine Costs 
  #add loop 
  ##discount vac 1 
  Total_cost_vac1_y1_dis<-Total_cost_vac1_y1*disc.2020
  Total_cost_vac1_y2_dis<-Total_cost_vac1_y2*disc.2021
  Total_cost_vac1_y3_dis<-Total_cost_vac1_y3*disc.2022
  Total_cost_vac1_y4_dis<-Total_cost_vac1_y4*disc.2023
  Total_cost_vac1_y5_dis<-Total_cost_vac1_y5*disc.2024
  
  ##discount vac 2
  Total_cost_vac2_y1_dis<-Total_cost_vac2_y1*disc.2020
  Total_cost_vac2_y2_dis<-Total_cost_vac2_y2*disc.2021
  Total_cost_vac2_y3_dis<-Total_cost_vac2_y3*disc.2022
  Total_cost_vac2_y4_dis<-Total_cost_vac2_y4*disc.2023
  Total_cost_vac2_y5_dis<-Total_cost_vac2_y5*disc.2024
  
  ##discount vac 3
  Total_cost_vac3_y1_dis<-Total_cost_vac3_y1*disc.2020
  Total_cost_vac3_y2_dis<-Total_cost_vac3_y2*disc.2021
  Total_cost_vac3_y3_dis<-Total_cost_vac3_y3*disc.2022
  Total_cost_vac3_y4_dis<-Total_cost_vac3_y4*disc.2023
  Total_cost_vac3_y5_dis<-Total_cost_vac3_y5*disc.2024
  
  ##discount vac 4
  Total_cost_vac4_y1_dis<-Total_cost_vac4_y1*disc.2020
  Total_cost_vac4_y2_dis<-Total_cost_vac4_y2*disc.2021
  Total_cost_vac4_y3_dis<-Total_cost_vac4_y3*disc.2022
  Total_cost_vac4_y4_dis<-Total_cost_vac4_y4*disc.2023
  Total_cost_vac4_y5_dis<-Total_cost_vac4_y5*disc.2024
  
  
  #prepare discounted parameters
  ##discount treatment cost
  Total_Cost_soc_dis <- array(NA, dim=c(4,5,1000))
  Total_Cost_gov_dis <- array(NA, dim=c(4,5,1000))
  
  Total_Cost_vac1_soc_dis <- array(NA, dim=c(4,5,1000)) 
  Total_Cost_vac1_gov_dis <- array(NA, dim=c(4,5,1000)) 
  
  Total_QALY_loss_novac_dis <- array(NA, dim=c(4,5,1000))
  Total_QALY_loss_vac1_dis <- array(NA, dim=c(4,5,1000)) 
  
  for (k in 1:1000) 
  {
    Total_Cost_soc_dis[,1,k]<-Total_Cost_soc[,1,k]*disc.2020
    Total_Cost_soc_dis[,2,k]<-Total_Cost_soc[,2,k]*disc.2021
    Total_Cost_soc_dis[,3,k]<-Total_Cost_soc[,3,k]*disc.2022
    Total_Cost_soc_dis[,4,k]<-Total_Cost_soc[,4,k]*disc.2023
    Total_Cost_soc_dis[,5,k]<-Total_Cost_soc[,5,k]*disc.2024
    
    Total_Cost_gov_dis[,1,k]<-Total_Cost_gov[,1,k]*disc.2020
    Total_Cost_gov_dis[,2,k]<-Total_Cost_gov[,2,k]*disc.2021
    Total_Cost_gov_dis[,3,k]<-Total_Cost_gov[,3,k]*disc.2022
    Total_Cost_gov_dis[,4,k]<-Total_Cost_gov[,4,k]*disc.2023
    Total_Cost_gov_dis[,5,k]<-Total_Cost_gov[,5,k]*disc.2024
    
    Total_Cost_vac1_soc_dis[,1,k]<-Total_Cost_vac1_soc[,1,k]*disc.2020
    Total_Cost_vac1_soc_dis[,2,k]<-Total_Cost_vac1_soc[,2,k]*disc.2021
    Total_Cost_vac1_soc_dis[,3,k]<-Total_Cost_vac1_soc[,3,k]*disc.2022
    Total_Cost_vac1_soc_dis[,4,k]<-Total_Cost_vac1_soc[,4,k]*disc.2023
    Total_Cost_vac1_soc_dis[,5,k]<-Total_Cost_vac1_soc[,5,k]*disc.2024
    
    Total_Cost_vac1_gov_dis[,1,k]<-Total_Cost_vac1_gov[,1,k]*disc.2020
    Total_Cost_vac1_gov_dis[,2,k]<-Total_Cost_vac1_gov[,2,k]*disc.2021
    Total_Cost_vac1_gov_dis[,3,k]<-Total_Cost_vac1_gov[,3,k]*disc.2022
    Total_Cost_vac1_gov_dis[,4,k]<-Total_Cost_vac1_gov[,4,k]*disc.2023
    Total_Cost_vac1_gov_dis[,5,k]<-Total_Cost_vac1_gov[,5,k]*disc.2024
    
    Total_QALY_loss_novac_dis[,1,k]<-Total_QALY_loss_novac[,1,k]*disc.2020
    Total_QALY_loss_novac_dis[,2,k]<-Total_QALY_loss_novac[,2,k]*disc.2021
    Total_QALY_loss_novac_dis[,3,k]<-Total_QALY_loss_novac[,3,k]*disc.2022
    Total_QALY_loss_novac_dis[,4,k]<-Total_QALY_loss_novac[,4,k]*disc.2023
    Total_QALY_loss_novac_dis[,5,k]<-Total_QALY_loss_novac[,5,k]*disc.2024
    
    Total_QALY_loss_vac1_dis[,1,k]<-Total_QALY_loss_vac1[,1,k]*disc.2020
    Total_QALY_loss_vac1_dis[,2,k]<-Total_QALY_loss_vac1[,2,k]*disc.2021
    Total_QALY_loss_vac1_dis[,3,k]<-Total_QALY_loss_vac1[,3,k]*disc.2022
    Total_QALY_loss_vac1_dis[,4,k]<-Total_QALY_loss_vac1[,4,k]*disc.2023
    Total_QALY_loss_vac1_dis[,5,k]<-Total_QALY_loss_vac1[,5,k]*disc.2024
  }
  
  
  
  Total_Cost_soc_dis
  Total_Cost_gov_dis
  Total_Cost_vac1_soc_dis
  Total_Cost_vac1_gov_dis
  
  ##chooss gov or soc
  Total_Cost_treat_dis<-Total_Cost_soc_dis
  ##chooss gov or soc
  Total_Cost_vac1_treat_dis<-Total_Cost_vac1_soc_dis
  #discount QALY without vac
  
  Total_QALY_loss_novac_dis
  Total_QALY_loss_vac1_dis
  
  ###Finish Setting up parameter distribution############
  Outcomes_withoutvac <- array(NA, dim=c(4,40,1000)) 
  Outcomes_vac1 <- array(0, dim=c(4,40,1000)) 
  Outcomes_vac2 <- array(0, dim=c(4,40,1000)) 
  Outcomes_vac3 <- array(0, dim=c(4,40,1000)) 
  Outcomes_vac4 <- array(0, dim=c(4,40,1000)) 
  
  
  for (k in 1:1000) # k is number of iterations
  {
    Outcomes_withoutvac[1:4,1:5,k]<-Total_case[1:4,1:5,k]
    Outcomes_withoutvac[1:4,6:10,k]<-Total_case_IPD[1:4,1:5,k]
    Outcomes_withoutvac[1:4,11:15,k]<-Total_case_OPD[1:4,1:5,k]
    Outcomes_withoutvac[1:4,16:20,k]<-Total_case_self[1:4,1:5,k]
    Outcomes_withoutvac[1:4,21:25,k]<-0
    Outcomes_withoutvac[1:4,26:30,k]<-Total_Cost_treat_dis[1:4,1:5,k]
    Outcomes_withoutvac[1:4,31:35,k]<-Total_Cost_treat_dis[1:4,1:5,k]+0
    Outcomes_withoutvac[1:4,36:40,k]<-Total_QALY_loss_novac_dis[1:4,1:5,k]
    
    Outcomes_vac1[1:4,1:5,k]<-Total_case_vac1[1:4,1:5,k]
    Outcomes_vac1[1:4,6:10,k]<-Total_case_IPD_vac1[1:4,1:5,k]
    Outcomes_vac1[1:4,11:15,k]<-Total_case_OPD_vac1[1:4,1:5,k]
    Outcomes_vac1[1:4,16:20,k]<-Total_case_self_vac1[1:4,1:5,k]
    Outcomes_vac1[1,21,k]<-Total_cost_vac1_y1_dis[k]
    Outcomes_vac1[1,22,k]<-Total_cost_vac1_y2_dis[k]
    Outcomes_vac1[1,23,k]<-Total_cost_vac1_y3_dis[k]
    Outcomes_vac1[1,24,k]<-Total_cost_vac1_y4_dis[k]
    Outcomes_vac1[1,25,k]<-Total_cost_vac1_y5_dis[k]
    Outcomes_vac1[1:4,26:30,k]<-Total_Cost_vac1_treat_dis[1:4,1:5,k]
    Outcomes_vac1[1:4,31:35,k]<-Total_Cost_vac1_treat_dis[1:4,1:5,k]+Outcomes_vac1[1:4,21:25,k] #to add vaccine cost
    Outcomes_vac1[1:4,36:40,k]<-Total_QALY_loss_vac1_dis[1:4,1:5,k]
    
    Outcomes_vac2[1:4,1:5,k]<-Total_case_vac1[1:4,1:5,k]
    Outcomes_vac2[1:4,6:10,k]<-Total_case_IPD_vac1[1:4,1:5,k]
    Outcomes_vac2[1:4,11:15,k]<-Total_case_OPD_vac1[1:4,1:5,k]
    Outcomes_vac2[1:4,16:20,k]<-Total_case_self_vac1[1:4,1:5,k]
    Outcomes_vac2[1,21,k]<-Total_cost_vac2_y1_dis[k]
    Outcomes_vac2[1,22,k]<-Total_cost_vac2_y2_dis[k]
    Outcomes_vac2[1,23,k]<-Total_cost_vac2_y3_dis[k]
    Outcomes_vac2[1,24,k]<-Total_cost_vac2_y4_dis[k]
    Outcomes_vac2[1,25,k]<-Total_cost_vac2_y5_dis[k]
    Outcomes_vac2[1:4,26:30,k]<-Total_Cost_vac1_treat_dis[1:4,1:5,k]
    Outcomes_vac2[1:4,31:35,k]<-Total_Cost_vac1_treat_dis[1:4,1:5,k]+Outcomes_vac2[1:4,21:25,k]
    Outcomes_vac2[1:4,36:40,k]<-Total_QALY_loss_vac1_dis[1:4,1:5,k]
    
    Outcomes_vac3[1:4,1:5,k]<-Total_case_vac1[1:4,1:5,k]
    Outcomes_vac3[1:4,6:10,k]<-Total_case_IPD_vac1[1:4,1:5,k]
    Outcomes_vac3[1:4,11:15,k]<-Total_case_OPD_vac1[1:4,1:5,k]
    Outcomes_vac3[1:4,16:20,k]<-Total_case_self_vac1[1:4,1:5,k]
    Outcomes_vac3[1,21,k]<-Total_cost_vac3_y1_dis[k]
    Outcomes_vac3[1,22,k]<-Total_cost_vac3_y2_dis[k]
    Outcomes_vac3[1,23,k]<-Total_cost_vac3_y3_dis[k]
    Outcomes_vac3[1,24,k]<-Total_cost_vac3_y4_dis[k]
    Outcomes_vac3[1,25,k]<-Total_cost_vac3_y5_dis[k]
    Outcomes_vac3[1:4,26:30,k]<-Total_Cost_vac1_treat_dis[1:4,1:5,k]
    Outcomes_vac3[1:4,31:35,k]<-Total_Cost_vac1_treat_dis[1:4,1:5,k]+Outcomes_vac3[1:4,21:25,k]
    Outcomes_vac3[1:4,36:40,k]<-Total_QALY_loss_vac1_dis[1:4,1:5,k]
    
    Outcomes_vac4[1:4,1:5,k]<-Total_case_vac1[1:4,1:5,k]
    Outcomes_vac4[1:4,6:10,k]<-Total_case_IPD_vac1[1:4,1:5,k]
    Outcomes_vac4[1:4,11:15,k]<-Total_case_OPD_vac1[1:4,1:5,k]
    Outcomes_vac4[1:4,16:20,k]<-Total_case_self_vac1[1:4,1:5,k]
    Outcomes_vac4[1,21,k]<-Total_cost_vac4_y1_dis[k]
    Outcomes_vac4[1,22,k]<-Total_cost_vac4_y2_dis[k]
    Outcomes_vac4[1,23,k]<-Total_cost_vac4_y3_dis[k]
    Outcomes_vac4[1,24,k]<-Total_cost_vac4_y4_dis[k]
    Outcomes_vac4[1,25,k]<-Total_cost_vac4_y5_dis[k]
    Outcomes_vac4[1:4,26:30,k]<-Total_Cost_vac1_treat_dis[1:4,1:5,k]
    Outcomes_vac4[1:4,31:35,k]<-Total_Cost_vac1_treat_dis[1:4,1:5,k]+Outcomes_vac4[1:4,21:25,k]
    Outcomes_vac4[1:4,36:40,k]<-Total_QALY_loss_vac1_dis[1:4,1:5,k]
  }
  
  rownames(Outcomes_withoutvac) <- c("age0-5","age6-14","age15-64","age>64")  
  colnames(Outcomes_withoutvac) <- c("Total_case2020","Total_case2021","Total_case2022","Total_case2023","Total_case2024","Total_case_IPD2020","Total_case_IPD2021","Total_case_IPD2022","Total_case_IPD2023","Total_case_IPD2024","Total_case_OPD2020","Total_case_OPD2021","Total_case_OPD2022","Total_case_OPD2023","Total_case_OPD2024","Total_case_self2020","Total_case_self2021","Total_case_self2022","Total_case_self2023","Total_case_self2024","cost_vac2020","cost_vac2021","cost_vac2022","cost_vac2023","cost_vac2024","cost_treat2020","cost_treat2021","cost_treat2022","cost_treat2023","cost_treat2024","Total_cost2020","Total_cost2021","Total_cost2022","Total_cost2023","Total_cost2024","QALY_loss2020","QALY_loss2021","QALY_loss2022","QALY_loss2023","QALY_loss2024" )
  rownames(Outcomes_vac1) <- c("age0-5","age6-14","age15-64","age>64")  
  colnames(Outcomes_vac1) <- c("Total_case2020","Total_case2021","Total_case2022","Total_case2023","Total_case2024","Total_case_IPD2020","Total_case_IPD2021","Total_case_IPD2022","Total_case_IPD2023","Total_case_IPD2024","Total_case_OPD2020","Total_case_OPD2021","Total_case_OPD2022","Total_case_OPD2023","Total_case_OPD2024","Total_case_self2020","Total_case_self2021","Total_case_self2022","Total_case_self2023","Total_case_self2024","cost_vac2020","cost_vac2021","cost_vac2022","cost_vac2023","cost_vac2024","cost_treat2020","cost_treat2021","cost_treat2022","cost_treat2023","cost_treat2024","Total_cost2020","Total_cost2021","Total_cost2022","Total_cost2023","Total_cost2024","QALY_loss2020","QALY_loss2021","QALY_loss2022","QALY_loss2023","QALY_loss2024" )
  rownames(Outcomes_vac2) <- c("age0-5","age6-14","age15-64","age>64")  
  colnames(Outcomes_vac2) <- c("Total_case2020","Total_case2021","Total_case2022","Total_case2023","Total_case2024","Total_case_IPD2020","Total_case_IPD2021","Total_case_IPD2022","Total_case_IPD2023","Total_case_IPD2024","Total_case_OPD2020","Total_case_OPD2021","Total_case_OPD2022","Total_case_OPD2023","Total_case_OPD2024","Total_case_self2020","Total_case_self2021","Total_case_self2022","Total_case_self2023","Total_case_self2024","cost_vac2020","cost_vac2021","cost_vac2022","cost_vac2023","cost_vac2024","cost_treat2020","cost_treat2021","cost_treat2022","cost_treat2023","cost_treat2024","Total_cost2020","Total_cost2021","Total_cost2022","Total_cost2023","Total_cost2024","QALY_loss2020","QALY_loss2021","QALY_loss2022","QALY_loss2023","QALY_loss2024" )
  rownames(Outcomes_vac3) <- c("age0-5","age6-14","age15-64","age>64")  
  colnames(Outcomes_vac3) <- c("Total_case2020","Total_case2021","Total_case2022","Total_case2023","Total_case2024","Total_case_IPD2020","Total_case_IPD2021","Total_case_IPD2022","Total_case_IPD2023","Total_case_IPD2024","Total_case_OPD2020","Total_case_OPD2021","Total_case_OPD2022","Total_case_OPD2023","Total_case_OPD2024","Total_case_self2020","Total_case_self2021","Total_case_self2022","Total_case_self2023","Total_case_self2024","cost_vac2020","cost_vac2021","cost_vac2022","cost_vac2023","cost_vac2024","cost_treat2020","cost_treat2021","cost_treat2022","cost_treat2023","cost_treat2024","Total_cost2020","Total_cost2021","Total_cost2022","Total_cost2023","Total_cost2024","QALY_loss2020","QALY_loss2021","QALY_loss2022","QALY_loss2023","QALY_loss2024" )
  rownames(Outcomes_vac4) <- c("age0-5","age6-14","age15-64","age>64")  
  colnames(Outcomes_vac4) <- c("Total_case2020","Total_case2021","Total_case2022","Total_case2023","Total_case2024","Total_case_IPD2020","Total_case_IPD2021","Total_case_IPD2022","Total_case_IPD2023","Total_case_IPD2024","Total_case_OPD2020","Total_case_OPD2021","Total_case_OPD2022","Total_case_OPD2023","Total_case_OPD2024","Total_case_self2020","Total_case_self2021","Total_case_self2022","Total_case_self2023","Total_case_self2024","cost_vac2020","cost_vac2021","cost_vac2022","cost_vac2023","cost_vac2024","cost_treat2020","cost_treat2021","cost_treat2022","cost_treat2023","cost_treat2024","Total_cost2020","Total_cost2021","Total_cost2022","Total_cost2023","Total_cost2024","QALY_loss2020","QALY_loss2021","QALY_loss2022","QALY_loss2023","QALY_loss2024" )
  
  Outcomes_withoutvac
  Outcomes_vac1
  Outcomes_vac2
  Outcomes_vac3
  Outcomes_vac4
  
  
  # comparison all vac compare with No vaccine (Outputs: Case averted, Incremental Cost, QALY gained)
  # calculate mean withour vac
  #####################
  
  #######################################
  #mean case of all strategies from 1000 iteration
  #######################################
  Outcomes_withoutvac_mean <- matrix(NA,4,40)
  Outcomes_vac1_mean <- matrix(NA,4,40)
  Outcomes_vac2_mean <- matrix(NA,4,40)
  Outcomes_vac3_mean <- matrix(NA,4,40)
  Outcomes_vac4_mean <- matrix(NA,4,40)
  for (i in 1:4) 
    for (j in 1:40) 
    {
      Outcomes_withoutvac_mean[i,j] <- mean(Outcomes_withoutvac[i,j,1:1000])
      Outcomes_vac1_mean[i,j] <- mean(Outcomes_vac1[i,j,1:1000])
      Outcomes_vac2_mean[i,j] <- mean(Outcomes_vac2[i,j,1:1000])
      Outcomes_vac3_mean[i,j] <- mean(Outcomes_vac3[i,j,1:1000])
      Outcomes_vac4_mean[i,j] <- mean(Outcomes_vac4[i,j,1:1000])
    }
  
  Outcomes_withoutvac_mean 
  Outcomes_vac1_mean
  Outcomes_vac2_mean
  Outcomes_vac3_mean
  Outcomes_vac4_mean
  #write.csv(Outcomes_withoutvac_mean, file = "Outcomes_all_novac_mean_1_byage_150919.csv")
  #write.csv(Outcomes_vac1_mean, file = "Outcomes_all_vac1mean_1_byage_150919.csv")
  #write.csv(Outcomes_vac2_mean, file = "Outcomes_all_vac2mean_1_byage_150919.csv")
  #write.csv(Outcomes_vac3_mean, file = "Outcomes_all_vac3mean_1_byage_150919.csv")
  #write.csv(Outcomes_vac4_mean, file = "Outcomes_all_vac4mean_1_byage_150919.csv")
  
  
  colSums(Outcomes_withoutvac_mean )
  colSums(Outcomes_vac1_mean)
  
  colSums(Outcomes_withoutvac_mean )-colSums(Outcomes_vac2_mean )
  
  
  #mean outcome table vac 1-4
  summary_outcome_mean <- matrix(NA,5,35)
  summary_outcome_mean[1,1:20] <- colSums(Outcomes_withoutvac_mean)[1:20]
  summary_outcome_mean[1,21:35] <- colSums(Outcomes_withoutvac_mean)[26:40]
  summary_outcome_mean[2,1:20] <- colSums(Outcomes_vac1_mean)[1:20]
  summary_outcome_mean[2,21:35] <- colSums(Outcomes_vac1_mean)[26:40]
  summary_outcome_mean[3,1:20] <- colSums(Outcomes_vac2_mean)[1:20]
  summary_outcome_mean[3,21:35] <- colSums(Outcomes_vac2_mean)[26:40]
  summary_outcome_mean[4,1:20] <- colSums(Outcomes_vac3_mean)[1:20]
  summary_outcome_mean[4,21:35] <- colSums(Outcomes_vac3_mean)[26:40]
  summary_outcome_mean[5,1:20] <- colSums(Outcomes_vac4_mean)[1:20]
  summary_outcome_mean[5,21:35] <- colSums(Outcomes_vac4_mean)[26:40]
  summary_outcome_mean
  
  rownames(summary_outcome_mean) <- c("no_vac","vac1","vac2","vac3","vac4")  
  colnames(summary_outcome_mean) <- c("Total_case2020","Total_case2021","Total_case2022","Total_case2023","Total_case2024",
                                      "Total_case_IPD2020","Total_case_IPD2021","Total_case_IPD2022","Total_case_IPD2023","Total_case_IPD2024",
                                      "Total_case_OPD2020","Total_case_OPD2021","Total_case_OPD2022","Total_case_OPD2023","Total_case_OPD2024",
                                      "Total_case_self2020","Total_case_self2021","Total_case_self2022","Total_case_self2023","Total_case_self2024",
                                      "cost2020","cost2021","cost2022","cost2023","cost2024",
                                      "Total_cost2020","Total_cost2021","Total_cost2022","Total_cost2023","Total_cost2024",
                                      "QALY_loss2020","QALY_loss2021","QALY_loss2022","QALY_loss2023","QALY_loss2024")
  summary_outcome_mean #basecase
  #write.csv(summary_outcome_mean, file = "Outcomes_all_mean_1_gov_BIA130919.csv")
  #write.csv(summary_outcome_mean, file = "Outcomes_all_mean_BIA170919.csv")
  #mean over 5 years
  mean_cost_novac<-mean(summary_outcome_mean[1,26:30])   #mean cost no vac
  mean_QALY_novac<-mean(summary_outcome_mean[1,31:35])   #mean QALY no vac
  
  mean_cost_vac1_base<-mean(summary_outcome_mean[2,26:30])   #mean cost vac
  mean_QALY_vac1_base<-mean(summary_outcome_mean[2,31:35])   #mean QALY vac
  
  mean_cost_vac2_base<-mean(summary_outcome_mean[3,26:30])   #mean cost vac
  mean_QALY_vac2_base<-mean(summary_outcome_mean[3,31:35])   #mean QALY vac
  
  mean_cost_vac3_base<-mean(summary_outcome_mean[4,26:30])   #mean cost vac
  mean_QALY_vac3_base<-mean(summary_outcome_mean[4,31:35])   #mean QALY vac
  
  mean_cost_vac4_base<-mean(summary_outcome_mean[5,26:30])   #mean cost vac
  mean_QALY_vac4_base<-mean(summary_outcome_mean[5,31:35])   #mean QALY vac
  
  ##ICER_summary mean vac 1-4
  summary_ICER_mean <- matrix(NA,4,35)
  summary_ICER_mean[1,1:20] <- summary_outcome_mean[1,1:20]-summary_outcome_mean[2,1:20]
  summary_ICER_mean[1,21:25] <- summary_outcome_mean[2,26:30]-summary_outcome_mean[1,26:30]  #cost vac- cost no vac
  summary_ICER_mean[1,26:30] <- summary_outcome_mean[1,31:35]-summary_outcome_mean[2,31:35]  #QALY no vac- QALY vac
  summary_ICER_mean[1,31:35] <- summary_ICER_mean[1,21:25]/summary_ICER_mean[1,26:30]
  
  summary_ICER_mean[2,1:20] <- summary_outcome_mean[1,1:20]-summary_outcome_mean[3,1:20]
  summary_ICER_mean[2,21:25] <- summary_outcome_mean[3,26:30]-summary_outcome_mean[1,26:30]  #cost vac- cost no vac
  summary_ICER_mean[2,26:30] <- summary_outcome_mean[1,31:35]-summary_outcome_mean[3,31:35]  #QALY no vac- QALY vac
  summary_ICER_mean[2,31:35] <- summary_ICER_mean[2,21:25]/summary_ICER_mean[2,26:30]
  
  summary_ICER_mean[3,1:20] <- summary_outcome_mean[1,1:20]-summary_outcome_mean[4,1:20]
  summary_ICER_mean[3,21:25] <- summary_outcome_mean[4,26:30]-summary_outcome_mean[1,26:30]  #cost vac- cost no vac
  summary_ICER_mean[3,26:30] <- summary_outcome_mean[1,31:35]-summary_outcome_mean[4,31:35]  #QALY no vac- QALY vac
  summary_ICER_mean[3,31:35] <- summary_ICER_mean[3,21:25]/summary_ICER_mean[3,26:30]
  
  summary_ICER_mean[4,1:20] <- summary_outcome_mean[1,1:20]-summary_outcome_mean[5,1:20]
  summary_ICER_mean[4,21:25] <- summary_outcome_mean[5,26:30]-summary_outcome_mean[1,26:30]  #cost vac- cost no vac
  summary_ICER_mean[4,26:30] <- summary_outcome_mean[1,31:35]-summary_outcome_mean[5,31:35]  #QALY no vac- QALY vac
  summary_ICER_mean[4,31:35] <- summary_ICER_mean[4,21:25]/summary_ICER_mean[4,26:30]
  rownames(summary_ICER_mean) <- c("vac1","vac2","vac3","vac4") 
  colnames(summary_ICER_mean) <- c("Averted_case2020","Averted_case2021","Averted_case2022","Averted_case2023","Averted_case2024",
                                   "Averted_case_IPD2020","Averted_case_IPD2021","Averted_case_IPD2022","Averted_case_IPD2023","Averted_case_IPD2024",
                                   "Averted_case_OPD2020","Averted_case_OPD2021","Averted_case_OPD2022","Averted_case_OPD2023","Averted_case_OPD2024",
                                   "Averted_case_self2020","Averted_case_self2021","Averted_case_self2022","Averted_case_self2023","Averted_case_self2024",
                                   "In_cost2020","In_cost2021","In_cost2022","In_cost2023","In_cost2024",
                                   "In_QALY2020","In_QALY2021","In_QALY2022","In_QALY2023","In_QALY2024",
                                   "ICER2020","ICER2021","ICER2022","ICER2023","ICER2024")
  summary_ICER_mean
  
  
  ##mean 5 years
  mean_Incost_vac1_base<-mean(summary_ICER_mean[1,21:25])  #vac1
  mean_InQALY_vac1_base<-mean(summary_ICER_mean[1,26:30])  #vac1
  mean_ICER_vac1_base<-mean(summary_ICER_mean[1,31:35])  #vac1
  
  mean_Incost_vac2_base<-mean(summary_ICER_mean[2,21:25])  #vac2
  mean_InQALY_vac2_base<-mean(summary_ICER_mean[2,26:30])  #vac2
  mean_ICER_vac2_base<-mean(summary_ICER_mean[2,31:35])  #vac2
  # This is what we need at the end
  mean_Incost_vac3_base<-mean(summary_ICER_mean[3,21:25])  #vac3
  mean_InQALY_vac3_base<-mean(summary_ICER_mean[3,26:30])  #vac3
  mean_ICER_vac3_base<-mean(summary_ICER_mean[3,31:35])  #vac3
  
  mean_Incost_vac4_base<-mean(summary_ICER_mean[4,21:25])  #vac4
  mean_InQALY_vac4_base<-mean(summary_ICER_mean[4,26:30])  #vac4
  mean_ICER_vac4_base<-mean(summary_ICER_mean[4,31:35])  #vac4
  
  #compare ROTAVAC with No vaccine (may do another one with Rotarix or Rotateq)
  #CQICER<-c(mean_Incost_vac3_base,mean_InQALY_vac3_base,mean_ICER_vac3_base)
  CQICER<-c(mean_Incost_vac4_base,mean_InQALY_vac4_base,mean_ICER_vac4_base)
  
  #change to vac4 too for SIIL
  return(CQICER)
}
#par1=Vaccine Cost
#par2=VE 90-30
TWSA(par1[1],par2[[1]])
TWSA(par1[1],par2[[2]])






