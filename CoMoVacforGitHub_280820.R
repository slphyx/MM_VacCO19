###########################################################################
## AGE-DEPENDANT SEIRS MODEL WITH 5-YEAR AGE CLASSES USING UN DEMOG DATA ##
## Modified for Thailand by sompob@tropmedres.ac (8 May 2020)
## Vaccine compartment added by nantasit@tropmedres.ac (22 July 2020)
###########################################################################

#ScenarioIII: Low burden peak at ~100-200 cases
rm(list=ls()) 
#setwd("~/Documents/MORU 2017/COVID19_modelling/CoMO/Model/TH_model")
setwd("~/Documents/MORU 2017/COVID19_modelling/CoMO/Model/TH_model")
library("deSolve")
library("dplyr")


#########  INCIDENCE DATA
incdata_X<-read.csv("Thcovidcases.csv")
incdata_X[,1]<-as.Date(incdata_X[,1],"%d-%m-%y")


########## POP Structure
######
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
startdate<-as.Date("2021-01-01") 
# stopdate<-Sys.Date() # today
stopdate<-as.Date("2025-12-31") #defult ("2020-12-31)
# stopdate<-as.Date("2020-03-18")
day_start <- as.numeric(startdate-startdate)
day_stop <- as.numeric(stopdate-startdate)
times <- seq(day_start, day_stop)
times.date <- startdate+times
tin<-as.numeric(startdate-as.Date("2021-01-01"))/365.25
initP<-sum(popstruc[,2])       # population size 
ageindcase<-20                 # age of index case (years)
aci <- floor((ageindcase/5)+1) # age class of index case


######    THESE ARE JUST VARIABLE DEFINITIONS - PLEASE DO NOT CHANGE   #################################
######     GO TO LINE 652 TO Choose Interventions' start dates         #################################
######
# date to start the self isolation intervention
#date_work_on<-as.Date("2121-12-15")
# date to start the school closure
#date_school_on<-as.Date("2020-03-23")
#date_school_on<-as.Date("2121-12-15")
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
date_lockdown_mid_on<-as.Date("2121-03-23")
#date_lockdown_mid_on<-as.Date("2021-12-18")

date_selfis_on<-as.Date("2121-08-16") #Selfisolation and social distancing as lockdown key intervention setting with different timing of implementation correspond to differnt magnitude

#######To adjust background intervention and vaccination timing########
date_hand_on<-as.Date("2021-1-01") #Hand washing as baseline intervnetion
date_dist_on<-as.Date("2021-08-16") #High:"2021-11-25", Medium:""2021-10-1" and Low:"2021-08-16"
date_work_on<-as.Date("2021-08-16")
# date to start the school closure
#date_school_on<-as.Date("2020-03-23")
date_school_on<-as.Date("2021-08-16")

date_vaccine_on<-as.Date("2021-05-01")
######

# DEFINE Vaccine Strategy
######

#Define vaccination strategy based on 21 age group
NoVac<-rep(0,21)    #No vaccine
AgeAll<-rep(1,21)   #Vac all age gr
AgeHighI<-rep(0,21)
AgeHighI[5:8]<-1    #Vac age 20-39
AgeAdult<-rep(1,21)
AgeAdult[1:3]<-0    #Vac age >15
AgeElder<-rep(0,21)  
AgeElder[13:21]<-1  #Vac age >60

#choose Vaccination strategy
AgeVac<-AgeAll   #5 Options; NoVac, AgeAll, AgeHighI, AgeAdult, AgeElder
######

# DEFINE PARAMETERS
#####
parameters <- c(
  # Transmission instrinsic
  p=0.042/2,           # probabilty of infection given a contact 
  rho = 25,          # relative infectiousness of incubation phase (%) min 0 max 100 step 0.5 
  omega=5,         # average duration of immunity (years) min 0.5 max 100 step 0.5  (default 200)
  omegav=1,        # average duration of vaccinated immunity (years) min 0.5 max 100 step 0.5
  gamma=3.5,         # average incubation period (days) min 1 max 7 step 0.5 
  nui=4.5,             # average duration of symptomatic infection period (days) min 1 max 7 step 0.5
  report=50,          # percentage of all asymptomatic infections that are reported (%) min 0 max 100 step 1
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
  pclin=60,          # default=15 probability upon infection of developing clinical symptoms
  prob_icu = 40,     # default=70 probability upon hospitalisation of requiring icu admission   
  prob_vent = 70,    # default=80 probability upon admission to the UCI of requiring a ventilator
  # INTERVENTIONS
  # vaccination
  vaccine_on= as.numeric(date_vaccine_on-startdate),
  vaccine_eff1=0,   # vaccine efficacy (%)- min 0 max 100 step 1; reduce infection
  vaccine_eff2=0,   # vaccine efficacy (%)- min 0 max 100 step 1; reduce transmission
  vaccine_eff3=0,   # vaccine efficacy (%)- min 0 max 100 step 1; reduce severity
  vaccine_cov=0,    # vaccine coverage (%)- min 0 max 100 step 1
  vac_campaign=8, # Number of weeks it takes to reach maximum coverage - min 1 max 8 step 1
  
  # self isolation
  selfis_on=as.numeric(date_selfis_on-startdate),
  selfis_dur=200,    # duration of self-isolation protocol (weeks) min 1 max 52 step 1
  selfis_cov=70,    # coverage of self isolation (%) min 0 max 100 step 1
  selfis_eff=70,    # adherence to self isolation (%) min 0 max 100 step 1
  # social distancing
  dist_on=as.numeric(date_dist_on-startdate),
  dist_dur=260,      # duration of social distancing protocol (weeks) min 1 max 52 step 1
  dist_cov=100,      # coverage of social distancing (%) min 0 max 100 step 1
  dist_eff=70,     # adherence to social distancing (%) min 0 max 100 step 1
  # hand washing
  hand_on=as.numeric(date_hand_on-startdate),
  hand_dur=300,      # duration of increased hand hygiene protocol (weeks) min 1 max 52 step 1
  hand_eff=30,       # efficacy of hand hygiene  (%) min 0 max 100 step 1 (default 5)
  
  
  # working at home
  work_on=as.numeric(date_work_on-startdate),
  work_dur=8,      # duration of working from home protocol (weeks) min 1 max 52 step 1
  work_cov=50,      # coverage of working from home (%) min 0 max 100 step 1
  work_eff=80,      # efficacy of working from home (%) min 0 max 100 step 1
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
  # imported cases 
  mean_imports = 0,           # user defined - mean number of infectious migrants per day (number) - min 0 max 500 step 1
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

# Scale parameters to percentages/ rates
####################
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

######
# Define the indices for each variable
###########################################################################

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

######
# MODEL INITIAL CONDITIONS
###########################################################################

initI<-0*popstruc[,2]  # Infected and symptomatic
initE<-0*popstruc[,2]  # Incubating
initE[aci]<-1          # place random index case in E compartment
initR<-0.035*popstruc[,2]  # Immune (seroprevalence from china 3.2-3.8% of the population already had immune)Xu X, Sun J, Nie S, et al. Seroprevalence of immunoglobulin M and G antibodies against SARS-CoV-2 in China [published online ahead of print, 2020 Jun 5]. Nat Med. 2020;10.1038/s41591-020-0949-6. doi:10.1038/s41591-020-0949-6
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


######
# MODEL SOLVING EQUARIONS
#####
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
             dist<-0.35 #default 0.35
             school<-0.85
             trvban_eff<-0
             quarantine_rate<-0.05
             work<-0.5
             cocoon<-0.95
             hand<-0.05 #defalt 0.05
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
         dSdt <- -S*lam-S*vaccinate*AgeVac+omega*R+ageing%*%S-mort*S+birth-quarantine_rate*S +(1/quarantine_days)*QS+omegav*SV+omega*RV
         dEdt <- S*lam-gamma*E+ageing%*%E-mort*E-quarantine_rate*E+(1/quarantine_days)*QE 
         dIdt <- gamma*(1-pclin)*(1-screen_eff)*(1-ihr[,2])*E-nui*I+ageing%*%I-mort*I + (1/quarantine_days)*QI - quarantine_rate*I
         dCLdt<- gamma*pclin*(1-selfis)*(1-ihr[,2])*E-nui*CL+ageing%*%CL-mort*CL + (1/quarantine_days)*QC
         dRdt <- nui*I-omega*R+nui*X+nui*CL+ageing%*%R-mort*R + (1/quarantine_days)*QR + nus*(1-pdeath_h*ifr[,2])*H + (1-pdeath_icu*ifr[,2])*nu_icu*ICU + (1-pdeath_icuc*ifr[,2])*nu_icuc*ICUC + (1-pdeath_hc*ifr[,2])*nusc*HC + (1-pdeath_vent*ifr[,2])*nu_vent*Vent+ (1-pdeath_ventc*ifr[,2])*nu_ventc*VentC - vaccinate*R*AgeVac + omegav*RV
         dXdt <- gamma*selfis*pclin*(1-ihr[,2])*E+gamma*(1-pclin)*screen_eff*(1-ihr[,2])*E-nui*X+ageing%*%X-mort*X 
         ############
         #dVdt <- vaccinate*S -(1-vaccine_eff)*lam*V +ageing%*%V - mort*V
         
         #to add terms (vaccine_eff1=%reduction in infection, vaccine_eff2=%reduction in duration of infection, vaccine_eff3=%reduction in risk of sevrerity, hospitalisation in this case, but we have hospi, ICU and ICUVent [need to clarify the effect of vac])
         #SV, EV, IV, CL,HV,ICUV, VentV, RV set initial condition
         #
         #if (vaccine){cf
         #Add vaccine compartment - add (AgeVac*)S*vaccinate to dSdt and dSVdt
         dSVdt <- S*vaccinate*AgeVac - (1-vaccine_eff1)*SV*lam + ageing%*%SV-mort*SV-omegav*SV #Assuming the lam is the same as general population
         
         dEVdt <- (1-vaccine_eff1)*SV*lam - gamma*EV +ageing%*%EV-mort*EV
         dIVdt <- gamma*(1-pclin)*(1-(ihr[,2]*(1-vaccine_eff3)))*EV-nui*(1/(1-vaccine_eff2))*IV+ageing%*%IV - mort*IV
         dCLVdt<- gamma*(pclin)*(1-(ihr[,2]*(1-vaccine_eff3)))*EV-nui*(1/(1-vaccine_eff2))*CLV+ageing%*%CLV - mort*CLV
         dHVdt <- gamma*ihr[,2]*(1-vaccine_eff3)*(1-prob_icu)*EV-nus*(1/(1-vaccine_eff2))*HV+ ageing%*%HV - mort*HV
         dICUVdt <- gamma*ihr[,2]*(1-vaccine_eff3)*prob_icu*(1-prob_vent)*EV-nu_icu*(1/(1-vaccine_eff2))*ICUV +ageing%*%ICUV - mort*ICUV
         dVentVdt <- gamma*ihr[,2]*(1-vaccine_eff3)*prob_icu*prob_vent*EV-nu_vent*(1/(1-vaccine_eff2))*VentV +ageing%*%VentV - mort*VentV
         dRVdt <- nui*(1/(1-vaccine_eff2))*IV + nui*(1/(1-vaccine_eff2))*CLV +nus*(1/(1-vaccine_eff2))*HV+nu_icu*(1+vaccine_eff2)*ICUV+nu_vent*(1/(1-vaccine_eff2))*VentV + ageing%*%RV-mort*RV
         + nus*(1-pdeath_h*ifr[,2])*HV + (1-pdeath_icu*ifr[,2])*nu_icu*ICU + (1-pdeath_icuc*ifr[,2])*nu_icuc*ICUC + (1-pdeath_hc*ifr[,2])*nusc*HC + (1-pdeath_vent*ifr[,2])*nu_vent*Vent+ (1-pdeath_ventc*ifr[,2])*nu_ventc*VentC-omega*RV + vaccinate*R*AgeVac - omegav*RV
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
         list(c(dSdt,dEdt,dIdt,dRdt,dXdt,dHdt,dHCdt,dCdt,dCMdt,dSVdt,dEVdt,dIVdt,dCLVdt,dHVdt,dICUVdt,dVentVdt,dRVdt,dQSdt,dQEdt,dQIdt,dQRdt,dCLdt,dQCdt,dICUdt,dICUCdt,dVentdt,dVentCdt,dCMCdt),lam)
       }
  ) 
}
#######Getting the Target Outcomes
process_ode_outcome <- function(out){
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
  
  inc1h<- parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*parameters["prob_icu"]+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*parameters["prob_icu"]+
    parameters["gamma"]*out[,(EVindex+1)]%*%ihr[,2]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(EVindex+1)]%*%ihr[,2]*parameters["prob_icu"]
  
  #Overall incidence
  dailyinc1<-rowSums(inc1)+rowSums(inc1h)      # daily incidence
  cuminc1<-colSums(inc1)+colSums(inc1h)        # cumulative incidence
  
  #incidence of each subgroup (by age group, n=21)
  #right approach with t(t(x)*y)
  incAsym<-parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*t(t(out[,(Eindex+1)])*(1-ihr[,2]))+
    parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*t(t(out[,(QEindex+1)])*(1-ihr[,2]))+
    parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*t(t(out[,(EVindex+1)])*(1-ihr[,2]))
  
  incSym<-parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*t(t(out[,(Eindex+1)])*(1-ihr[,2]))+
    parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*t(t(out[,(QEindex+1)])*(1-ihr[,2]))+
    parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*t(t(out[,(EVindex+1)])*(1-ihr[,2]))
  
  incHosp<-parameters["gamma"]*t(t(out[,(Eindex+1)])*ihr[,2])*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*t(t(out[,(QEindex+1)])*ihr[,2])*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*t(t(out[,(EVindex+1)])*ihr[,2])*(1-parameters["prob_icu"])
  
  incICU<-parameters["gamma"]*t(t(out[,(Eindex+1)])*ihr[,2])*parameters["prob_icu"]*(1-parameters["prob_vent"])+
    parameters["gamma"]*t(t(out[,(QEindex+1)])*ihr[,2])*parameters["prob_icu"]*(1-parameters["prob_vent"])+
    parameters["gamma"]*t(t(out[,(EVindex+1)])*ihr[,2])*parameters["prob_icu"]*(1-parameters["prob_vent"])
  
  incICUv<-parameters["gamma"]*t(t(out[,(Eindex+1)])*ihr[,2])*parameters["prob_icu"]*parameters["prob_vent"]+
    parameters["gamma"]*t(t(out[,(QEindex+1)])*ihr[,2])*parameters["prob_icu"]*parameters["prob_vent"]+
    parameters["gamma"]*t(t(out[,(EVindex+1)])*ihr[,2])*parameters["prob_icu"]*parameters["prob_vent"]
  
  
  DailyincAsym<-parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*out[,(Eindex+1)]%*%(1-ihr[,2])+
    parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*out[,(QEindex+1)]%*%(1-ihr[,2])+
    parameters["report"]*parameters["gamma"]*(1-parameters["pclin"])*out[,(EVindex+1)]%*%(1-ihr[,2])        # cumulative incidence
  
  DailyincSym<-parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*out[,(Eindex+1)]%*%(1-ihr[,2])+
    parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*out[,(QEindex+1)]%*%(1-ihr[,2])+
    parameters["reportc"]*parameters["gamma"]*parameters["pclin"]*out[,(EVindex+1)]%*%(1-ihr[,2])
  
  DailyincHosp<-parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*parameters["reporth"]*(1-parameters["prob_icu"])+
    parameters["gamma"]*out[,(EVindex+1)]%*%ihr[,2]*(1-parameters["prob_icu"])
  
  DailyincICU<-parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*parameters["prob_icu"]*(1-parameters["prob_vent"])+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*parameters["prob_icu"]*(1-parameters["prob_vent"])+
    parameters["gamma"]*out[,(EVindex+1)]%*%ihr[,2]*parameters["prob_icu"]*(1-parameters["prob_vent"])
  
  DailyincICUv<-parameters["gamma"]*out[,(Eindex+1)]%*%ihr[,2]*parameters["prob_icu"]*parameters["prob_vent"]+
    parameters["gamma"]*out[,(QEindex+1)]%*%ihr[,2]*parameters["prob_icu"]*parameters["prob_vent"]+
    parameters["gamma"]*out[,(EVindex+1)]%*%ihr[,2]*parameters["prob_icu"]*parameters["prob_vent"]
  
  
  Vac_pop<-sum(AgeVac*popstruc[,2])*parameters["vaccine_cov"]
  #parameters["vaccinate"]*out[,(Sindex+1)]*AgeVac + parameters["vaccinate"]*out[,(Rindex+1)]*AgeVac  #S*vaccinate*AgeVac+R*vaccinate*AgeVac
  
  
  cmortality1<-rowSums(out[,(CMindex+1)])      # cumulative mortality
  ccases1<-rowSums(out[,(Cindex+1)])           # cumulative cases
  
  ##########################    CALCULATE MORTALITY 
  cinc_mort_H1 <- cumsum(rowSums(parameters["nus"]*parameters["pdeath_h"]*(out[,(Hindex+1)]%*%ifr[,2])+ out[,(Hindex+1)]%*%mort))+
    cumsum(rowSums(parameters["nus"]*parameters["pdeath_h"]*(out[,(HVindex+1)]%*%ifr[,2])+ out[,(HVindex+1)]%*%mort))
  cinc_mort_HC1 <- cumsum(rowSums(parameters["nusc"]*parameters["pdeath_hc"]*(out[,(HCindex+1)]%*%ifr[,2]) + out[,(HCindex+1)]%*%mort))
  cinc_mort_ICU1 <- cumsum(rowSums(parameters["nu_icu"]*parameters["pdeath_icu"]*out[,(ICUindex+1)]%*%ifr[,2] + out[,(ICUindex+1)]%*%mort))+
    cumsum(rowSums(parameters["nu_icu"]*parameters["pdeath_icu"]*out[,(ICUVindex+1)]%*%ifr[,2] + out[,(ICUVindex+1)]%*%mort))
  cinc_mort_ICUC1 <- cumsum(rowSums(parameters["nu_icuc"]*parameters["pdeath_icuc"]*out[,(ICUCindex+1)]%*%ifr[,2] + out[,(ICUCindex+1)]%*%mort))
  cinc_mort_Vent1 <- cumsum(rowSums(parameters["nu_vent"]*parameters["pdeath_vent"]*out[,(Ventindex+1)]%*%ifr[,2] + out[,(Ventindex+1)]%*%mort))+
    cumsum(rowSums(parameters["nu_vent"]*parameters["pdeath_vent"]*out[,(VentVindex+1)]%*%ifr[,2] + out[,(VentVindex+1)]%*%mort))
  cinc_mort_VentC1 <- cumsum(rowSums(parameters["nu_ventc"]*parameters["pdeath_ventc"]*out[,(VentCindex+1)]%*%ifr[,2] + out[,(VentCindex+1)]%*%mort))
  #(new) death by age group overtime
  incDeath<-parameters["nus"]*parameters["pdeath_h"]*t(t(out[,(Hindex+1)])*ifr[,2])+
    parameters["nus"]*parameters["pdeath_h"]*t(t(out[,(HVindex+1)])*ifr[,2])+ 
    parameters["nusc"]*parameters["pdeath_hc"]*t(t(out[,(HCindex+1)])*ifr[,2])+ 
    parameters["nu_icu"]*parameters["pdeath_icu"]*t(t(out[,(ICUindex+1)])*ifr[,2])+
    parameters["nu_icu"]*parameters["pdeath_icu"]*t(t(out[,(ICUVindex+1)])*ifr[,2])+
    parameters["nu_icuc"]*parameters["pdeath_icuc"]*t(t(out[,(ICUCindex+1)])*ifr[,2])+
    parameters["nu_vent"]*parameters["pdeath_vent"]*t(t(out[,(Ventindex+1)])*ifr[,2])+
    parameters["nu_vent"]*parameters["pdeath_vent"]*t(t(out[,(VentVindex+1)])*ifr[,2])+
    parameters["nu_ventc"]*parameters["pdeath_ventc"]*t(t(out[,(VentCindex+1)])*ifr[,2])
  
  #(new) death combined all age groups
  DailyincDeath<-rowSums(parameters["nus"]*parameters["pdeath_h"]*t(t(out[,(Hindex+1)])*ifr[,2])+
                           parameters["nus"]*parameters["pdeath_h"]*t(t(out[,(HVindex+1)])*ifr[,2])+ 
                           parameters["nusc"]*parameters["pdeath_hc"]*t(t(out[,(HCindex+1)])*ifr[,2])+ 
                           parameters["nu_icu"]*parameters["pdeath_icu"]*t(t(out[,(ICUindex+1)])*ifr[,2])+
                           parameters["nu_icu"]*parameters["pdeath_icu"]*t(t(out[,(ICUVindex+1)])*ifr[,2])+
                           parameters["nu_icuc"]*parameters["pdeath_icuc"]*t(t(out[,(ICUCindex+1)])*ifr[,2])+
                           parameters["nu_vent"]*parameters["pdeath_vent"]*t(t(out[,(Ventindex+1)])*ifr[,2])+
                           parameters["nu_vent"]*parameters["pdeath_vent"]*t(t(out[,(VentVindex+1)])*ifr[,2])+
                           parameters["nu_ventc"]*parameters["pdeath_ventc"]*t(t(out[,(VentCindex+1)])*ifr[,2]))
  
  
  # Export in a cohesive format ----
  results <- list()
  results$time <- startdate + times  # dates
  results$cum_mortality <- round(cmortality1)  # cumulative mortality
  results$pct_total_pop_infected <- round(100 * tail(cumsum(rowSums(parameters["gamma"]*out[,(Eindex+1)])),1)/sum(popstruc[,2]), 1) +round(100 * tail(cumsum(rowSums(parameters["gamma"]*out[,(EVindex+1)])),1)/sum(popstruc[,2]), 1) # proportion of the  population that has been infected at the end of the simulation
  results$daily_incidence <- round(dailyinc1)  # daily incidence (Reported)
  results$daily_total_cases <- round(rowSums(parameters["gamma"]*out[,(Eindex+1)]+parameters["gamma"]*out[,(QEindex+1)]+parameters["gamma"]*out[,(EVindex+1)])) # daily incidence (Reported + Unreported)  # daily incidence (Reported + Unreported)
  results$cum_cases <- ccases1
  results$asymp_cases <- round(incAsym) #round(previcureq0)
  results$symp_cases <-  round(incSym) #round(previcureq01) 
  results$hospital <- round(incHosp)   #round(previcureq1)
  results$icu <- round(incICU)         #round(previcureq21)
  results$icuvent <-  round(incICUv)  #round(previcureq31)
  results$death <- round(incDeath)
  
  
  #total case wo age gr
  results$daily_asym <-round(DailyincAsym) #total case each day without age group
  results$daily_symp <-round(DailyincSym)
  results$daily_hosp <-round(DailyincHosp)
  results$daily_icu <-round(DailyincICU)
  results$daily_icuvent <-round(DailyincICUv)
  results$daily_death <-round(DailyincDeath)
  
  results$death_treated_hospital <- round(cinc_mort_H1)
  results$death_treated_icu <- round(cinc_mort_ICU1)
  results$death_treated_ventilator <- round(cinc_mort_Vent1)
  results$death_untreated_hospital <- round(cinc_mort_HC1)
  results$death_untreated_icu <- round(cinc_mort_ICUC1)
  results$death_untreated_ventilator <- round(cinc_mort_VentC1)
  results$total_deaths <- results$death_treated_hospital + results$death_treated_icu + results$death_treated_ventilator + results$death_untreated_hospital + results$death_untreated_icu + results$death_untreated_ventilator
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
  #
  totbase1<-as.data.frame(basemort_H1+basemort_HC1+basemort_ICU1+basemort_ICUC1+basemort_Vent1+basemort_VentC1+basemort_HV+basemort_ICUV+basemort_VentV)
  
  
  #Outcomes to export
  results$Outcomes_mean <- matrix(NA,21,31) #  Vac_pop + 6 outputs (vac_pop, asym, sym, hosp, icu, icuven, death) * 5 years (2020-2024) = 31 in total
  
  results$Outcomes_mean[1,1] <- Vac_pop
  results$Outcomes_mean[,2] <- rowSums(t(round(incAsym))[,2:321]) #round(previcureq0)
  results$Outcomes_mean[,3] <- rowSums(t(round(incSym))[,2:321]) #round(previcureq01) 
  results$Outcomes_mean[,4] <- rowSums(t(round(incHosp))[,2:321])   #round(previcureq1)
  results$Outcomes_mean[,5] <- rowSums(t(round(incICU))[,2:321])        #round(previcureq21)
  results$Outcomes_mean[,6] <- rowSums(t(round(incICUv))[,2:321])  #round(previcureq31)
  results$Outcomes_mean[,7] <- rowSums(t(round(incDeath))[,2:321])
  
  results$Outcomes_mean[,8] <- rowSums(t(round(incAsym))[,322:(322+365)]) #round(previcureq0)
  results$Outcomes_mean[,9] <- rowSums(t(round(incSym))[,322:(322+365)]) #round(previcureq01) 
  results$Outcomes_mean[,10] <- rowSums(t(round(incHosp))[,322:(322+365)])   #round(previcureq1)
  results$Outcomes_mean[,11] <- rowSums(t(round(incICU))[,322:(322+365)])        #round(previcureq21)
  results$Outcomes_mean[,12] <- rowSums(t(round(incICUv))[,322:(322+365)])  #round(previcureq31)
  results$Outcomes_mean[,13] <- rowSums(t(round(incDeath))[,322:(322+365)])
  
  results$Outcomes_mean[,14] <- rowSums(t(round(incAsym))[,(322+365):(322+365+365)]) #round(previcureq0)
  results$Outcomes_mean[,15] <- rowSums(t(round(incSym))[,(322+365):(322+365+365)]) #round(previcureq01) 
  results$Outcomes_mean[,16] <- rowSums(t(round(incHosp))[,(322+365):(322+365+365)])   #round(previcureq1)
  results$Outcomes_mean[,17] <- rowSums(t(round(incICU))[,(322+365):(322+365+365)])        #round(previcureq21)
  results$Outcomes_mean[,18] <- rowSums(t(round(incICUv))[,(322+365):(322+365+365)])  #round(previcureq31)
  results$Outcomes_mean[,19] <- rowSums(t(round(incDeath))[,(322+365):(322+365+365)])
  
  results$Outcomes_mean[,20] <- rowSums(t(round(incAsym))[,(322+365+365):(322+365+365+365)]) #round(previcureq0)
  results$Outcomes_mean[,21] <- rowSums(t(round(incSym))[,(322+365+365):(322+365+365+365)]) #round(previcureq01) 
  results$Outcomes_mean[,22] <- rowSums(t(round(incHosp))[,(322+365+365):(322+365+365+365)])   #round(previcureq1)
  results$Outcomes_mean[,23] <- rowSums(t(round(incICU))[,(322+365+365):(322+365+365+365)])        #round(previcureq21)
  results$Outcomes_mean[,24] <- rowSums(t(round(incICUv))[,(322+365+365):(322+365+365+365)])  #round(previcureq31)
  results$Outcomes_mean[,25] <- rowSums(t(round(incDeath))[,(322+365+365):(322+365+365+365)])
  
  results$Outcomes_mean[,26] <- rowSums(t(round(incAsym))[,(322+365+365+365):(322+365+365+365+365)]) #round(previcureq0)
  results$Outcomes_mean[,27] <- rowSums(t(round(incSym))[,(322+365+365+365):(322+365+365+365+365)]) #round(previcureq01) 
  results$Outcomes_mean[,28] <- rowSums(t(round(incHosp))[,(322+365+365+365):(322+365+365+365+365)])   #round(previcureq1)
  results$Outcomes_mean[,29] <- rowSums(t(round(incICU))[,(322+365+365+365):(322+365+365+365+365)])        #round(previcureq21)
  results$Outcomes_mean[,30] <- rowSums(t(round(incICUv))[,(322+365+365+365):(322+365+365+365+365)])  #round(previcureq31)
  results$Outcomes_mean[,31] <- rowSums(t(round(incDeath))[,(322+365+365+365):(322+365+365+365+365)])
  
  colnames(results$Outcomes_mean) <- c("Vac_pop","Asym_Y1","Sym_Y1","Hosp_Y1","ICU_Y1","ICUV_Y1","Death_Y1",
                                       "Asym_Y2","Sym_Y2","Hosp_Y2","ICU_Y2","ICUV_Y2","Death_Y2",
                                       "Asym_Y3","Sym_Y3","Hosp_Y3","ICU_Y3","ICUV_Y3","Death_Y3",
                                       "Asym_Y4","Sym_Y4","Hosp_Y4","ICU_Y4","ICUV_Y4","Death_Y4",
                                       "Asym_Y5","Sym_Y5","Hosp_Y5","ICU_Y5","ICUV_Y5","Death_Y5")
  rownames(results$Outcomes_mean) <- c("0-4 y.o.","5-9 y.o.","10-14 y.o.","15-19 y.o.","20-24 y.o.","25-29 y.o.","30-34 y.o.","35-39 y.o.","40-44 y.o.","45-49 y.o.","50-54 y.o.","55-59 y.o.","60-64 y.o.","65-69 y.o.","70-74 y.o.","75-79 y.o.","80-84 y.o.","85-89 y.o.","90-94 y.o.","95-99 y.o.",">=100 y.o.")
  
  
  #results$Outcomes_PSA <- array(NA, dim=c(4,3,1000)) #  Vac_pop + 6 outputs (vac_pop, asym, sym, hosp, icu, icuven, death) * 5 years (2020-2024) = 31 in total
  
  return(results)
}




########### 

###########    RUN BASELINE MODEL - start time for interventions is set to day 1e5, i.e. interventions are always off
Y<-c(initS,initE,initI,initR,initX,initH,initHC,initC,initCM,initSV,initEV,initIV,initCLV,initHV,initICUV,initVentV,initRV, initQS, initQE, initQI,initQR, initCL, initQC, initICU, initICUC, initVent, initVentC, initCMC) # initial conditions for the main solution vector

# INITIAL RESULTS AND PLOT
out<- ode(y = Y, times = times, func = covid, parms = parameters, method = 'euler',hini=0.05)
results.vaceff00 <- process_ode_outcome(out)
results.vaceff01 <- process_ode_outcome(out)
#results.vaceff00$Outcomes_mean

plot(x=times.date, y=results.vaceff00$daily_incidence,type = 'l',col='black',ylim = c(0,500), xlab="Year", ylab="Daily Incidence")
lines(x=times.date,y=rowSums(results.vaceff01$asymp_cases),type = 'l',col='brown')
last(results.vaceff00$cum_cases)
results.vaceff00$total_deaths_end




plot(x=times.date, y=results.vaceff00$daily_incidence,type = 'l',col='black',ylim = c(0,1000), xlab="Year", ylab="Daily Incidence")
plot(x=times.date, y=results.vaceff00$daily_incidence,type = 'l',col='black',ylim = c(0,10000), xlab="Year", ylab="Daily Incidence")

lines(x=times.date,y=rowSums(results.vaceff01$symp_cases),type = 'l',col='blue')
lines(x=times.date,y=rowSums(results.vaceff00$hospital),type = 'l',col='orange')
lines(x=times.date,y=rowSums(results.vaceff00$icu),type = 'l',col='green')
lines(x=times.date,y=rowSums(results.vaceff00$icuvent),type = 'l',col='red')
lines(x=times.date,y=rowSums(results.vaceff00$death),type = 'l',col='black')

