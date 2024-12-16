##########################Loading R packages#############################
library(mrgsolve)
library(dplyr)
library(zoo)   
library(ggplot2)
library(patchwork)
library(FME)
library(Metrics)
library(tidyr)
library(magrittr)    
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(MESS) 
library(RColorBrewer)
library(tidyverse) # Easily Install and Load the 'Tidyverse'

##########################Rat-Gd-IV-PBPK Model#############################

# Define a mrgsolve-based PBPK model for rats with intravenous Gd administration

rat_gd_iv_PBPK_code <- '

$PARAM @annotated

// Cardiac Output and Blood Flow Parameters
QCC            : 14.0    : L/h/kg^0.75,                   Cardiac output (Brown, 1997, Table 22)
QLC            : 0.183   : %,                             Fraction blood flow to liver (Brown, 1997, Table 23)
QLuC           : 1.00    : %,                             Fraction blood flow to lungs (Brown, 1997, Table 23)
QKC            : 0.141   : %,                             Fraction blood flow to kidneys (Brown, 1997, Table 23)
QBrC           : 0.02    : %,                             Fraction blood flow to brain (Brown, 1997, Table 23)
QSC            : 0.0146  : %,                             Fraction blood flow to spleen (Davies and Morries, 1993)

// Tissue Volume Parameters
BW             : 0.3     : kg,                            Body weight of the rat (Brown, 1997）
VLC            : 0.035   : L/kg BW,                       Fraction liver tissue (Brown, 1997, Table 21)
VLuC           : 0.005   : L/kg BW,                       Fraction lung tissue (Brown, 1997, Table 21)
VKC            : 0.007   : L/kg BW,                       Fraction kidney tissue (Brown, 1997, Table 21)
VBrC           : 0.006   : L/kg BW,                       Fraction brain tissue (Brown, 1997, Table 21)
VSC            : 0.031   : L/kg BW,                       Fraction spleen tissue (Davies and Morries, 1993)
VPlasC         : 0.0312  : L/kg BW,                       Fraction plasma (Davies and Morries, 1993)

// Blood Volume Fractions in Organs and Tissues
BVL            : 0.21    : %,                             Fraction blood volume in liver (Brown, 1997, Table 30)
BVBr           : 0.03    : %,                             Fraction blood volume in brain (Brown, 1997, Table 30)
BVK            : 0.16    : %,                             Fraction blood volume in kidneys (Brown, 1997, Table 30)
BVS            : 0.22    : %,                             Fraction blood volume in spleen (Brown, 1997, Table 30)
BVLu           : 0.36    : %,                             Fraction blood volume in lungs (Brown, 1997, Table 30)
BVR            : 0.04    : %,                             Fraction blood volume in rest of the body

// Partition Coefficients (Tissue to Plasma)
PL             : 0.209   : unitless,                      Partition coefficient for liver (Li, 2016)
PK             : 0.209   : unitless,                      Partition coefficient for kidneys (Li, 2016)
PBr            : 0.209   : unitless,                      Partition coefficient for brain (Li, 2016)
PS             : 0.209   : unitless,                      Partition coefficient for spleen (Li, 2016)
PLu            : 0.209   : unitless,                      Partition coefficient for lungs (Li, 2016)
PR             : 0.209   : unitless,                      Partition coefficient for rest of the body (Li, 2016)

// Membrane-Limited Permeability Coefficient Constants
PALC           : 0.776   : mg/L,                          Permeability constant for liver (Li, 2016)
PABrC          : 0.000000675 : mg/L,                      Permeability constant for brain (Li, 2016)
PAKC           : 0.776   : mg/L,                          Permeability constant for kidneys (Li, 2016)
PASC           : 0.776   : mg/L,                          Permeability constant for spleen (Li, 2016)
PALuC          : 0.776   : mg/L,                          Permeability constant for lungs (Li, 2016)
PARC           : 0.0171  : mg/L,                          Permeability constant for rest of the body (Li, 2016)

// Endocytic Parameters
KLRESrelease   : 5.30e-19 : 1/h,                          Release rate constant for phagocytic cells in liver (Mathieu,2021)
KLRESabsorb    : 1.45     : 1/h,                          Uptake rate constant for phagocytic cells in liver (Mathieu,2021)

KSRESrelease   : 5.30e-19 : 1/h,                          Release rate constant for phagocytic cells in spleen (Mathieu,2021)
KSRESabsorb    : 0.518    : 1/h,                          Uptake rate constant for phagocytic cells in spleen (Mathieu,2021)

KKRESrelease   : 5.30e-19 :1/h,                           Release rate constant of phyagocytic cells in kidney (Mathieu,2021)
KKRESabsorb    : 1.45     :1/h,                           Uptake rate constant for phagocytic cells in kidney (Mathieu,2021)

KLuRESrelease  : 5.30e-19 :1/h,                           Release rate constant of phyagocytic cells in lung (Mathieu,2021)
KLuRESabsorb   : 1.45     :1/h,                           Uptake rate constant for phagocytic cells in lung (Mathieu,2021)

KBrRESrelease  : 5.30e-19 :1/h,                           Release rate constant of phyagocytic cells in Brain (Mathieu,2021)
KBrRESabsorb   : 1.45     :1/h,                           Uptake rate constant for phagocytic cells in Brain (Mathieu,2021)

KRRESrelease   : 5.30e-19 :1/h,                           Release rate constant of phyagocytic cells in Brain (Mathieu,2021)
KRRESabsorb    : 1.45     :1/h,                           Uptake rate constant for phagocytic cells in Brain (Mathieu,2021)

// Excretion Parameters
KbileC         : 0.00152  : L/h/kg^0.75,                  Bile clearance
KurineC        : 0.141    : L/h/kg^0.75,                  Urine clearance

$MAIN
double QRC     = 1 - (QLC + QKC  + QSC  + QBrC );                          //Fraction of blood flow to rest of body
double VRC     = 1 - (VLC + VLuC + VKC  + VSC + VBrC + VPlasC);                          //Tissue volume of rest of body

double QC      = QCC * pow(BW, 0.75);                                                    //L/h, Cardiac output (adjusted for plasma)
double QL      = QC * QLC;                                                           //L/h, Blood flow to liver
double QBr     = QC * QBrC;                                                         //L/h, Blood flow to brain
double QK      = QC * QKC;                                                           //L/h, Blood flow to kidney
double QS      = QC * QSC;                                                          //L/h, Blood flow to spleen
double QR      = QC * QRC;                                                           //L/h, Blood flow to the rest of body


//Tissue volumes
double VL      = BW * VLC;                                                           //L, Liver volume
double VBr     = BW * VBrC;                                                          //L, Brain volume
double VK      = BW * VKC;                                                           //L, Kidney volume
double VS      = BW * VSC;                                                          //L, spleen volume
double VLu     = BW * VLuC;                                                          //L, lung volume
double VR      = BW * VRC;                                                           //L, Rest of body
double VPlasma = BW * VPlasC;                                                    //L, Plasma

double VLb     = VL * BVL; 							                                            //Weight/volume of capillary blood in liver compartment
double VLt     = VL - VLb; 							                                            //Weight/volume of tissue in liver compartment
double VBrb    = VBr * BVBr; 						                                            //Weight/volume of capillary blood in brain compartment
double VBrt    = VBr - VBrb; 						                                            //Weight/volume of tissue in brain compartment
double VKb     = VK * BVK; 							                                            //Weight/volume of capillary blood in kidney compartment
double VKt     = VK - VKb; 							                                            //Weight/volume of tissue in kidney compartment
double VSb     = VS * BVS; 							                                          //Weight/volume of capillary blood in spleen compartment
double VSt     = VS - VSb; 							                                          //Weight/volume of tissue in spleen compartment
double VLub    = VLu * BVLu; 						                                            //Weight/volume of capillary blood in lung compartment
double VLut    = VLu - VLub; 						                                            //Weight/volume of tissue in lung compartment
double VRb     = VR * BVR; 							                                            //Weight/volume of capillary blood in rest of body compartment
double VRt     = VR - VRb; 							                                            //Weight/volume of tissue in rest of body compartment

//Permeability coefficient-surface area cross-product (mg/h)
double PAL     = PALC  * QL; 						                                            //mg/h, Liver
double PABr    = PABrC * QBr; 						                                          //mg/h, Brain
double PAK     = PAKC  * QK; 						                                            //mg/h, Kidneys
double PAS     = PASC  * QS; 						                                          //mg/h, Spleen
double PALu    = PALuC * QC;  						                                          //mg/h, Lungs
double PAR     = PARC  * QR; 						                                            //mg/h, Rest of body

//Endocytosis rate (1/h)
double KLRESUP =KLRESabsorb;
double KSRESUP = KSRESabsorb;
double KKRESUP = KKRESabsorb;
double KLuRESUP =KLuRESabsorb;
double KRRESUP  = KRRESabsorb;
double KBrRESUP  = KBrRESabsorb;

double Kurine  = KurineC * pow(BW, 0.75);
double Kbile   = KbileC * pow(BW, 0.75);   

$CMT AA AV ALub ALut ALuRES ABrb ABrt  ARb ARt ARRES Abile ABrRES AKb AKt AKRES Aurine ASb ASt ASRES ALb ALt ALRES ADOSE

$ODE

double APlasma  = AA + AV;
double ALung    = ALub+ALut+ALuRES;
double ALungt   = ALut+ALuRES;
double Arestall = ARb+ARt+ARRES;
double Aresttissue = ARt+ARRES;
double AKidney  = AKb+AKt+AKRES;
double AKidneyt = AKt+AKRES;
double ABrain   = ABrb + ABrt+ABrRES;
double ABraint   = ABrt+ABrRES;
double ASpleen  = ASb+ASt+ASRES;
double ASpleent = ASt+ASRES;
double ALiver   = ALb+ALt+ALRES;
double ALivert  = ALt+ALRES;



double CA       = AA/(VPlasma * 0.2);
double CV       = AV / (VPlasma * 0.8);
double CPlasma  = APlasma/VPlasma;
double CVLu     = ALub / VLub;
double CLut     = ALut / VLut;
double CLung    = (ALub + ALut + ALuRES)/VLu;
double CLungt   = (ALut+ALuRES)/VLut;
double CVBr     = ABrb/VBrb;
double CBrt     = (ABrt)/VBrt;
double CBrain   = (ABrain+ABrRES+ ABrb)/ VBr;
double CBraint  = (ABrt+ABrRES)/VBrt;
double CVR      = ARb/VRb;
double CRt      = ARt/VRt;
double Crestall = (ARb+ARt+ARRES)/VR;
double Cresttissue = (ARt+ARRES)/VRt;
double CVK      = AKb/VKb;
double CKt      = AKt/VKt;
double CKidney  = (AKb+AKt+AKRES)/VK;
double CKidneyt = (AKt+AKRES)/VKt;
double CVS      = ASb/VSb;
double CSt      = ASt/VSt;
double CSpleen  = (ASb+ASt+ASRES)/VS;
double CSpleent = (ASt+ASRES)/VSt;
double CVL      = ALb/VLb;
double CLt      = ALt/VLt;
double CLiver = (ALb+ALt+ALRES)/VL;
double CLivert = (ALt+ALRES)/VLt;

double Rbile    = Kbile*ALt  ; 
double RA      = QC * CVLu - QC * CA; 
double RV      = QL * CVL + QBr  *CVBr + QK * CVK +  QR * CVR - QC * CV;
double RLub    =  QC * (CV - CVLu) - PALu * CVLu + (PALu * CLut)/ PLu; 
double RLut    = PALu * CVLu - (PALu * CLut)/ PLu - KLuRESUP * ALut + KLuRESrelease * ALuRES;
double RLuRES  = KLuRESUP * ALut - KLuRESrelease * ALuRES;
double RBrb    = QBr *(CA - CVBr) - PABr * CVBr + (PABr * CBrt)/ PBr;
double RBrt    = PABr * CVBr - (PABr * CBrt)/ PBr+ KBrRESrelease*ABrRES;
double RBrRES  = KBrRESUP*ABrt-KBrRESrelease*ABrRES;
double RRb      = QR*(CA-CVR) - PAR*CVR + (PAR*CRt)/PR;
double RRt      = PAR*CVR - (PAR*CRt)/PR - KRRESUP*ARt + KRRESrelease*ARRES;
double RRRES    = KRRESUP*ARt-KRRESrelease*ARRES;
double Rurine   = Kurine*AKt;
double RKb      = QK*(CA-CVK) - PAK*CVK + (PAK*CKt)/PK;
double RKt      = PAK*CVK - (PAK*CKt)/PK - KKRESUP*AKt + KKRESrelease*AKRES - Rurine; 
double RKRES    = KKRESUP*AKt-KKRESrelease*AKRES;
double RSb      = QS*(CA-CVS) - PAS*CVS + (PAS*CSt)/PS; 
double RSt      = PAS*CVS - (PAS*CSt)/PS - KSRESUP*ASt + KSRESrelease*ASRES;
double RSRES    = KSRESUP*ASt-KSRESrelease*ASRES;
double RLb      = QL*(CA-CVL) + QS*CVS - PAL*CVL + (PAL*CLt)/PL  ;
double RLt      = PAL*CVL - (PAL*CLt)/PL - Rbile- KLRESUP*ALb + KLRESrelease*ALRES;
double RLRES    = KLRESUP*ALb-KLRESrelease*ALRES;



//Concentration of the chemical in blood compartment
//CA: Arterial blood concentration (mg/L or ug/ml)
dxdt_AA     = RA;
dxdt_AV     = RV;
dxdt_ALub   = RLub;
dxdt_ALut   = RLut;
dxdt_ALuRES = RLuRES;
dxdt_ABrb   = RBrb;
dxdt_ABrt   = RBrt;
dxdt_ABrRES = RBrRES;
dxdt_ARb    = RRb;
dxdt_ARt    = RRt;
dxdt_ARRES  = RRRES;
dxdt_AKb    = RKb;
dxdt_AKt    = RKt;
dxdt_AKRES  = RKRES;
dxdt_Aurine = Rurine;
dxdt_ASb    = RSb;
dxdt_ASt    = RSt;
dxdt_ASRES  = RSRES;
dxdt_Abile  = Rbile;
dxdt_ALb    = RLb;
dxdt_ALt    = RLt;
dxdt_ALRES  = RLRES;
dxdt_ADOSE  = 0;

// {Mass balance equations}
double Tmass = APlasma + ALiver + ABrain + AKidney + ALung + Arestall  + ASpleen  +Abile+ Aurine  ;
double Bal   =Tmass;



$TABLE
capture Lung     = CLung;
capture Liver    = CLiver;
capture Kidney   = CKidney;
capture Spleen   = CSpleen;
capture Plasma   = CPlasma;
capture Brain   = CBrain;
capture BAL      = Bal;

'

# Compile the PBPK model
mod <- mcode("rat_gd_iv_PBPK_code", rat_gd_iv_PBPK_code)

##########################Rat PBPK Model:Intravenous injection#############################

# Define the function for intravenous injection simulation
pred.A <- function(Gpars, n_bootstrap = 1000, conf_level = 0.95) {
  
  # Transform parameters from the log scale to the normal scale
  Gpars <- exp(Gpars)  # Return the parameters in normal scale
  
  # Define the exposure scenario for intravenous injection
  GBW <- 0.3              # Body weight of the rat (kg)
  tinterval <- 0          # Time interval between doses (hours)
  GTDOSE <- 1             # Number of doses
  DOSE <- 26.35           # Dose in mg/kg
  GDOSEiv <- DOSE * GBW   # Calculate the intravenous dose
  
  # Define the event for intravenous injection
  Gex.oral.B <- ev(
    ID = 1,         # Subject ID
    time = 0,       # Start time
    amt = GDOSEiv,  # Amount of dose
    ii = tinterval, # Time interval between doses
    addl = GTDOSE - 1,  # Additional doses (none in this case)
    cmt = "AV",     # Compartment for the dose (assumed to be 'AV')
    replicate = FALSE  # No replication of the event
  )
  
  # Set up the time grid for simulation output (0 to 2 hours, with steps of 0.005 hours)
  Gtsamp <- tgrid(0, tinterval * (GTDOSE - 1) + 2, 0.005)
  
  # Define the function to run a single simulation
  run_simulation <- function(pars) {
    mod %>%
      param(pars) %>%
      update(atol = 1E-6, maxsteps = 500000) %>%
      mrgsim_d(data = Gex.oral.B, tgrid = Gtsamp) %>%
      as.data.frame() %>%
      filter(time > 0)  # Filter out data where time <= 0
  }
  
  # Run the simulation with the original parameters
  Gout.B <- run_simulation(Gpars)
  
  # Perform bootstrapping to calculate confidence intervals
  bootstrap_results <- replicate(n_bootstrap, {
    sampled_pars <- Gpars * exp(rnorm(length(Gpars), 0, 0.1))  # Add small perturbations to the parameters
    run_simulation(sampled_pars)
  }, simplify = FALSE)
  
  # Combine bootstrap results into a data frame
  combined_results <- do.call(rbind, bootstrap_results)
  
  # Calculate confidence intervals for each time point
  conf_intervals <- combined_results %>%
    group_by(time) %>%
    summarise(across(c(Plasma, Liver, Kidney, Spleen, Lung, Brain, BAL),
                     list(mean = ~mean(.),
                          lower = ~quantile(., (1 - conf_level) / 2),
                          upper = ~quantile(., 1 - (1 - conf_level) / 2))))
  
  # Extract relevant columns from the output data
  Goutdf.B <- Gout.B %>%
    select(time, Plasma, Liver, Kidney, Spleen, Lung, Brain, BAL)
  
  # Return both the simulation data and the confidence intervals
  return(list("Goutdf.B" = Goutdf.B, "conf_intervals" = conf_intervals))
}

# Run the model with the initial parameter set
result <- pred.A(Gtheta.int)

# Write completed data to CSV,
write.csv(result, "Rat-iv-simulated-data.csv", row.names = FALSE)
##########################Find Sensitive Parameters#############################

# Define the cost function based on the observed data and model predictions
Gcost <- function(pars) { 
  
  # Run the model with the given parameter set
  outdf <- pred.A(Gpars = pars)
  
  # Calculate the cost (error) between model predictions and observed plasma concentrations
  cost <- modCost(model = outdf$Goutdf.B, obs = OBS.CPlas, x = "time", weight = "mean")
  
  # Add the cost for kidney, liver, and spleen concentrations
  cost <- modCost(model = outdf$Goutdf.B, obs = OBS.CK, x = "time", weight = "mean", cost = cost)
  cost <- modCost(model = outdf$Goutdf.B, obs = OBS.CL, x = "time", weight = "mean", cost = cost)
  cost <- modCost(model = outdf$Goutdf.B, obs = OBS.CS, x = "time", weight = "mean", cost = cost)
  
  # Return the total cost
  return(cost)
}

# Initialize parameters (log-scale)
Gtheta.int <- log(c(
  PL = 0.209,
  PK = 0.209, 
  PBr = 0.209,   
  PS = 0.209,                       
  PLu = 0.209,             
  PR = 0.209,                     
  PALC = 0.776,                     
  PABrC = 0.000000675,                     
  PAKC = 0.776,                     
  PASC = 0.776,                     
  PALuC = 0.776,                      
  PARC = 0.171,
  KLRESrelease = 5.30 * 10^-19,                       
  KLRESabsorb = 1.45,                  
  KSRESrelease = 5.30 * 10^-19,    
  KSRESabsorb = 0.518,     
  KKRESrelease = 5.30 * 10^-19,                   
  KKRESabsorb = 1.45, 
  KBrRESrelease = 5.30 * 10^-19,                   
  KBrRESabsorb = 1.45,   
  KLuRESrelease = 5.30 * 10^-19,                              
  KLuRESabsorb = 1.45,
  KRRESrelease = 5.30 * 10^-19,                
  KRRESabsorb = 1.45,     
  KurineC = 0.141,  
  Kbile = 0.00152
))

# Perform local sensitivity analysis
SnsPlasma <- sensFun(func = Gcost, parms = Gtheta.int, varscale = 1)

# Summarize and plot the sensitivity results
Sen <- summary(SnsPlasma)
plot(Sen)

## Selected Parameters to change

# Initialize the parameter set (log-scale) for the next step
Gtheta.int <- log(c(
  #PL = 0.209,
  #PK = 0.209,
  #PBr = 0.209,
  #PS = 0.209,                     
  #PLu = 0.209,             
  #PR = 0.209,             
  #PALC = 0.776,             
  #PABrC = 0.000000675,             
  #PAKC = 0.776,             
  #PASC = 0.776,             
  #PALuC = 0.776,             
  #PARC = 0.171,             
  #KLRESrelease = 5.30 * 10^-19,             
  #KLRESabsorb = 1.45,             
  #KSRESrelease = 5.30 * 10^-19,             
  KSRESabsorb = 0.518,  # Only the spleen absorption rate is selected
  #KKRESrelease = 5.30 * 10^-19,             
  KKRESabsorb = 1.45,  # Kidney absorption rate
  #KLuRESrelease = 5.30 * 10^-19,             
  #KLuRESabsorb = 1.45,             
  #KBrRESrelease = 5.30 * 10^-19,             
  #KBrRESabsorb = 1.45,             
  #KRRESrelease = 5.30 * 10^-19,             
  KurineC = 0.141,  # Urine clearance rate
  #Kbile = 0.00152,             
  KRRESabsorb = 1.45  # Absorption rate for the rest of the body
))

# Perform parameter fitting
# Generate upper and lower bounds (bounds are 20x and 1/20x of initial values)
lower_limits <- Gtheta.int - log(20)
upper_limits <- Gtheta.int + log(20)

# Perform the fitting process using the Nelder-Mead method
GFit <- modFit(f = Gcost, p = Gtheta.int, method = "Nelder-Mead", lower = lower_limits, upper = upper_limits)

# Extract the new parameters from the fitting result (returning to normal scale)
NewPar <- exp(GFit$par)

# Print the result to verify
print(NewPar)

##########################Sensitivity Analysis#############################

# Define initial parameter values for the rat model
rat.theta.G <- log(c(
  PL = 0.209,
  PK = 0.209,
  PBr = 0.209,
  PS = 0.209,                   
  PLu = 0.209,          
  PR = 0.209,            
  PALC = 0.776,                  
  PABrC = 0.000000675,                   
  PAKC = 0.776,             
  PASC = 15.519976041206,                 
  PALuC = 0.776,             
  PARC = 0.171,
  KLRESrelease = 5.30e-19,                     
  KLRESabsorb = 14.058133222604,           
  KSRESrelease = 5.30e-19,  
  KSRESabsorb = 10.359965835688, 
  KKRESrelease = 5.30e-19,                   
  KKRESabsorb = 0.072500589917,              
  KLuRESrelease = 5.30e-19,                           
  KLuRESabsorb = 1.45,
  KBrRESrelease = 5.30e-19,                       
  KBrRESabsorb = 1.45,
  KRRESrelease = 5.30e-19,         
  KurineC = 0.007050035498,
  Kbile = 0.00152,
  KRRESabsorb = 0.072500093305 
))

# Define the model prediction function

pred.rat <- function(Gpars) { 
  # Convert parameters from log domain to normal
  Gpars <- exp(Gpars)
  
  # Define exposure scenario
  GBW = 0.3
  tinterval = 0
  GTDOSE = 1
  DOSE = 26.35
  GDOSEiv = DOSE * GBW
  
  Gex.oral.B <- ev(ID = 1,        
                   time = 0,      
                   amt = GDOSEiv,
                   ii = tinterval,  
                   addl = GTDOSE - 1, 
                   cmt = "AV", 
                   replicate = FALSE) 
  
  Gtsamp = tgrid(0, tinterval * (GTDOSE - 1) + 2, 0.01)
  
  # Simulate exposure scenario
  Gout.B <- mod %>% 
    param(Gpars) %>% 
    update(atol = 1E-6, maxsteps = 500000) %>% 
    mrgsim_d(data = Gex.oral.B, tgrid = Gtsamp)
  
  # Extract relevant columns and calculate AUC
  Goutdf.B <- cbind.data.frame(
    time = Gout.B$time, 
    CPlas = Gout.B$Plasma,     
    CL = Gout.B$Liver,     
    CK = Gout.B$Kidney,    
    CS = Gout.B$Spleen,    
    CLu = Gout.B$Lung, 
    CBr = Gout.B$Brain,
    Bal = Gout.B$BAL
  )
  
  # Calculate AUC for each organ
  Goutdf.B <- Goutdf.B %>% 
    mutate(
      AUC_CPlas = auc(time, CPlas),
      AUC_Liver = auc(time, CL),
      AUC_Kidney = auc(time, CK),
      AUC_Spleen = auc(time, CS),
      AUC_Lung = auc(time, CLu),
      AUC_Brain = auc(time, CBr)
    )
  
  # Filter out values with time greater than 0
  Goutdf.B <- Goutdf.B %>% filter(time > 0)
  
  return(list("G" = Goutdf.B))
}

# Call the prediction function
result <- pred.rat(GFit$par)

# Check the results
print(result)

par.G <- GFit$par
rat.theta.G[names(par.G)] <- as.numeric(par.G)

# Sensitivity analysis function
NSC_func <- function(Gpars, Pred) {
  nG <- length(Gpars)
  NSC_GCA <- matrix(NA, nrow = length(Gpars), ncol = 6)
  
  for (i in 1:nG) {
    Gpars.new <- Gpars %>% replace(i, log(exp(Gpars[i]) * 1.01))
    Mnew.G <- Pred(Gpars.new)
    M.G <- Pred(Gpars)
    delta.Gpars <- exp(Gpars[i]) / (exp(Gpars[i]) * 0.01)
    
    # Extract AUC values
    AUC.GPlas.new <- as.numeric(Mnew.G$G$AUC_CPlas[nrow(Mnew.G$G)])
    AUC.GPlas.ori <- as.numeric(M.G$G$AUC_CPlas[nrow(M.G$G)])
    AUC.GLiver.new <- as.numeric(Mnew.G$G$AUC_Liver[nrow(Mnew.G$G)])
    AUC.GLiver.ori <- as.numeric(M.G$G$AUC_Liver[nrow(M.G$G)])
    AUC.GKidney.new <- as.numeric(Mnew.G$G$AUC_Kidney[nrow(Mnew.G$G)])
    AUC.GKidney.ori <- as.numeric(M.G$G$AUC_Kidney[nrow(M.G$G)])
    AUC.GSpleen.new <- as.numeric(Mnew.G$G$AUC_Spleen[nrow(Mnew.G$G)])
    AUC.GSpleen.ori <- as.numeric(M.G$G$AUC_Spleen[nrow(M.G$G)])
    AUC.GLung.new <- as.numeric(Mnew.G$G$AUC_Lung[nrow(Mnew.G$G)])
    AUC.GLung.ori <- as.numeric(M.G$G$AUC_Lung[nrow(M.G$G)])
    AUC.GBrain.new <- as.numeric(Mnew.G$G$AUC_Brain[nrow(Mnew.G$G)])
    AUC.GBrain.ori <- as.numeric(M.G$G$AUC_Brain[nrow(M.G$G)])
    
    # Calculate delta
    delta.AUC.GPlas <- AUC.GPlas.new - AUC.GPlas.ori
    delta.AUC.GLiver <- AUC.GLiver.new - AUC.GLiver.ori
    delta.AUC.GKidney <- AUC.GKidney.new - AUC.GKidney.ori
    delta.AUC.GSpleen <- AUC.GSpleen.new - AUC.GSpleen.ori
    delta.AUC.GLung <- AUC.GLung.new - AUC.GLung.ori
    delta.AUC.GBrain <- AUC.GBrain.new - AUC.GBrain.ori
    
    NSC_GCA[i, 1] <- (delta.AUC.GPlas / AUC.GPlas.ori) * delta.Gpars
    NSC_GCA[i, 2] <- (delta.AUC.GLiver / AUC.GLiver.ori) * delta.Gpars
    NSC_GCA[i, 3] <- (delta.AUC.GKidney / AUC.GKidney.ori) * delta.Gpars
    NSC_GCA[i, 4] <- (delta.AUC.GSpleen / AUC.GSpleen.ori) * delta.Gpars
    NSC_GCA[i, 5] <- (delta.AUC.GLung / AUC.GLung.ori) * delta.Gpars
    NSC_GCA[i, 6] <- (delta.AUC.GBrain / AUC.GBrain.ori) * delta.Gpars
  }
  
  return(NSC_GCA)
}

# Perform sensitivity analysis
A <- NSC_func(rat.theta.G, pred.rat)
A <- as.data.frame(A)

# Create ID column
ID <- 1:nrow(A)

# Extract parameter names
parameter <- gsub("SA_", "", names(rat.theta.G))

# Add new columns to the dataframe
A <- cbind(ID = ID, parameter = parameter, A)

# Rename columns
colnames(A)[3:ncol(A)] <- c("NSC_CPlas", "NSC_CLiver", "NSC_CKidney", "NSC_CSpleen", "NSC_CLung", "NSC_CBrain")
NSC_GCA_M <- data.frame(A)

# Output results to CSV
write.csv(NSC_GCA_M, "Rat-sensitive analysis-data", row.names = FALSE)


##########################rat-oral-PBPK Model##########################

# Define a mrgsolve-based PBPK model for oral dosing in rats
rat-oral-PBPK.code <- '

$PARAM @annotated

// Cardiac Output and Blood Flow
QCC            : 14.0    :L/h/kg^0.75,                    Cardiac output, (Brown, 1997, Table 22)
QLC            : 0.183   :%                       Fraction of blood flow to the liver, (Brown, 1997, Table 23)
QLuC           : 1       :%,                       Fraction of blood flow to the lung, (Brown, 1997, Table 23)
QKC            : 0.141   :%,                       Fraction of blood flow to the kidney, (Brown, 1997, Table 23)
QBrC           : 0.02    :%,                       Fraction of blood flow to the brain, (Brown, 1997, Table 23)
QSC            : 0.0146  :%,                       Fraction of blood flow to the spleen (Davies and Morries, 1993)*

// Tissue Volume
BW             : 0.3     :kg,                       Body weight
VLC            : 0.035   :L/kg BW,                  Fraction of liver tissue, (Brown, 1997, Table 21)
VLuC           : 0.005   :L/kg BW,                  Fraction of lung tissue, (Brown, 1997, Table 21)
VKC            : 0.007   :L/kg BW,                  Fraction of kidney tissue, (Brown, 1997, Table 21)
VBrC           : 0.006   :L/kg BW,                  Fraction of brain tissue, (Brown, 1997, Table 21)
VSC            : 0.031   :L/kg BW                   Fraction of spleen tissue (Davies and Morries, 1993)*
VPlasC         : 0.0312  :L/kg BW,                  Fraction of plasma, (Davies and Morries, 1993)

// Blood Volume Fraction in Organs and Tissues
BVL            : 0.21    :%,                        Liver blood volume, (Brown, 1997, Table 30)
BVBr           : 0.03    :%,                        Brain blood volume, (Brown, 1997, Table 30)
BVK            : 0.16    :%,                        Kidney blood volume, (Brown, 1997, Table 30)
BVS            : 0.22    :%,                        Spleen blood volume, (Brown, 1997, Table 30)
BVLu           : 0.36    :%,                        Lung blood volume, (Brown, 1997, Table 30)
BVR            : 0.04    :%,                        Rest of body blood volume

// Partition Coefficients (Tissue to Plasma)
PL             : 0.209   :unitless,                 Liver partition coefficient, (Li, 2016)
PK             : 0.209   :unitless,                 Kidney partition coefficient, (Li, 2016)
PBr            : 0.209   :unitless,                 Brain partition coefficient, (Li, 2016)
PS             : 0.209   :unitless,                 Spleen partition coefficient, (Li, 2016)
PLu            : 0.209   :unitless,                 Lung partition coefficient, (Li, 2016)
PR             : 0.209   :unitless,                 Rest of body partition coefficient, (Li, 2016)

// Membrane-limited Permeability Coefficients
PALC           : 0.776   :mg/L,                     Liver permeability coefficient, (Li, 2016)
PABrC          : 0.000000675 :mg/L,                Brain permeability coefficient, (Li, 2016)
PAKC           : 0.776   :mg/L,                     Kidney permeability coefficient, (Li, 2016)
PASC           : 15.5199 :mg/L,                     Spleen permeability coefficient, (Li, 2016)
PALuC          : 0.776   :mg/L,                     Lung permeability coefficient, (Li, 2016)
PARC           : 0.0171  :mg/L,                     Rest of body permeability coefficient, (Li, 2016)

// Endocytic Parameters for Phagocytic Cells
KLRESrelease   : 0.00000000000000000053 :1/h,       Liver release rate constant for phagocytic cells
KLRESabsorb    : 14.0581 :1/h,                      Liver uptake rate constant for phagocytic cells
KSRESrelease   : 0.00000000000000000053 :1/h,       Spleen release rate constant for phagocytic cells
KSRESabsorb    : 10.3599 :1/h,                      Spleen uptake rate constant for phagocytic cells
KKRESrelease   : 0.00000000000000000053 :1/h,       Kidney release rate constant for phagocytic cells
KKRESabsorb    : 0.0725  :1/h,                      Kidney uptake rate constant for phagocytic cells
KLuRESrelease  : 0.00000000000000000053 :1/h,       Lung release rate constant for phagocytic cells
KLuRESabsorb   : 1.45    :1/h,                      Lung uptake rate constant for phagocytic cells
KBrRESrelease  : 0.00000000000000000053 :1/h,       Brain release rate constant for phagocytic cells
KBrRESabsorb   : 1.45    :1/h,                      Brain uptake rate constant for phagocytic cells
KRRESrelease   : 0.00000000000000000053 :1/h,       Rest of body release rate constant for phagocytic cells
KRRESabsorb    : 0.0725  :1/h,                      Rest of body uptake rate constant for phagocytic cells

// Excretion Parameters
KbileC         : 0.00152 :L/hr/kg^0.75,             Bile clearance
KurineC        : 0.00705 :L/hr/kg^0.75,             Urine clearance

$MAIN
double QRC     = 1 - (QLC + QKC  + QSC  + QBrC );                          //Fraction of blood flow to rest of body
double VRC     = 1 - (VLC + VLuC + VKC  + VSC + VBrC + VPlasC);                          //Tissue volume of rest of body

double QC      = QCC * pow(BW, 0.75);                                                    //L/h, Cardiac output (adjusted for plasma)
double QL      = QC * QLC;                                                           //L/h, Blood flow to liver
double QBr     = QC * QBrC;                                                         //L/h, Blood flow to brain
double QK      = QC * QKC;                                                           //L/h, Blood flow to kidney
double QS      = QC * QSC;                                                          //L/h, Blood flow to spleen
double QR      = QC * QRC;                                                           //L/h, Blood flow to the rest of body


//Tissue volumes
double VL      = BW * VLC;                                                           //L, Liver volume
double VBr     = BW * VBrC;                                                          //L, Brain volume
double VK      = BW * VKC;                                                           //L, Kidney volume
double VS      = BW * VSC;                                                          //L, spleen volume
double VLu     = BW * VLuC;                                                          //L, lung volume
double VR      = BW * VRC;                                                           //L, Rest of body
double VPlasma = BW * VPlasC;                                                    //L, Plasma

double VLb     = VL * BVL; 							                                            //Weight/volume of capillary blood in liver compartment
double VLt     = VL - VLb; 							                                            //Weight/volume of tissue in liver compartment
double VBrb    = VBr * BVBr; 						                                            //Weight/volume of capillary blood in brain compartment
double VBrt    = VBr - VBrb; 						                                            //Weight/volume of tissue in brain compartment
double VKb     = VK * BVK; 							                                            //Weight/volume of capillary blood in kidney compartment
double VKt     = VK - VKb; 							                                            //Weight/volume of tissue in kidney compartment
double VSb     = VS * BVS; 							                                          //Weight/volume of capillary blood in spleen compartment
double VSt     = VS - VSb; 							                                          //Weight/volume of tissue in spleen compartment
double VLub    = VLu * BVLu; 						                                            //Weight/volume of capillary blood in lung compartment
double VLut    = VLu - VLub; 						                                            //Weight/volume of tissue in lung compartment
double VRb     = VR * BVR; 							                                            //Weight/volume of capillary blood in rest of body compartment
double VRt     = VR - VRb; 							                                            //Weight/volume of tissue in rest of body compartment

//Permeability coefficient-surface area cross-product (mg/h)
double PAL     = PALC  * QL; 						                                            //mg/h, Liver
double PABr    = PABrC * QBr; 						                                          //mg/h, Brain
double PAK     = PAKC  * QK; 						                                            //mg/h, Kidneys
double PAS     = PASC  * QS; 						                                          //mg/h, Spleen
double PALu    = PALuC * QC;  						                                          //mg/h, Lungs
double PAR     = PARC  * QR; 						                                            //mg/h, Rest of body

//Endocytosis rate (1/h)
double KLRESUP =KLRESabsorb;
double KSRESUP = KSRESabsorb;
double KKRESUP = KKRESabsorb;
double KLuRESUP =KLuRESabsorb;
double KRRESUP  = KRRESabsorb;
double KBrRESUP  = KBrRESabsorb;

double Kurine  = KurineC * pow(BW, 0.75);
double Kbile   = KbileC * pow(BW, 0.75);   

$CMT AA AV ALub ALut ALuRES ABrb ABrt  ARb ARt ARRES Abile ABrRES AKb AKt AKRES Aurine ASb ASt ASRES ALb ALt ALRES ADOSE

$ODE

double APlasma  = AA + AV;
double ALung    = ALub+ALut+ALuRES;
double ALungt   = ALut+ALuRES;
double Arestall = ARb+ARt+ARRES;
double Aresttissue = ARt+ARRES;
double AKidney  = AKb+AKt+AKRES;
double AKidneyt = AKt+AKRES;
double ABrain   = ABrb + ABrt+ABrRES;
double ABraint   = ABrt+ABrRES;
double ASpleen  = ASb+ASt+ASRES;
double ASpleent = ASt+ASRES;
double ALiver   = ALb+ALt+ALRES;
double ALivert  = ALt+ALRES;



double CA       = AA/(VPlasma * 0.2);
double CV       = AV / (VPlasma * 0.8);
double CPlasma  = APlasma/VPlasma;
double CVLu     = ALub / VLub;
double CLut     = ALut / VLut;
double CLung    = (ALub + ALut + ALuRES)/VLu;
double CLungt   = (ALut+ALuRES)/VLut;
double CVBr     = ABrb/VBrb;
double CBrt     = (ABrt)/VBrt;
double CBrain   = (ABrain+ABrRES+ ABrb)/ VBr;
double CBraint  = (ABrt+ABrRES)/VBrt;
double CVR      = ARb/VRb;
double CRt      = ARt/VRt;
double Crestall = (ARb+ARt+ARRES)/VR;
double Cresttissue = (ARt+ARRES)/VRt;
double CVK      = AKb/VKb;
double CKt      = AKt/VKt;
double CKidney  = (AKb+AKt+AKRES)/VK;
double CKidneyt = (AKt+AKRES)/VKt;
double CVS      = ASb/VSb;
double CSt      = ASt/VSt;
double CSpleen  = (ASb+ASt+ASRES)/VS;
double CSpleent = (ASt+ASRES)/VSt;
double CVL      = ALb/VLb;
double CLt      = ALt/VLt;
double CLiver = (ALb+ALt+ALRES)/VL;
double CLivert = (ALt+ALRES)/VLt;

double Rbile    = Kbile*ALt  ; 
double RA      = QC * CVLu - QC * CA; 
double RV      = QL * CVL + QBr  *CVBr + QK * CVK +  QR * CVR - QC * CV;
double RLub    =  QC * (CV - CVLu) - PALu * CVLu + (PALu * CLut)/ PLu; 
double RLut    = PALu * CVLu - (PALu * CLut)/ PLu - KLuRESUP * ALut + KLuRESrelease * ALuRES;
double RLuRES  = KLuRESUP * ALut - KLuRESrelease * ALuRES;
double RBrb    = QBr *(CA - CVBr) - PABr * CVBr + (PABr * CBrt)/ PBr;
double RBrt    = PABr * CVBr - (PABr * CBrt)/ PBr+ KBrRESrelease*ABrRES;
double RBrRES  = KBrRESUP*ABrt-KBrRESrelease*ABrRES;
double RRb      = QR*(CA-CVR) - PAR*CVR + (PAR*CRt)/PR;
double RRt      = PAR*CVR - (PAR*CRt)/PR - KRRESUP*ARt + KRRESrelease*ARRES;
double RRRES    = KRRESUP*ARt-KRRESrelease*ARRES;
double Rurine   = Kurine*AKt;
double RKb      = QK*(CA-CVK) - PAK*CVK + (PAK*CKt)/PK;
double RKt      = PAK*CVK - (PAK*CKt)/PK - KKRESUP*AKt + KKRESrelease*AKRES - Rurine; 
double RKRES    = KKRESUP*AKt-KKRESrelease*AKRES;
double RSb      = QS*(CA-CVS) - PAS*CVS + (PAS*CSt)/PS; 
double RSt      = PAS*CVS - (PAS*CSt)/PS - KSRESUP*ASt + KSRESrelease*ASRES;
double RSRES    = KSRESUP*ASt-KSRESrelease*ASRES;
double RLb      = QL*(CA-CVL) + QS*CVS - PAL*CVL + (PAL*CLt)/PL  ;
double RLt      = PAL*CVL - (PAL*CLt)/PL - Rbile- KLRESUP*ALb + KLRESrelease*ALRES;
double RLRES    = KLRESUP*ALb-KLRESrelease*ALRES;



//Concentration of the chemical in blood compartment
//CA: Arterial blood concentration (mg/L or ug/ml)
dxdt_AA     = RA;
dxdt_AV     = RV;
dxdt_ALub   = RLub;
dxdt_ALut   = RLut;
dxdt_ALuRES = RLuRES;
dxdt_ABrb   = RBrb;
dxdt_ABrt   = RBrt;
dxdt_ABrRES = RBrRES;
dxdt_ARb    = RRb;
dxdt_ARt    = RRt;
dxdt_ARRES  = RRRES;
dxdt_AKb    = RKb;
dxdt_AKt    = RKt;
dxdt_AKRES  = RKRES;
dxdt_Aurine = Rurine;
dxdt_ASb    = RSb;
dxdt_ASt    = RSt;
dxdt_ASRES  = RSRES;
dxdt_Abile  = Rbile;
dxdt_ALb    = RLb;
dxdt_ALt    = RLt;
dxdt_ALRES  = RLRES;
dxdt_ADOSE  = 0;

// {Mass balance equations}
double Tmass = APlasma + ALiver + ABrain + AKidney + ALung + Arestall  + ASpleen  +Abile+ Aurine  ;
double Bal   =Tmass;



$TABLE
capture Lung     = CLung;
capture Liver    = CLiver;
capture Kidney   = CKidney;
capture Spleen   = CSpleen;
capture Plasma   = CPlasma;
capture Brain   = CBrain;
capture BAL      = Bal;

'

# Load the model and set initial values
mod <- mcode("rat-oral-PBPK.code", rat-oral-PBPK.code)

##########################Rat Oral Exposure Setup##########################

# Initialize parameters
Gtheta.int <- log(c(
  PL = 0.209,
  PK = 0.209, 
  PBr = 0.209,  
  PS = 0.209,                       
  PLu = 0.209,            
  PR = 0.209,                    
  PALC = 0.776,                     
  PABrC = 0.000000675,                     
  PAKC = 0.776,                     
  PASC = 15.519935855126,                     
  PALuC = 0.776,                      
  PARC = 0.171,
  KLRESrelease = 5.30 * 10^-19,                       
  KLRESabsorb = 14.058133222604,                  
  KSRESrelease = 5.30 * 10^-19,    
  KSRESabsorb = 10.359965835688,     
  KKRESrelease = 5.30 * 10^-19,                   
  KKRESabsorb = 0.072500589917, 
  KBrRESrelease = 5.30 * 10^-19,                   
  KBrRESabsorb = 1.45,   
  KLuRESrelease = 5.30 * 10^-19,                             
  KLuRESabsorb = 1.45,
  KRRESrelease = 5.30 * 10^-19,                
  KRRESabsorb = 0.072500093305,     
  KurineC = 0.007050035498, 
  Kbile = 0.00152
))

# Define a function for predictions
pred.A <- function(Gpars, doses, n_bootstrap = 1000, conf_level = 0.95) {
  # Gpars: input parameters, doses: list of doses, n_bootstrap: number of bootstrap iterations, conf_level: confidence level
  
  # Convert parameters from the log domain to the normal domain
  Gpars <- exp(Gpars)  # Return a list of gestational model parameters from log scale
  
  # Parameters for gestational exposure scenario
  GBW = 0.3  # Body weight (kg)
  tinterval = 24  # Dosing interval (hours)
  GTDOSE = 28  # Total number of doses
  kabsorb = 0.0001  # Absorption rate
  kface = 0.9999  # Adjustment factor for dose
  
  results <- list()  # List to store results for each dose
  
  # Function to run a single simulation
  run_simulation <- function(pars, DOSE) {
    GDOSEiv = DOSE * GBW * kabsorb  # Calculate intravenous dose
    
    # Define oral dose event
    Gex.oral.B <- ev(
      ID = 1,        
      time = 0,      
      amt = GDOSEiv,
      ii = tinterval,  
      addl = GTDOSE - 1, 
      cmt = "ALb", 
      replicate = FALSE
    ) 
    
    # Define dose event for Abile compartment
    GDOSEiv_Abeli = DOSE * GBW * kface  # Adjust dose
    Gex.oral.Abile <- ev(
      ID = 1,        
      time = 0,      
      amt = GDOSEiv_Abeli,
      ii = tinterval,  
      addl = GTDOSE - 1, 
      cmt = "Abile", 
      replicate = FALSE
    ) 
    
    # Combine both events
    Gex.oral <- Gex.oral.B + Gex.oral.Abile
    
    # Set output time
    Gtsamp = tgrid(0, tinterval * (GTDOSE - 1) + 24, 1)  # Set output time
    
    # Simulate the exposure scenario
    Gout.B <- 
      mod %>%  # Gestational PBPK model
      param(pars) %>%  # Update parameter list
      update(atol = 1E-6, maxsteps = 500000) %>%  # Set simulation parameters
      mrgsim_d(data = Gex.oral, tgrid = Gtsamp)  # Run the simulation
    
    # Extract variables of interest
    Goutdf.B <- cbind.data.frame(
      time = Gout.B$time, 
      CPlas = Gout.B$Plasma,     
      CL = Gout.B$Liver,     
      CK = Gout.B$Kidney,     
      CS = Gout.B$Spleen,    
      CLu = Gout.B$Lung, 
      CBr = Gout.B$Brain,
      Bal = Gout.B$BAL,
      Abile = Gout.B$Abile  # Assuming 'Abile' is a compartment name in the model
    )
    Goutdf.B <- Goutdf.B %>% filter(time > 0)  # Filter out data with time = 0
    
    return(Goutdf.B)
  }
  
  # Iterate over each dose
  for (DOSE in doses) {
    # Original simulation
    original_simulation <- run_simulation(Gpars, DOSE)
    
    # Bootstrap to calculate confidence intervals
    bootstrap_results <- replicate(n_bootstrap, {
      sampled_pars <- Gpars * exp(rnorm(length(Gpars), 0, 0.1))  # Adding small perturbations
      run_simulation(sampled_pars, DOSE)
    }, simplify = FALSE)
    
    # Combine results into a data frame
    combined_results <- do.call(rbind, bootstrap_results)
    
    # Calculate confidence intervals
    conf_intervals <- combined_results %>%
      group_by(time) %>%
      summarise(across(c(CPlas, CL, CK, CS, CLu, CBr, Bal, Abile), 
                       list(mean = ~mean(.), 
                            lower = ~quantile(., (1 - conf_level) / 2), 
                            upper = ~quantile(., 1 - (1 - conf_level) / 2))))
    
    # Store results
    results[[as.character(DOSE)]] <- list(
      "simulation" = original_simulation,
      "conf_intervals" = conf_intervals
    )
  }
  
  return(results)  # Return the results list
}

# Define a list of doses
doses <- c(50, 500, 5000)

# Calculate results
result <- pred.A(Gtheta.int, doses)

# Print results to confirm
print(result)

# Write completed data to CSV,
write.csv(result, "Rat-validation-oral-predicted-data.csv", row.names = FALSE)


