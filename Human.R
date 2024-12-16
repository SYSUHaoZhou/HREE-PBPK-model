####################### Install Required Packages #######################
library(ggplot2)
library(grid)
library(reshape2)
library(scales)
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
####################### Parameter Conversion #######################

# Body Weights
W_rats <- 0.3     # Rat body weight (kg)
W_human <- 70     # Human body weight (kg)

# Rat Parameters
RAT <- c(
  KLRESrelease   = 5.30 * 10^-19,                       
  KLRESabsorb    = 14.058133222604,                   
  KSRESrelease   = 5.30 * 10^-19,     
  KSRESabsorb    = 10.359965835688,      
  KKRESrelease   = 5.30 * 10^-19,                    
  KKRESabsorb    = 0.072500589917, 
  KBrRESrelease  = 5.30 * 10^-19,                    
  KBrRESabsorb   = 1.45,   
  KLuRESrelease  = 5.30 * 10^-19,                              
  KLuRESabsorb   = 1.45,
  KRRESrelease   = 5.30 * 10^-19,                 
  KRRESabsorb    = 0.072500093305,      
  KurineC        = 0.007050035498, 
  Kbile          = 0.00152,
  kabsorb        = 0.00002 
)

# Power Exponent
b <- -0.25

# Calculate Human Parameters Based on Rat Parameters
human <- RAT * (W_human / W_rats)^b

# Output the Calculated Human Parameters
print(human)

#######################Human PBPK Model#######################
humanPBPK.code <- '
$PARAM: Annotated Parameters

# Cardiac Output and Blood Flow
QCC            : 12.500    :L/h/kg^0.75,                    // Cardiac output (Brown, 1997, Table 22)
QLC            : 0.046     :%,                              // Fraction blood flow to liver (Brown, 1997, Table 23)
QLuC           : 1         :%,                              // Fraction blood flow to lungs (Brown, 1997, Table 23)
QKC            : 0.175     :%,                              // Fraction blood flow to kidney (Brown, 1997, Table 23)
QBrC           : 0.114     :%,                              // Fraction blood flow to brain (Brown, 1997, Table 23)
QSC            : 0.036     :%,                              // Fraction blood flow to spleen (Elliot Offman, 2015)

# Tissue Volumes
BW             : 70        :kg,                             // Body weight
VLC            : 0.026     :L/kg BW,                        // Fraction liver tissue (Brown, 1997, Table 21)
VLuC           : 0.008     :L/kg BW,                        // Fraction lung tissue (Brown, 1997, Table 21)
VKC            : 0.004     :L/kg BW,                        // Fraction kidney tissue (Brown, 1997, Table 21)
VBrC           : 0.02      :L/kg BW,                        // Fraction brain tissue (Brown, 1997, Table 21)
VSC            : 0.0027    :L/kg BW,                        // Fraction spleen (Davies and Morries, 1993)
VPlasC         : 0.0428    :L/kg BW,                        // Fraction plasma (Davies and Morries, 1993)

# Blood Volume Fraction in Organs and Tissues
BVL            : 0.11      :,                               // Blood volume in liver (Brown, 1997, Table 30)
BVBr           : 0.04      :,                               // Blood volume in brain (Brown, 1997, Table 30)
BVK            : 0.36      :,                               // Blood volume in kidney (Brown, 1997, Table 30)
BVS            : 0.3       :,                               // Blood volume in spleen (Zhe-Yi Hu, 2014)
BVLu           : 0.33      :,                               // Blood volume in lungs (Zhe-Yi Hu, 2014)
BVR            : 0.04      :,                               // Blood volume in rest of body (equals bone)

# Partition Coefficients (PC, tissue:plasma)
PL             : 0.209     :unitless,                       // Liver (Li, 2016)
PK             : 0.209     :unitless,                       // Kidney (Li, 2016)
PBr            : 0.209     :unitless,                       // Brain (Li, 2016)
PS             : 0.209     :unitless,                       // Spleen (Li, 2016)
PLu            : 0.209     :unitless,                       // Lungs (Li, 2016)
PR             : 0.209     :unitless,                       // Rest of body (Li, 2016)

# Membrane-limited Permeability Coefficient Constants
PALC           : 0.776     :mg/L,                           // Liver (Li, 2016)
PABrC          : 0.000000675 :mg/L,                         // Brain (Li, 2016)
PAKC           : 0.776     :mg/L,                           // Kidney (Li, 2016)
PASC           : 15.52     :mg/L,                           // Spleen (Li, 2016)
PALuC          : 0.776     :mg/L,                           // Lungs (Li, 2016)
PARC           : 0.0171    :mg/L,                           // Rest of body (Li, 2016)

# Endocytic Parameters; RES Represent Endocytic/Phagocytic Cells
KLRESrelease   : 1.356069e-19 :1/h,                         // Liver: Release rate constant of phagocytic cells
KLRESabsorb    : 3.596942     :1/h,                         // Liver: Uptake rate constant of phagocytic cells
KSRESrelease   : 1.356069e-19 :1/h,                         // Spleen: Release rate constant of phagocytic cells
KSRESabsorb    : 2.650722     :1/h,                         // Spleen: Uptake rate constant of phagocytic cells
KKRESrelease   : 1.356069e-19 :1/h,                         // Kidney: Release rate constant of phagocytic cells
KKRESabsorb    : 1.855015e-02 :1/h,                         // Kidney: Uptake rate constant of phagocytic cells
KLuRESrelease  : 1.356069e-19 :1/h,                         // Lung: Release rate constant of phagocytic cells
KLuRESabsorb   : 3.709999e-01 :1/h,                         // Lung: Uptake rate constant of phagocytic cells
KBrRESrelease  : 1.356069e-19 :1/h,                         // Brain: Release rate constant of phagocytic cells
KBrRESabsorb   : 3.709999e-01 :1/h,                         // Brain: Uptake rate constant of phagocytic cells
KRRESrelease   : 1.356069e-19 :1/h,                         // Rest of body: Release rate constant of phagocytic cells
KRRESabsorb    : 1.855002e-02 :1/h,                         // Rest of body: Uptake rate constant of phagocytic cells

# Excretion Parameters
KbileC         : 3.889103e-04 :L/hr/kg^0.75,                // Bile clearance
KurineC        : 1.803836e-03 :L/hr/kg^0.75,                // Urine clearance

####################### $MAIN #######################

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
double Cbile= Abile/0.2;
double Curine= Aurine/1.5;

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
double RKb      = QK*(CA-CVK) - PAK*CVK + (PAK*CKt)/PK ;
double RKt      = PAK*CVK - (PAK*CKt)/PK - KKRESUP*AKt + KKRESrelease*AKRES- Rurine; 
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
capture urine    = Aurine;
capture CONurine    = Curine;
capture CONbile    = Cbile;

'
# Load Model and Set Initial Values

humanmod <- mcode("humanPBPK", humanPBPK.code)


####################### Human PBPK Model: Oral Administration #######################

# Parameter verification (log-transformed)
Gtheta.int <- log(c(
  BW = 70,                            # Body weight (kg)
  PL = 0.209,                         # Liver partition coefficient
  PK = 0.209,                         # Kidney partition coefficient
  PBr = 0.209,                        # Brain partition coefficient
  PS = 0.209,                         # Spleen partition coefficient
  PLu = 0.209,                        # Lung partition coefficient
  PR = 0.209,                         # Rest of body partition coefficient
  PALC = 0.776,                       # Liver permeability coefficient
  PABrC = 0.000000675,                # Brain permeability coefficient
  PAKC = 0.776,                       # Kidney permeability coefficient
  PASC = 15.519935855126,             # Spleen permeability coefficient
  PALuC = 0.776,                      # Lung permeability coefficient
  PARC = 0.171,                       # Rest of body permeability coefficient
  KLRESrelease = 1.356069e-19,        # Liver: RES release rate
  KLRESabsorb = 3.596942,             # Liver: RES absorption rate
  KSRESrelease = 1.356069e-19,        # Spleen: RES release rate
  KSRESabsorb = 2.650722,             # Spleen: RES absorption rate
  KKRESrelease = 1.356069e-19,        # Kidney: RES release rate
  KKRESabsorb = 1.855015e-02,         # Kidney: RES absorption rate
  KLuRESrelease = 1.356069e-19,       # Lung: RES release rate
  KLuRESabsorb = 3.709999e-01,        # Lung: RES absorption rate
  KBrRESrelease = 1.356069e-19,       # Brain: RES release rate
  KBrRESabsorb = 3.709999e-01,        # Brain: RES absorption rate
  KRRESrelease = 1.356069e-19,        # Rest of body: RES release rate
  KRRESabsorb = 1.855002e-02,         # Rest of body: RES absorption rate
  KbileC = 3.889103e-04,              # Bile clearance
  KurineC = 1.803836e-03              # Urine clearance
))

# Function to predict results for individuals
pred.A <- function(Gpars, individuals, n_bootstrap = 100, conf_level = 0.95) {
  # Gpars: Input parameters, individuals: List of individual info including weight, age, dose
  # n_bootstrap: Number of bootstrap simulations, conf_level: Confidence level for intervals
  
  Gpars <- exp(Gpars)  # Convert parameters from log-domain back to normal domain
  
  # Fixed parameters for exposure scenario
  tinterval = 24        # Dosing interval in hours
  kabsorb = 0.00002     # Absorption rate
  kface = 0.99998       # Dose adjustment factor
  
  results <- list()         # To store results for each individual
  summary_data <- list()    # To store summary data
  
  # Function to run a single simulation
  run_simulation <- function(pars, individual, id) {
    GTDOSE <- individual$age * 365                          # Calculate total dose based on age
    GDOSEiv <- individual$dose * individual$weight * kabsorb # IV dose calculation
    
    # Oral dose event entering liver
    Gex.oral.B <- ev(
      ID = 1, 
      time = 0, 
      amt = GDOSEiv, 
      ii = tinterval, 
      addl = GTDOSE - 1, 
      cmt = "ALb", 
      replicate = FALSE
    )
    
    # Unabsorbed rare earth directly excreted in bile
    GDOSEiv_Abeli <- individual$dose * individual$weight * kface # Adjusted dose
    Gex.oral.Abile <- ev(
      ID = 1, 
      time = 0, 
      amt = GDOSEiv_Abeli, 
      ii = tinterval, 
      addl = GTDOSE - 1, 
      cmt = "Abile", 
      replicate = FALSE
    )
    
    # Combine the two events
    Gex.oral <- Gex.oral.B + Gex.oral.Abile
    
    # Set output time
    Gtsamp <- tgrid(0, tinterval * (GTDOSE - 1) + 24, 500)
    
    # Simulate exposure scenario
    Gout.B <- humanmod %>%
      param(pars) %>%
      update(atol = 1E-6, maxsteps = 500000) %>%
      mrgsim_d(data = Gex.oral, tgrid = Gtsamp)
    
    # Extract the variables of interest and add individual ID column
    Goutdf.B <- cbind.data.frame(
      ID = id,
      time = Gout.B$time,
      CPlas = Gout.B$Plasma,
      Bal = Gout.B$BAL,
      Aurine = Gout.B$Aurine,
      Curine = Gout.B$CONurine
    )
    Goutdf.B <- Goutdf.B %>% filter(time > 0)  # Filter out time = 0
    
    return(Goutdf.B)
  }
  
  # Loop through each individual
  for (i in seq_along(individuals)) {
    individual <- individuals[[i]]
    GTDOSE <- individual$age * 365  # Calculate total dose based on age
    
    # Original simulation
    original_simulation <- run_simulation(Gpars, individual, id = i)
    
    # Bootstrap to calculate confidence intervals
    bootstrap_results <- replicate(n_bootstrap, {
      sampled_pars <- Gpars * exp(rnorm(length(Gpars), 0, 0.1))  # Perturb parameters
      run_simulation(sampled_pars, individual, id = i)
    }, simplify = FALSE)
    
    # Combine bootstrap results
    combined_results <- do.call(rbind, bootstrap_results)
    
    # Calculate confidence intervals
    conf_intervals <- combined_results %>%
      group_by(time) %>%
      summarise(across(c(Curine), 
                       list(mean = ~mean(.), 
                            lower = ~quantile(., (1 - conf_level) / 2), 
                            upper = ~quantile(., 1 - (1 - conf_level) / 2))))
    
    # Store results
    results[[paste0("Individual_", i)]] <- list(
      "simulation" = original_simulation,
      "conf_intervals" = conf_intervals
    )
    
    # Store summary data for the individual
    summary_data[[paste0("Individual_", i)]] <- list(
      "GTDOSE" = GTDOSE,
      "dose" = individual$dose,
      "weight" = individual$weight
    )
  }
  
  return(list("results" = results, "summary" = summary_data))  # Return results and summary
}

# Define individual information list
individuals <- list(
  list(weight = 70, age = 50, dose = 0.000004104553),
  list(weight = 70, age = 50, dose = 0.000003830916),
  list(weight = 70, age = 50, dose = 0.000000621139)
  # Add more individuals if necessary
)

# Calculate results
resulthuman <- pred.A(Gtheta.int, individuals)

# View summary data
summary_data <- resulthuman$summary
print(resulthuman)

# Write completed data to CSV, manually extract the concentration value of the target organ
write.csv(resulthuman, "human-50years-oral-predicted-data.csv", row.names = FALSE)


#######################Calculate the Uncertain Factor(UF)#######################
##########calculate concentration of human liver
# Load Model and Set Initial Values

humanmod <- mcode("humanPBPK", humanPBPK.code)


# Parameter verification (log-transformed)
Gtheta.int <- log(c(
  BW = 70,                            # Body weight (kg)
  PL = 0.209,                         # Liver partition coefficient
  PK = 0.209,                         # Kidney partition coefficient
  PBr = 0.209,                        # Brain partition coefficient
  PS = 0.209,                         # Spleen partition coefficient
  PLu = 0.209,                        # Lung partition coefficient
  PR = 0.209,                         # Rest of body partition coefficient
  PALC = 0.776,                       # Liver permeability coefficient
  PABrC = 0.000000675,                # Brain permeability coefficient
  PAKC = 0.776,                       # Kidney permeability coefficient
  PASC = 15.519935855126,             # Spleen permeability coefficient
  PALuC = 0.776,                      # Lung permeability coefficient
  PARC = 0.171,                       # Rest of body permeability coefficient
  KLRESrelease = 1.356069e-19,        # Liver: RES release rate
  KLRESabsorb = 3.596942,             # Liver: RES absorption rate
  KSRESrelease = 1.356069e-19,        # Spleen: RES release rate
  KSRESabsorb = 2.650722,             # Spleen: RES absorption rate
  KKRESrelease = 1.356069e-19,        # Kidney: RES release rate
  KKRESabsorb = 1.855015e-02,         # Kidney: RES absorption rate
  KLuRESrelease = 1.356069e-19,       # Lung: RES release rate
  KLuRESabsorb = 3.709999e-01,        # Lung: RES absorption rate
  KBrRESrelease = 1.356069e-19,       # Brain: RES release rate
  KBrRESabsorb = 3.709999e-01,        # Brain: RES absorption rate
  KRRESrelease = 1.356069e-19,        # Rest of body: RES release rate
  KRRESabsorb = 1.855002e-02,         # Rest of body: RES absorption rate
  KbileC = 3.889103e-04,              # Bile clearance
  KurineC = 1.803836e-03              # Urine clearance
))

# Function to predict results for individuals
pred.A <- function(Gpars, individuals, n_bootstrap = 100, conf_level = 0.95) {
  # Gpars: Input parameters, individuals: List of individual info including weight, age, dose
  # n_bootstrap: Number of bootstrap simulations, conf_level: Confidence level for intervals
  
  Gpars <- exp(Gpars)  # Convert parameters from log-domain back to normal domain
  
  # Fixed parameters for exposure scenario
  tinterval = 24        # Dosing interval in hours
  kabsorb = 0.00002     # Absorption rate
  kface = 0.99998       # Dose adjustment factor
  
  results <- list()         # To store results for each individual
  summary_data <- list()    # To store summary data
  
  # Function to run a single simulation
  run_simulation <- function(pars, individual, id) {
    GTDOSE <- individual$age * 365                          # Calculate total dose based on age
    GDOSEiv <- individual$dose * individual$weight * kabsorb # IV dose calculation
    
    # Oral dose event entering liver
    Gex.oral.B <- ev(
      ID = 1, 
      time = 0, 
      amt = GDOSEiv, 
      ii = tinterval, 
      addl = GTDOSE - 1, 
      cmt = "ALb", 
      replicate = FALSE
    )
    
    # Unabsorbed rare earth directly excreted in bile
    GDOSEiv_Abeli <- individual$dose * individual$weight * kface # Adjusted dose
    Gex.oral.Abile <- ev(
      ID = 1, 
      time = 0, 
      amt = GDOSEiv_Abeli, 
      ii = tinterval, 
      addl = GTDOSE - 1, 
      cmt = "Abile", 
      replicate = FALSE
    )
    
    # Combine the two events
    Gex.oral <- Gex.oral.B + Gex.oral.Abile
    
    # Set output time
    Gtsamp <- tgrid(0, tinterval * (GTDOSE - 1) + 24, 500)
    
    # Simulate exposure scenario
    Gout.B <- humanmod %>%
      param(pars) %>%
      update(atol = 1E-6, maxsteps = 500000) %>%
      mrgsim_d(data = Gex.oral, tgrid = Gtsamp)
    
    # Extract the variables of interest and add individual ID column
    Goutdf.B <- cbind.data.frame(
      ID = id,
      time = Gout.B$time,
      CPlas = Gout.B$Plasma,
      Bal = Gout.B$BAL,
      Aurine = Gout.B$Aurine,
      Curine = Gout.B$CONurine
    )
    Goutdf.B <- Goutdf.B %>% filter(time > 0)  # Filter out time = 0
    
    return(Goutdf.B)
  }
  
  # Loop through each individual
  for (i in seq_along(individuals)) {
    individual <- individuals[[i]]
    GTDOSE <- individual$age * 365  # Calculate total dose based on age
    
    # Original simulation
    original_simulation <- run_simulation(Gpars, individual, id = i)
    
    # Bootstrap to calculate confidence intervals
    bootstrap_results <- replicate(n_bootstrap, {
      sampled_pars <- Gpars * exp(rnorm(length(Gpars), 0, 0.1))  # Perturb parameters
      run_simulation(sampled_pars, individual, id = i)
    }, simplify = FALSE)
    
    # Combine bootstrap results
    combined_results <- do.call(rbind, bootstrap_results)
    
    # Calculate confidence intervals
    conf_intervals <- combined_results %>%
      group_by(time) %>%
      summarise(across(c(Curine), 
                       list(mean = ~mean(.), 
                            lower = ~quantile(., (1 - conf_level) / 2), 
                            upper = ~quantile(., 1 - (1 - conf_level) / 2))))
    
    # Store results
    results[[paste0("Individual_", i)]] <- list(
      "simulation" = original_simulation,
      "conf_intervals" = conf_intervals
    )
    
    # Store summary data for the individual
    summary_data[[paste0("Individual_", i)]] <- list(
      "GTDOSE" = GTDOSE,
      "dose" = individual$dose,
      "weight" = individual$weight
    )
  }
  
  return(list("results" = results, "summary" = summary_data))  # Return results and summary
}

# Define individual information list
individuals <- list(
  list(weight = 70, age = 50, dose = 1)
  # Add more individuals if necessary
)

# Calculate results
resulthuman <- pred.A(Gtheta.int, individuals)

# View summary data
summary_data <- resulthuman$summary
print(resulthuman)

# Write completed data to CSV, manually extract the concentration value of the target organ
write.csv(resulthuman, "human-oral-50years-1mg-kg-bw.csv", row.names = FALSE)

##########calculate concentration of mice liver
micegdPBPK.code <- '
$PARAM @annotated
// Cardiac Output and Blood flow
QCC            : 16.5    :L/h/kg^0.75,                    // Cardiac output (Brown, 1997, Table 22)
QLC            : 0.02    :%                               // Fraction blood flow to liver (Brown, 1997, Table 23)
QLuC           : 1       :%                               // Fraction blood flow to lung (Brown, 1997, Table 23)
QKC            : 0.091   :%                               // Fraction blood flow to kidney (Brown, 1997, Table 23)
QBrC           : 0.033   :%                               // Fraction blood flow to brain (Brown, 1997, Table 23)
QSC            : 0.011   :%                               // Fraction blood flow to spleen (Davies and Morries, 1993)

// Tissue Volume
BW             : 0.025   :kg                              // Body weight
VLC            : 0.055   :L/kg BW                         // Fraction liver tissue (Brown, 1997, Table 21)
VLuC           : 0.007   :L/kg BW                         // Fraction lung tissue (Brown, 1997, Table 21)
VKC            : 0.017   :L/kg BW                         // Fraction kidney tissue (Brown, 1997, Table 21)
VBrC           : 0.017   :L/kg BW                         // Fraction brain tissue (Brown, 1997, Table 21)
VSC            : 0.005   :L/kg BW                         // Fraction spleen (Davies and Morries, 1993)
VPlasC         : 0.0355  :L/kg BW                         // Fraction plasma (Davies and Morries, 1993)

// Blood volume fraction in organs and tissues
BVL            : 0.31    :%,                              // Liver (Brown, 1997, Table 30)
BVBr           : 0.03    :%,                              // Brain (Brown, 1997, Table 30)
BVK            : 0.24    :%,                              // Kidney (Brown, 1997, Table 30)
BVS            : 0.17    :%,                              // Spleen (Brown, 1997, Table 30)
BVLu           : 0.5     :%,                              // Lungs (Brown, 1997, Table 30)
BVR            : 0.04    :%,                              // Rest of body

// Partition coefficients (PC, tissue:plasma)
PL             : 0.209   :unitless                        // Liver (Li, 2016)
PK             : 0.209   :unitless                        // Kidney (Li, 2016)
PBr            : 0.209   :unitless                        // Brain (Li, 2016)
PS             : 0.209   :unitless                        // Spleen (Li, 2016)
PLu            : 0.209   :unitless                        // Lungs (Li, 2016)
PR             : 0.209   :unitless                        // Rest of body (Li, 2016)

// Membrane-limited permeability coefficient constants
PALC           : 0.776   :mg/L                            // Liver (Li, 2016)
PABrC          : 0.000000675 :mg/L                        // Brain (Li, 2016)
PAKC           : 0.776   :mg/L                            // Kidney (Li, 2016)
PASC           : 15.52   :mg/L                            // Spleen (Li, 2016)
PALuC          : 0.776   :mg/L                            // Lung (Li, 2016)
PARC           : 0.0171  :mg/L                            // Rest of body (Li, 2016)

// Endocytic parameters; RES represent endocytic/phagocytic cells
KLRESrelease   : 9.864412e-19  :1/h                       // Liver release rate constant of phagocytic cells
KLRESabsorb    : 26.16513      :1/h                       // Liver maximum uptake rate
KSRESrelease   : 9.864412e-19  :1/h                       // Spleen release rate constant of phagocytic cells
KSRESabsorb    : 19.28207      :1/h                       // Spleen maximum uptake rate
KKRESrelease   : 9.864412e-19  :1/h                       // Kidney release rate constant of phagocytic cells
KKRESabsorb    : 0.1349388     :1/h                       // Kidney maximum uptake rate
KLuRESrelease  : 9.864412e-19  :1/h                       // Lung release rate constant of phagocytic cells
KLuRESabsorb   : 2.698754      :1/h                       // Lung maximum uptake rate
KBrRESrelease  : 9.864412e-19  :1/h                       // Brain release rate constant of phagocytic cells
KBrRESabsorb   : 2.698754      :1/h                       // Brain maximum uptake rate
KRRESrelease   : 9.864412e-19  :1/h                       // Rest of body release rate constant of phagocytic cells
KRRESabsorb    : 0.1349379     :1/h                       // Rest of body maximum uptake rate

// Excretion parameters
KbileC         : 2.829039e-03  :L/hr/kg^0.75              // Bile clearance
KurineC        : 1.312159e-02  :L/hr/kg^0.75              // Urine clearance

$MAIN
// Blood flow fractions
double QRC = 1 - (QLC + QKC + QSC + QBrC);                // Fraction of blood flow to rest of body
double VRC = 1 - (VLC + VLuC + VKC + VSC + VBrC + VPlasC); // Tissue volume of rest of body

// Cardiac output and tissue-specific blood flows
double QC  = QCC * pow(BW, 0.75);                         // Cardiac output
double QL  = QC * QLC;                                    // Blood flow to liver
double QBr = QC * QBrC;                                   // Blood flow to brain
double QK  = QC * QKC;                                    // Blood flow to kidney
double QS  = QC * QSC;                                    // Blood flow to spleen
double QR  = QC * QRC;                                    // Blood flow to the rest of body

// Tissue volumes
double VL  = BW * VLC;                                    // Liver volume
double VBr = BW * VBrC;                                   // Brain volume
double VK  = BW * VKC;                                    // Kidney volume
double VS  = BW * VSC;                                    // Spleen volume
double VLu = BW * VLuC;                                   // Lung volume
double VR  = BW * VRC;                                    // Rest of body volume
double VPlasma = BW * VPlasC;                             // Plasma volume

// Capillary blood volumes and tissue volumes in organs
double VLb = VL * BVL;
double VLt = VL - VLb;
double VBrb = VBr * BVBr;
double VBrt = VBr - VBrb;
double VKb = VK * BVK;
double VKt = VK - VKb;
double VSb = VS * BVS;
double VSt = VS - VSb;
double VLub = VLu * BVLu;
double VLut = VLu - VLub;
double VRb = VR * BVR;
double VRt = VR - VRb;

// Permeability coefficient-surface area cross-product (mg/h)
double PAL   = PALC * QL;                                 // Liver
double PABr  = PABrC * QBr;                               // Brain
double PAK   = PAKC * QK;                                 // Kidneys
double PAS   = PASC * QS;                                 // Spleen
double PALu  = PALuC * QC;                                // Lungs
double PAR   = PARC * QR;                                 // Rest of body

// Endocytosis rates
double KLRESUP = KLRESabsorb;
double KSRESUP = KSRESabsorb;
double KKRESUP = KKRESabsorb;
double KLuRESUP = KLuRESabsorb;
double KRRESUP = KRRESabsorb;
double KBrRESUP = KBrRESabsorb;

double Kurine  = KurineC * pow(BW, 0.75);
double Kbile   = KbileC * pow(BW, 0.75);

$CMT AA AV ALub ALut ALuRES ABrb ABrt ARb ARt ARRES Abile ABrRES AKb AKt AKRES Aurine ASb ASt ASRES ALb ALt ALRES ADOSE

$ODE
// Mass balance equations
double APlasma = AA + AV;
double ALung = ALub + ALut + ALuRES;
double ALungt = ALut + ALuRES;
double Arestall = ARb + ARt + ARRES;
double Aresttissue = ARt + ARRES;
double AKidney = AKb + AKt + AKRES;
double AKidneyt = AKt + AKRES;
double ABrain = ABrb + ABrt + ABrRES;
double ABraint = ABrt + ABrRES;
double ASpleen = ASb + ASt + ASRES;
double ASpleent = ASt + ASRES;
double ALiver = ALb + ALt + ALRES;
double ALivert = ALt + ALRES;

double CA = AA / (VPlasma * 0.2);
double CV = AV / (VPlasma * 0.8);
double CPlasma = APlasma / VPlasma;
double CVLu = ALub / VLub;
double CLut = ALut / VLut;
double CLung = (ALub + ALut + ALuRES) / VLu;
double CLungt = (ALut + ALuRES) / VLut;
double CVBr = ABrb / VBrb;
double CBrt = ABrt / VBrt;
double CBrain = (ABrain + ABrRES + ABrb) / VBr;
double CBraint = (ABrt + ABrRES) / VBrt;
double CVR = ARb / VRb;
double CRt = ARt / VRt;
double Crestall = (ARb + ARt + ARRES) / VR;
double Cresttissue = (ARt + ARRES) / VRt;
double CVK = AKb / VKb;
double CKt = AKt / VKt;
double CKidney = (AKb + AKt + AKRES) / VK;
double CKidneyt = (AKt + AKRES) / VKt;
double CVS = ASb / VSb;
double CSt = ASt / VSt;
double CSpleen = (ASb + ASt + ASRES) / VS;
double CSpleent = (ASt + ASRES) / VSt;
double CVL = ALb / VLb;
double CLt = ALt / VLt;
double CLiver = (ALb + ALt + ALRES) / VL;
double CLivert = (ALt + ALRES) / VLt;

double Rbile = Kbile * ALt;
double RA = QC * CVLu - QC * CA;
double RV = QL * CVL + QBr * CVBr + QK * CVK + QR * CVR - QC * CV;
double RLub = QC * (CV - CVLu) - PALu * CVLu + (PALu * CLut) / PLu;
double RLut = PALu * CVLu - (PALu * CLut) / PLu - KLuRESUP * ALut + KLuRESrelease * ALuRES;
double RLuRES = KLuRESUP * ALut - KLuRESrelease * ALuRES;
double RBrb = QBr * (CA - CVBr) - PABr * CVBr + (PABr * CBrt) / PBr;
double RBrt = PABr * CVBr - (PABr * CBrt) / PBr + KBrRESrelease * ABrRES;
double RBrRES = KBrRESUP * ABrt - KBrRESrelease * ABrRES;
double RRb = QR * (CA - CVR) - PAR * CVR + (PAR * CRt) / PR;
double RRt = PAR * CVR - (PAR * CRt) / PR - KRRESUP * ARt + KRRESrelease * ARRES;
double RRRES = KRRESUP * ARt - KRRESrelease * ARRES;
double Rurine = Kurine * AKt;
double RKb = QK * (CA - CVK) - PAK * CVK + (PAK * CKt) / PK ;
double RKt = PAK * CVK - (PAK * CKt) / PK - KKRESUP * AKt + KKRESrelease * AKRES- Rurine;
double RKRES = KKRESUP * AKt - KKRESrelease * AKRES;
double RSb = QS * (CA - CVS) - PAS * CVS + (PAS * CSt) / PS;
double RSt = PAS * CVS - (PAS * CSt) / PS - KSRESUP * ASt + KSRESrelease * ASRES;
double RSRES = KSRESUP * ASt - KSRESrelease * ASRES;
double RLb = QL * (CA - CVL) + QS * CVS - PAL * CVL + (PAL * CLt) / PL;
double RLt = PAL * CVL - (PAL * CLt) / PL - Rbile - KLRESUP * ALb + KLRESrelease * ALRES;
double RLRES = KLRESUP * ALb - KLRESrelease * ALRES;

$TABLE
capture Lung = CLung;
capture Liver = CLiver;
capture Kidney = CKidney;
capture Spleen = CSpleen;
capture Plasma = CPlasma;
capture Brain = CBrain;
capture BAL = Bal;
capture CONbile = Abile / 0.002;
capture CONurine = Aurine / 0.002;
'

# Load model and set initial values
mod <- mcode("micegdPBPK.code", micegdPBPK.code)

Gtheta.int <- log(c(
  PL = 0.209, PK = 0.209, PBr = 0.209, PS = 0.209, PLu = 0.209, PR = 0.209,
  PALC = 0.776, PABrC = 0.000000675, PAKC = 0.776, PASC = 15.519935855126,
  PALuC = 0.776, PARC = 0.171,
  KLRESrelease = 9.864412e-19, KLRESabsorb = 2.616513e+01,
  KSRESrelease = 9.864412e-19, KSRESabsorb = 1.928207e+01,
  KKRESrelease = 9.864412e-19, KKRESabsorb = 1.349388e-01,
  KLuRESrelease = 9.864412e-19, KLuRESabsorb = 2.698754e+00,
  KBrRESrelease = 9.864412e-19, KBrRESabsorb = 2.698754e+00,
  KRRESrelease = 9.864412e-19, KRRESabsorb = 1.349379e-01,
  KbileC = 2.829039e-03, KurineC = 1.312159e-02
))

# Simulation Function 
pred.A <- function(Gpars, doses, n_bootstrap = 1000, conf_level = 0.95) {
  # Gpars: Model parameters
  # doses: List of doses
  # n_bootstrap: Number of bootstrap simulations
  # conf_level: Confidence level for intervals
  
  # Convert parameters from log-domain to normal domain
  Gpars <- exp(Gpars)
  
  # Fixed parameters for the exposure scenario
  GBW = 0.025  # Body weight (kg)
  tinterval = 24  # Dosing interval (hours)
  GTDOSE = 28  # Total number of doses
  kabsorb = 0.00002  # Absorption rate
  kface = 0.99998  # Dose adjustment factor
  
  results <- list()  # To store the results for each dose
  
  # Function to run a single simulation
  run_simulation <- function(pars, DOSE) {
    GDOSEiv = DOSE * GBW * kabsorb  # Calculate intravenous dose
    
    # Oral dose event for liver compartment
    Gex.oral.B <- ev(
      ID = 1, time = 0, amt = GDOSEiv, ii = tinterval,
      addl = GTDOSE - 1, cmt = "ALb", replicate = FALSE
    )
    
    # Adjusted dose for bile compartment
    GDOSEiv_Abeli = DOSE * GBW * kface
    Gex.oral.Abile <- ev(
      ID = 1, time = 0, amt = GDOSEiv_Abeli, ii = tinterval,
      addl = GTDOSE - 1, cmt = "Abile", replicate = FALSE
    )
    
    # Combine both events
    Gex.oral <- Gex.oral.B + Gex.oral.Abile
    
    # Set output time
    Gtsamp <- tgrid(0, tinterval * (GTDOSE - 1) + 24, 1)
    
    # Simulate the exposure scenario
    Gout.B <- mod %>%
      param(pars) %>%
      update(atol = 1E-6, maxsteps = 500000) %>%
      mrgsim_d(data = Gex.oral, tgrid = Gtsamp)
    
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
      Abile = Gout.B$Abile,
      Aurine = Gout.B$Aurine,
      Cbile = Gout.B$CONbile,
      Curine = Gout.B$CONurine
    )
    Goutdf.B <- Goutdf.B %>% filter(time > 0)  # Remove data where time = 0
    return(Goutdf.B)
  }
  
  # Loop through each dose and run the simulation
  for (DOSE in doses) {
    # Original simulation
    original_simulation <- run_simulation(Gpars, DOSE)
    
    # Bootstrap to calculate confidence intervals
    bootstrap_results <- replicate(n_bootstrap, {
      sampled_pars <- Gpars * exp(rnorm(length(Gpars), 0, 0.1))  # Add small perturbations
      run_simulation(sampled_pars, DOSE)
    }, simplify = FALSE)
    
    # Combine results into a data frame
    combined_results <- do.call(rbind, bootstrap_results)
    
    # Calculate confidence intervals
    conf_intervals <- combined_results %>%
      group_by(time) %>%
      summarise(across(c(CPlas, CL, CK, CS, CLu, CBr, Bal, Abile, Aurine, Cbile, Curine),
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

# Dose List and Execution
doses <- c(1)  # Define the list of doses

# Run the simulation for the defined doses
resultmice <- pred.A(Gtheta.int, doses)

# Print the simulation results
print(resultmice)

# Write completed data to CSV,
write.csv(resultmice, "mice-oral-28days-1mg-kg-bw.csv", row.names = FALSE)

##########read data
data <- read_csv("human-oral-50years-1mg-kg-bw.csv")
time <- data$results.Individual_1.simulation.time
concentration <- data$results.Individual_1.simulation.CL

# calculate AUC
calculate_auc <- function(time, concentration) {
  auc <- 0
  for (i in 1:(length(time) - 1)) {
    auc <- auc + 0.5 * (concentration[i] + concentration[i+1]) * (time[i+1] - time[i])
  }
  return(auc)
}

auc_value_human <- calculate_auc(time, concentration)

# print AUC
print(paste("The AUC value is:", auc_value_human))


data <- read_csv("mice-oral-28days-1mg-kg-bw.csv")

# read data
time <- data$X1.simulation.time
concentration <- data$X1.simulation.CL

# calculate AUC
calculate_auc <- function(time, concentration) {
  auc <- 0
  for (i in 1:(length(time) - 1)) {
    auc <- auc + 0.5 * (concentration[i] + concentration[i+1]) * (time[i+1] - time[i])
  }
  return(auc)
}

auc_value_mice <- calculate_auc(time, concentration)

# print AUC
print(paste("The AUC value is:", auc_value_mice))

# calculate UF
UF <- auc_value_mice / auc_value_human
print(UF)