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

############################# Parameter Conversion #############################

# Body Weights
W_rats <- 0.3    # Rat body weight (kg)
W_mice <- 0.025  # Mouse body weight (kg)

# Rat Parameters
RAT <- c(
  KLRESrelease = 5.30 * 10^-19,    # Liver release rate constant
  KLRESabsorb = 14.058133222604,   # Liver absorption rate constant
  KSRESrelease = 5.30 * 10^-19,    # Spleen release rate constant
  KSRESabsorb = 10.359965835688,   # Spleen absorption rate constant
  KKRESrelease = 5.30 * 10^-19,    # Kidney release rate constant
  KKRESabsorb = 0.072500589917,    # Kidney absorption rate constant
  KBrRESrelease = 5.30 * 10^-19,   # Brain release rate constant
  KBrRESabsorb = 1.45,             # Brain absorption rate constant
  KLuRESrelease = 5.30 * 10^-19,   # Lung release rate constant
  KLuRESabsorb = 1.45,             # Lung absorption rate constant
  KRRESrelease = 5.30 * 10^-19,    # Rest of body release rate constant
  KRRESabsorb = 0.072500093305,    # Rest of body absorption rate constant
  KurineC = 0.007050035498,        # Urine clearance rate
  Kbile = 0.00152                  # Bile clearance rate
)

# Power exponent
b <- -0.25

# Calculate mouse liver parameters based on rat parameters
mice <- RAT * (W_mice / W_rats)^b

# Output the calculated parameters for mice
print(mice)


############################# Mouse PBPK Model #############################
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
BW             : 0.025   :kg                              // Body weight (Brown, 1997)
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
double RKb = QK * (CA - CVK) - PAK * CVK + (PAK * CKt) / PK - Rurine;
double RKt = PAK * CVK - (PAK * CKt) / PK - KKRESUP * AKt + KKRESrelease * AKRES - Rurine;
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

############################# Mouse PBPK Model: Oral Administration  #############################
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
doses <- c(50, 100, 200)  # Define the list of doses

# Run the simulation for the defined doses
resultmice <- pred.A(Gtheta.int, doses)

# Print the simulation results
print(resultmice)

# Write completed data to CSV,
write.csv(resultmice, "mice-28-oral-predicted-data.csv", row.names = FALSE)



############################# Calculate external exposure based on internal exposure #############################
# Load model and set initial values
mod <- mcode("micegdPBPK.code", micegdPBPK.code)

# Initialize parameters
Gtheta.int <- log(c(
  PL = 0.209, PK = 0.209, PBr = 0.209, PS = 0.209, PLu = 0.209, PR = 0.209,
  PALC = 0.776, PABrC = 0.000000675, PAKC = 0.776, PASC = 15.519935855126,
  PALuC = 0.776, PARC = 0.171, KLRESrelease = 9.864412e-19, KLRESabsorb = 2.616513e+01,
  KSRESrelease = 9.864412e-19, KSRESabsorb = 1.928207e+01, KKRESrelease = 9.864412e-19,
  KKRESabsorb = 1.349388e-01, KLuRESrelease = 9.864412e-19, KLuRESabsorb = 2.698754e+00,
  KBrRESrelease = 9.864412e-19, KBrRESabsorb = 2.698754e+00, KRRESrelease = 9.864412e-19,
  KRRESabsorb = 1.349379e-01, KbileC = 2.829039e-03, KurineC = 1.312159e-02
))

# Simulation function
pred.A <- function(Gpars, doses, n_bootstrap = 1000, conf_level = 0.95) {
  # Convert parameters from log scale to normal scale
  Gpars <- exp(Gpars)
  
  # Define exposure scenario parameters
  GBW = 0.025  # Body weight (kg)
  tinterval = 24  # Dose interval (hours)
  GTDOSE = 28  # Total number of doses
  kabsorb = 0.00002  # Absorption rate
  kface = 0.99998  # Adjustment coefficient
  
  results <- list()
  
  # Function to run a single simulation
  run_simulation <- function(pars, DOSE) {
    GDOSEiv = DOSE * GBW * kabsorb  # Calculate intravenous dose
    
    # Define oral dose event
    Gex.oral.B <- ev(
      ID = 1, time = 0, amt = GDOSEiv, ii = tinterval, addl = GTDOSE - 1,
      cmt = "ALb", replicate = FALSE
    )
    
    # Define bile compartment event
    GDOSEiv_Abeli = DOSE * GBW * kface  # Adjusted dose
    Gex.oral.Abile <- ev(
      ID = 1, time = 0, amt = GDOSEiv_Abeli, ii = tinterval, addl = GTDOSE - 1,
      cmt = "Abile", replicate = FALSE
    )
    
    # Combine the two events
    Gex.oral <- Gex.oral.B + Gex.oral.Abile
    
    # Set output time grid
    Gtsamp = tgrid(0, tinterval * (GTDOSE - 1) + 24, 1)
    
    # Run the simulation
    Gout.B <- mod %>%
      param(pars) %>%
      update(atol = 1E-6, maxsteps = 500000) %>%
      mrgsim_d(data = Gex.oral, tgrid = Gtsamp)
    
    # Extract variables of interest
    Goutdf.B <- cbind.data.frame(
      time = Gout.B$time,
      CPlas = Gout.B$Plasma, CL = Gout.B$Liver, CK = Gout.B$Kidney,
      CS = Gout.B$Spleen, CLu = Gout.B$Lung, CBr = Gout.B$Brain,
      Bal = Gout.B$BAL, Abile = Gout.B$Abile, Aurine = Gout.B$Aurine,
      Cbile = Gout.B$CONbile, Curine = Gout.B$CONurine
    )
    Goutdf.B <- Goutdf.B %>% filter(time > 0)
    
    return(Goutdf.B)
  }
  
  # Iterate through each dose
  for (DOSE in doses) {
    # Original simulation
    original_simulation <- run_simulation(Gpars, DOSE)
    
    # Bootstrap to calculate confidence intervals
    bootstrap_results <- replicate(n_bootstrap, {
      sampled_pars <- Gpars * exp(rnorm(length(Gpars), 0, 0.1))
      run_simulation(sampled_pars, DOSE)
    }, simplify = FALSE)
    
    # Combine results into a data frame
    combined_results <- do.call(rbind, bootstrap_results)
    
    # Calculate confidence intervals
    conf_intervals <- combined_results %>%
      group_by(time) %>%
      summarise(across(
        c(CPlas, CL, CK, CS, CLu, CBr, Bal, Abile, Aurine, Cbile, Curine),
        list(mean = ~mean(.), lower = ~quantile(., (1 - conf_level) / 2), upper = ~quantile(., 1 - (1 - conf_level) / 2))
      ))
    
    # Store results
    results[[as.character(DOSE)]] <- list(
      "simulation" = original_simulation,
      "conf_intervals" = conf_intervals
    )
  }
  
  return(results)
}


#############Brain
# Function to find the dose with brain concentration closest to the target concentration using a refining search approach
find_closest_brain_dose <- function(dose_range, target_concentration, refinement_steps = c(1, 0.1, 0.01, 0.001,0.0001,0.00001,0.000001,0.0000001)) {
  dose_brain_concentrations <- data.frame(Dose = numeric(), brainConcentration = numeric())
  
  # Start with initial coarse step and then refine
  for (step in refinement_steps) {
    refined_dose_brain_concentrations <- data.frame(Dose = numeric(), brainConcentration = numeric())
    
    for (current_dose in seq(dose_range[1], dose_range[2], by = step)) {
      # Run simulation with the current dose
      simulation_results <- pred.A(Gtheta.int, doses = c(current_dose), n_bootstrap = 100, conf_level = 0.95)
      
      # Extract brain concentration at the last time point
      brain_concentration <- tail(simulation_results[[as.character(current_dose)]]$simulation$CBr, n = 1)
      
      # Store the dose and brain concentration
      refined_dose_brain_concentrations <- rbind(refined_dose_brain_concentrations, data.frame(Dose = current_dose, brainConcentration = brain_concentration))
    }
    
    # Find the dose with brain concentration closest to the target concentration
    closest_row <- refined_dose_brain_concentrations[which.min(abs(refined_dose_brain_concentrations$brainConcentration - target_concentration)), ]
    
    # Set new search range around the closest dose for next refinement
    dose_range <- c(max(closest_row$Dose - step, 0), closest_row$Dose + step)
    
    # Store the results from this refinement step
    dose_brain_concentrations <- rbind(dose_brain_concentrations, refined_dose_brain_concentrations)
  }
  
  # Find the final closest dose after all refinements
  closest_dose <- dose_brain_concentrations[which.min(abs(dose_brain_concentrations$brainConcentration - target_concentration)), ]
  
  return(closest_dose)
}

# Parameters for the brain dose search
dose_range <- c(0, 5)  # Single dose range to search
multiple_target_concentrations_brain <- c("targer concentrarion in brain")  # List of target concentrations for brain

# Find the closest dose for each target concentration for brain
closest_brain_dose_results <- lapply(multiple_target_concentrations_brain, function(target_concentration) {
  find_closest_brain_dose(dose_range, target_concentration)
})

# Display the results
print(closest_brain_dose_results)

#############Liver
# Function to find the dose with liver concentration closest to the target concentration using a refining search approach
find_closest_liver_dose <- function(dose_range, target_concentration, refinement_steps = c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001,0.00000001)) {
  dose_liver_concentrations <- data.frame(Dose = numeric(), LiverConcentration = numeric())
  
  # Start with initial coarse step and then refine
  for (step in refinement_steps) {
    refined_dose_liver_concentrations <- data.frame(Dose = numeric(), LiverConcentration = numeric())
    
    for (current_dose in seq(dose_range[1], dose_range[2], by = step)) {
      # Run simulation with the current dose
      simulation_results <- pred.A(Gtheta.int, doses = c(current_dose), n_bootstrap = 100, conf_level = 0.95)
      
      # Extract liver concentration at the last time point
      liver_concentration <- tail(simulation_results[[as.character(current_dose)]]$simulation$CL, n = 1)
      
      # Store the dose and liver concentration
      refined_dose_liver_concentrations <- rbind(refined_dose_liver_concentrations, data.frame(Dose = current_dose, LiverConcentration = liver_concentration))
    }
    
    # Find the dose with liver concentration closest to the target concentration
    closest_row <- refined_dose_liver_concentrations[which.min(abs(refined_dose_liver_concentrations$LiverConcentration - target_concentration)), ]
    
    # Set new search range around the closest dose for next refinement
    dose_range <- c(max(closest_row$Dose - step, 0), closest_row$Dose + step)
    
    # Store the results from this refinement step
    dose_liver_concentrations <- rbind(dose_liver_concentrations, refined_dose_liver_concentrations)
  }
  
  # Find the final closest dose after all refinements
  closest_dose <- dose_liver_concentrations[which.min(abs(dose_liver_concentrations$LiverConcentration - target_concentration)), ]
  
  return(closest_dose)
}

# Parameters for the liver dose search
dose_range <- c(1, 100)  # Single dose range to search
multiple_target_concentrations_liver <- c("targer concentrarion in liver")  # List of target concentrations for liver

# Find the closest dose for each target concentration for liver
closest_liver_dose_results <- lapply(multiple_target_concentrations_liver, function(target_concentration) {
  find_closest_liver_dose(dose_range, target_concentration)
})

# Display the results
print(closest_liver_dose_results)

#############Kidney
# Function to find the dose with kidney concentration closest to the target concentration using a refining search approach
find_closest_kidney_dose <- function(dose_range, target_concentration, refinement_steps = c(1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001,0.00000001)) {
  dose_kidney_concentrations <- data.frame(Dose = numeric(), kidneyConcentration = numeric())
  
  # Start with initial coarse step and then refine
  for (step in refinement_steps) {
    refined_dose_kidney_concentrations <- data.frame(Dose = numeric(), kidneyConcentration = numeric())
    
    for (current_dose in seq(dose_range[1], dose_range[2], by = step)) {
      # Run simulation with the current dose
      simulation_results <- pred.A(Gtheta.int, doses = c(current_dose), n_bootstrap = 100, conf_level = 0.95)
      
      # Extract kidney concentration at the last time point
      kidney_concentration <- tail(simulation_results[[as.character(current_dose)]]$simulation$CK, n = 1)
      
      # Store the dose and kidney concentration
      refined_dose_kidney_concentrations <- rbind(refined_dose_kidney_concentrations, data.frame(Dose = current_dose, kidneyConcentration = kidney_concentration))
    }
    
    # Find the dose with kidney concentration closest to the target concentration
    closest_row <- refined_dose_kidney_concentrations[which.min(abs(refined_dose_kidney_concentrations$kidneyConcentration - target_concentration)), ]
    
    # Set new search range around the closest dose for next refinement
    dose_range <- c(max(closest_row$Dose - step, 0), closest_row$Dose + step)
    
    # Store the results from this refinement step
    dose_kidney_concentrations <- rbind(dose_kidney_concentrations, refined_dose_kidney_concentrations)
  }
  
  # Find the final closest dose after all refinements
  closest_dose <- dose_kidney_concentrations[which.min(abs(dose_kidney_concentrations$kidneyConcentration - target_concentration)), ]
  
  return(closest_dose)
}

# Parameters for the kidney dose search
dose_range <- c(0, 66)  # Single dose range to search
multiple_target_concentrations_kidney <- c("targer concentrarion in kidney")  # List of target concentrations for kidney

# Find the closest dose for each target concentration for kidney
closest_kidney_dose_results <- lapply(multiple_target_concentrations_kidney, function(target_concentration) {
  find_closest_kidney_dose(dose_range, target_concentration)
})

# Display the results
print(closest_kidney_dose_results)