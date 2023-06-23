## BioGeoBEARS results ###
## BTTW project July 2023

#------------------------
# Load libraries
#------------------------
library(BioGeoBEARS)
library(tidyverse)
library(ape)

# Quick AIC function
AIC <- function(loglik, nparams){
  2*nparams - 2*loglik
}

#---------------------------------------
# Load results from BioGeoBEARS analyses
#---------------------------------------
# DEC
#-------
# resDEC = all taxa, impossible states removed
load("biogeography/outputs/pinnipeds-all-DEC_9areas_impossible.Rdata")
# resDEC_fossil = fossil taxa, impossible states removed
load("biogeography/outputs/pinnipeds-fossil-DEC_9areas_impossible.Rdata")
# resDEC_extant = extant taxa, impossible states removed
load("biogeography/outputs/pinnipeds-extant-DEC_9areas_impossible.Rdata")

# resDEC1 = all taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-all-DEC_9areas_unlikely.Rdata")
# resDEC1_fossil = fossil taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-fossil-DEC_9areas_unlikely.Rdata")
# resDEC1_extant = extant taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-extant-DEC_9areas_unlikely.Rdata")

#-------
# DECj
#-------
# resDECj = all taxa, impossible states removed
load("biogeography/outputs/pinnipeds-all-DECj_9areas_impossible.Rdata")
# resDECj_fossil = fossil taxa, impossible states removed
load("biogeography/outputs/pinnipeds-fossil-DECj_9areas_impossible.Rdata")
# resDECj_extant = extant taxa, impossible states removed
load("biogeography/outputs/pinnipeds-extant-DECj_9areas_impossible.Rdata")

# resDECj1 = all taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-all-DECj_9areas_unlikely.Rdata")
# resDECj1_fossil = fossil taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-fossil-DECj_9areas_unlikely.Rdata")
# resDECj1_extant = extant taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-extant-DECj_9areas_unlikely.Rdata")

#----------------------------
# Summary stats DEC, DEC+J
#----------------------------
# Extract parameters
#----------------------------
# All
params_DEC <- extract_params_from_BioGeoBEARS_results_object(results_object = resDEC, returnwhat = "table", 
                                                          paramsstr_digits = 4)
params_DECj <- extract_params_from_BioGeoBEARS_results_object(results_object = resDECj, returnwhat = "table", 
                                                           addl_params=c("j"), paramsstr_digits = 4)
params_DEC1 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDEC1, returnwhat = "table", 
                                                           paramsstr_digits = 4)
params_DECj1 <- extract_params_from_BioGeoBEARS_results_object(results_object = resDECj1, returnwhat = "table", 
                                                               addl_params=c("j"), paramsstr_digits = 4)
# Fossil
params_DEC_fossil <- extract_params_from_BioGeoBEARS_results_object(results_object = resDEC_fossil, returnwhat = "table", 
                                                                 paramsstr_digits = 4)
params_DECj_fossil <- extract_params_from_BioGeoBEARS_results_object(results_object = resDECj_fossil, returnwhat = "table", 
                                                                     addl_params=c("j"), paramsstr_digits = 4)
params_DEC1_fossil <- extract_params_from_BioGeoBEARS_results_object(results_object = resDEC1_fossil, returnwhat = "table", 
                                                                  paramsstr_digits = 4)
params_DECj1_fossil <- extract_params_from_BioGeoBEARS_results_object(results_object = resDECj1_fossil, returnwhat = "table", 
                                                                      addl_params=c("j"), paramsstr_digits = 4)
# Extant
params_DEC_extant <- extract_params_from_BioGeoBEARS_results_object(results_object = resDEC_extant, returnwhat = "table", 
                                                                 paramsstr_digits = 4)
params_DECj_extant <- extract_params_from_BioGeoBEARS_results_object(results_object = resDECJ_extant, returnwhat = "table", 
                                                                     addl_params=c("j"), paramsstr_digits = 4)
params_DEC1_extant <- extract_params_from_BioGeoBEARS_results_object(results_object = resDEC1_extant, returnwhat = "table", 
                                                                  paramsstr_digits = 4)
params_DECj1_extant <- extract_params_from_BioGeoBEARS_results_object(results_object = resDECJ1_extant, returnwhat = "table", 
                                                                      addl_params=c("j"), paramsstr_digits = 4)

#-----------------------------------------
# Stick results tables together
# Need to add J column to DEC models first
#-----------------------------------------
results <-
  rbind(cbind(params_DECx, data.frame(j = c(0))), params_DECj, 
        cbind(params_DEC1, data.frame(j = c(0))), params_DECj1,
        cbind(params_DEC_fossil, data.frame(j = c(0))), params_DECj_fossil, 
        cbind(params_DEC1_fossil, data.frame(j = c(0))), params_DECj1_fossil,
        cbind(params_DEC_extant, data.frame(j = c(0))), params_DECj_extant, 
        cbind(params_DEC1_extant, data.frame(j = c(0))), params_DECj1_extant)

# Add other useful info...

results2 <- results %>%
  # Add model names
  mutate(model = rep(c("DEC", "DEC+J"), 6)) %>%
  # Add states
  mutate(states = rep(c(rep("impossible", 2), rep("unlikely" , 2)), 3)) %>%
  # Add taxa
  mutate(taxa = c(rep("all", 4), rep("fossil" , 4), rep("extant" , 4))) %>%
  # Add AIC
  mutate(AIC = AIC(results$LnL, results$numparams)) %>%
  # Reorder
  dplyr::select(taxa, states, model, LnL, AIC, d, e, j)

# Save outputs
write_csv(results2, file = "biogeography/outputs/results-DEC-DECJ_9areas.csv")