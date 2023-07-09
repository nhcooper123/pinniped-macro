### BioGeoBEARS analysis ###
## BTTW project July 2023
## FOSSILS only. 
## Remove impossible and unlikely states

# This R-script is a modified version of the original BioGeoBEARS
# R script that can be found at http://phylo.wikidot.com/biogeobears
# Copyright Nicholas J. Matzke

# BEFORE THIS YOU NEED TO RUN 01 to get state lists
source("biogeography/analyses/01-pinniped-state-list-fix.R")

# Load the packages.
library(cladoRcpp)
library(BioGeoBEARS)

#-----------------------------------------------------
# Locate the tree and geography files
#------------------------------------------------------
trfn <- "biogeography/data/pinniped-tree-fossil_9areas.tre"
geogfn <- "biogeography/data/pinniped-fossil-geography_9areas.txt"

#----------------------------------
# Set tip ranges and max range size
#----------------------------------
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
# Set the maximum number of areas any species may occupy (max observed + 1 = 4)
max_range_size <- 4

#---------------------------------
# DEC ANALYSIS SETUP
#---------------------------------
# Initialize a default model (DEC model)
BioGeoBEARS_run_object <- define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn <- trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn <- geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size <- max_range_size

BioGeoBEARS_run_object$min_branchlength <- 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range <- TRUE    

# Set up bins for time-stratified analyses
BioGeoBEARS_run_object$timesfn <- "biogeography/data/timeperiods-all_9areas.txt"

# Speed options and multicore processing if desired
# returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$on_NaN_error <- -1e50    
# shortcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$speedup <- TRUE          
BioGeoBEARS_run_object$use_optimx <- "GenSA"
BioGeoBEARS_run_object$num_cores_to_use <- 1
# force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object$force_sparse <- FALSE    

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object <- section_the_tree(inputs = BioGeoBEARS_run_object, 
                                           make_master_table = TRUE, 
                                           plot_pieces = FALSE, 
                                           cut_fossils = FALSE)
# The stratified tree is described in this table:
BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run

# Setup state lists
# Impossible and unlikely removed
BioGeoBEARS_run_object$lists_of_states_lists_0based[[1]] = states_list_0based_1B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[2]] = states_list_0based_2B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[3]] = states_list_0based_3B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[4]] = states_list_0based_4B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[5]] = states_list_0based_5B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[6]] = states_list_0based_6B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[7]] = states_list_0based_7B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[8]] = states_list_0based_8B

# Check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#----------------------------------
# Run DEC model and save results
#----------------------------------
# Impossible and unlikely removed
resfn1_fossil <- "biogeography/outputs/pinnipeds-fossil-DEC_9areas_neither.Rdata"
resDEC1_fossil <- bears_optim_run(BioGeoBEARS_run_object)
resDEC1_fossil    
save(resDEC1_fossil, file = resfn1_fossil)

#----------------------------------
# SETUP DEC+J
#----------------------------------
# Load DEC results
load("biogeography/outputs/pinnipeds-fossil-DEC_9areas_unlikely.Rdata")

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE    

# time-stratified analysis:
BioGeoBEARS_run_object$timesfn <- "biogeography/data/timeperiods-all_9areas.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "time_strat_disp.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error <- -1e50    
BioGeoBEARS_run_object$speedup <- TRUE          
BioGeoBEARS_run_object$use_optimx <- "GenSA"    
BioGeoBEARS_run_object$num_cores_to_use <- 1
BioGeoBEARS_run_object$force_sparse <- FALSE    

# Read in files and check
BioGeoBEARS_run_object <- readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata 
BioGeoBEARS_run_object <- section_the_tree(inputs = BioGeoBEARS_run_object, 
                                           make_master_table = TRUE, 
                                           plot_pieces = FALSE, 
                                           cut_fossils = FALSE)

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart <- resDEC1_fossil$outputs@params_table["d","est"]
estart <- resDEC1_fossil$outputs@params_table["e","est"]
jstart <- 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Setup state lists
# Impossible and unlikely removed
BioGeoBEARS_run_object$lists_of_states_lists_0based[[1]] = states_list_0based_1B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[2]] = states_list_0based_2B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[3]] = states_list_0based_3B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[4]] = states_list_0based_4B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[5]] = states_list_0based_5B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[6]] = states_list_0based_6B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[7]] = states_list_0based_7B
BioGeoBEARS_run_object$lists_of_states_lists_0based[[8]] = states_list_0based_8B

# Check object
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#----------------------------------
# Run DECJ model and save results
#----------------------------------
# Impossible and unlikely removed
resfnDECJ1_fossil <- "biogeography/outputs/pinnipeds-fossil-DECJ_9areas_neither.Rdata"
resDECJ1_fossil <- bears_optim_run(BioGeoBEARS_run_object)
resDECJ1_fossil    
save(resDECJ1_fossil, file = resfnDECJ1_fossil)