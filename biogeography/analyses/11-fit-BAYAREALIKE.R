### BioGeoBEARS analysis ###
## BTTW project July 2023
## ALL TAXA. 
### NOTE THIS CODE DOES NOT WORK ###
### FAILS EITHER TO FIT DUE TO OVER CONSTRAINED ANALYSIS IF TIME PERIODS ARE INCLUDED
### OR FAILS DUE TO MEMORY ALLOCATION MAX BEING REACHED IF ANALYSIS IS NOT CONSTRAINED

# This R-script is a modified version of the original BioGeoBEARS
# R script that can be found at http://phylo.wikidot.com/biogeobears
# Copyright Nicholas J. Matzke

# Load the packages.
library(cladoRcpp)
library(BioGeoBEARS)

#-----------------------------------------------------
# Locate the tree and geography files
#------------------------------------------------------
trfn <- "biogeography/data/pinniped-tree-all_9areas.tre"
geogfn <- "biogeography/data/pinniped-all-geography_26areas.txt"

#----------------------------------
# Set tip ranges and max range size
#----------------------------------
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
# Set the maximum number of areas any species may occupy (max observed + 1 = 9)
max_range_size <- 9

#---------------------------------
# BAYAREALIKE ANALYSIS SETUP
#---------------------------------
# Initialize a default model
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
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs; fixing any initial ("init") values outside min/max
BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)

# Check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#----------------------------------------
# Run BAYAREALIKE model and save results
#----------------------------------------
resfn <- "biogeography/outputs/pinnipeds-all-BAYAREALIKE_26areas_.Rdata"
resBA <- bears_optim_run(BioGeoBEARS_run_object)
resBA    
save(resBA, file = resfn)

#----------------------------------
# SETUP BAYAREALIKE+J
#----------------------------------
# Load BAYAREALIKE results
load("biogeography/outputs/pinnipeds-all-BAYAREALIKE_26areas.Rdata")

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$min_branchlength <- 0.000001
BioGeoBEARS_run_object$include_null_range <- TRUE    

# time-stratified analysis:
BioGeoBEARS_run_object$timesfn <- "biogeography/data/timeperiods-all_9areas.txt"

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

# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart <- resBA$outputs@params_table["d","est"]
estart <- resBA$outputs@params_table["e","est"]
jstart <- 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check object
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

#----------------------------------------
# Run BAYAREALIKE model and save results
#----------------------------------------
# Impossible removed only
resfnBAJ <- "biogeography/outputs/pinnipeds-all-BAJ_26areas.Rdata"
resBAJ <- bears_optim_run(BioGeoBEARS_run_object)
resBAJ    
save(resBAJ, file = resfnBAJ)

#----------------------------------
# Extract AIC
#----------------------------------
# We have to extract the log-likelihood differently, depending on the 
# version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBA)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAJ)
numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)
stats