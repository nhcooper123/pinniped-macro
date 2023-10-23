# Get csv file of all probabilities
# relprobs_matrix_for_internal_states = pies
# MLstates = blobs
# This R-script is a modified version of the original BioGeoBEARS results
# R script that can be found at http://phylo.wikidot.com/biogeobears
# Copyright Nicholas J. Matzke

# Load libraries
library(BioGeoBEARS)
library(ggtree)
library(tidyverse)
library(ape)
library(cladoRcpp)
source("biogeography/analyses/function-extract-results.R")

# BEFORE THIS YOU NEED TO RUN 01 to get state lists
source("biogeography/analyses/01-pinniped-state-list-fix.R")

#-----------------------------------------------
# Read in the geography files and get tip ranges
#-----------------------------------------------
geogfn <- "biogeography/data/pinniped-all-geography_9areas.txt"
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
#------------------
# Read in outputs
#------------------
load("biogeography/outputs/pinnipeds-all-DEC_9areas_impossible.Rdata")
load("biogeography/outputs/pinnipeds-all-DECJ_9areas_impossible.Rdata")
#-------------------
# Read in the trees
#------------------
tree1 <- read.tree("biogeography/data/pinniped-tree-all_9areas.tre")
#----------------------------------------------
# Extract results
#----------------------------------------------
DEC_all <- extract_results(tree1, tipranges, resDEC)
DECJ_all <- extract_results(tree1, tipranges, resDECJ)
BA_all <- extract_results(tree1, tipranges, resBA)
BAJ_all <- extract_results(tree1, tipranges, resBAJ)
#-------------------------------------------
# Save all results for all nodes and states
#-------------------------------------------
DEC_allx <- DEC_all %>%
  select(node, ML, everything())

write_csv(DEC_allx, file = "supplemental/DEC-all-node-probabilities.csv")

DECJ_allx <- DECJ_all %>%
  select(node, ML, everything())

write_csv(DECJ_allx, file = "supplemental/DECJ-all-node-probabilities.csv")

BA_allx <- BA_all %>%
  select(node, ML, everything())

write_csv(BA_allx, file = "supplemental/BA-all-node-probabilities.csv")

BAJ_allx <- BAJ_all %>%
  select(node, ML, everything())

write_csv(BAJ_allx, file = "supplemental/BAJ-all-node-probabilities.csv")
