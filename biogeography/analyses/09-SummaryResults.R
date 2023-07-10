# Extract data for plotting from BGB
# relprobs_matrix_for_internal_states = pies
# MLstates = blobs
# This R-script is a modified version of the original BioGeoBEARS results
# R script that can be found at http://phylo.wikidot.com/biogeobears
# Copyright Nicholas J. Matzke

# Load libraries
library(BioGeoBEARS)
library(ape)
library(cladoRcpp)

# BEFORE THIS YOU NEED TO RUN 01 to get state lists
source("biogeography/analyses/01-pinniped-state-list-fix.R")

#-------------------
# Read in the trees
#------------------
tree1 <- read.tree("biogeography/data/pinniped-tree-all_9areas.tre")
tree_fossil <- read.tree("biogeography/data/pinniped-tree-fossil_9areas.tre")
tree_extant <- read.tree("biogeography/data/pinniped-tree-extant_9areas.tre")

#-----------------------------------------------
# Read in the geography files and get tip ranges
#-----------------------------------------------
geogfn <- "biogeography/data/pinniped-all-geography_9areas.txt"
tipranges <- getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)

# Fossil
geogfn_fossil <- "biogeography/data/pinniped-fossil-geography_9areas.txt"
tipranges_fossil <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn_fossil)

# Extant
geogfn_extant <- "biogeography/data/pinniped-extant-geography_9areas.txt"
tipranges_extant <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn_extant)

#------------------
# Read in outputs
#------------------
load("biogeography/outputs/pinnipeds-all-DEC_9areas_impossible.Rdata")
load("biogeography/outputs/pinnipeds-all-DECJ_9areas_impossible.Rdata")
load("biogeography/outputs/pinnipeds-fossil-DEC_9areas_impossible.Rdata")
load("biogeography/outputs/pinnipeds-all-DECJ_9areas_impossible.Rdata")
load("biogeography/outputs/pinnipeds-extant-DEC_9areas_impossible.Rdata")
load("biogeography/outputs/pinnipeds-extant-DECJ_9areas_impossible.Rdata")

#-----------------------------------------------
# Get colours 
#-----------------------------------------------
source("biogeography/analyses/messing-about-with-colours.R")
colors_list_for_states <- colour_list_impossible


#-----------------------------------------------
# Get area and state lists
#-----------------------------------------------
get_results <- function(tree, tipranges, results_object){

tr_pruningwise <- reorder(tree, "pruningwise")
tips <- 1:length(tr_pruningwise$tip.label)
nodes <- (length(tr_pruningwise$tip.label) + 1):(length(tr_pruningwise$tip.label) + 
                                                   tr_pruningwise$Nnode)

areas <- getareas_from_tipranges_object(tipranges)
numareas <- length(areas)
max_range_size <- 4

numstates <- numstates_from_numareas(numareas = length(areas), 
                                     maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)

states_list_areaLetters <- areas_list_to_states_list_new(areas, 
                                                         maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)

states_list_0based_index <- rcpp_areas_list_to_states_list(areas, 
                                                           maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
#-----------------------------------------------
# Set up matrix for node values
#-----------------------------------------------
leftright_nodes_matrix <- get_leftright_nodes_matrix_from_results(tr_pruningwise)

# Get marginal probs
marprobs <- results_object$ML_marginal_prob_each_state_at_branch_bottom_below_node

left_ML_marginals_by_node <- marprobs[leftright_nodes_matrix[, 2], ]
right_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 1], ]

#-----------------------------------------------
# Extract relative probabilities
#-----------------------------------------------
relprobs_matrix <- results_object$ML_marginal_prob_each_state_at_branch_top_AT_node

relprobs_matrix_for_internal_states <- relprobs_matrix[nodes, ]

#-----------------------------------------------
# Get list of ranges
#-----------------------------------------------
ranges_list <- states_list_0based_to_ranges_txt_list(state_indices_0based = states_list_0based_index, 
                                                     areanames = areas)
statenames <- unlist(ranges_list)

#-----------------------------------------------
# Extract ML probabilities at nodes
#-----------------------------------------------
MLprobs <- get_ML_probs(relprobs_matrix)
MLstates <- get_ML_states_from_relprobs(relprobs_matrix, statenames, 
                                        returnwhat = "states", if_ties = "takefirst")

#-----------------------------------------------
# Sort colours
#-----------------------------------------------
possible_ranges_list_txt <- ranges_list    

cols_byNode <- rangestxt_to_colors(possible_ranges_list_txt, 
                                   colors_list_for_states, MLstates)

#----------------------------------------------
# Set up dataframe for plotting
#----------------------------------------------
dd2 <- data.frame(relprobs_matrix_for_internal_states)
colnames(dd2) <- ranges_list
dd2$node <- (length(tree$tip.label)+1):((length(tree$edge)/2)+1) 
dd2$ML <- MLstates[(length(tree$tip.label)+1):((length(tree$edge)/2)+1) ] # First node is 121
dd2$colour <- cols_byNode[(length(tree$tip.label)+1):((length(tree$edge)/2)+1) ]

return(dd2)
}
#----------------------------------------------
# Run function
#----------------------------------------------
DEC_all <- get_results(tree1, tipranges, resDEC)
DECJ_all <- get_results(tree1, tipranges, resDECJ)
DEC_fossil <- get_results(tree_fossil, tipranges_fossil, resDEC_fossil)
DECJ_fossil <- get_results(tree_fossil, tipranges_fossil, resDECJ_fossil)
DEC_extant <- get_results(tree_extant, tipranges_extant, resDEC_extant)
DECJ_extant <- get_results(tree_extant, tipranges_extant, resDECJ_extant)


# Get node numbers. Check using
# ggtree(tree, branch.length = "none") + geom_text2(aes(subset=!isTip, label=node), size =2,  hjust=-.3) + geom_tiplab(size = 2)

dd3 <- filter(dd2, node == 121 | node == 130|  node == 131
              | node == 132 | node == 133 | node == 186)

colnames(DEC_all)[which(DEC_all[121,] > 0.1)]

#----------------------
# DEC all
#----------------------
a <- 
  DEC_all %>% 
  filter(node == 121) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

b <-
  DEC_all %>% 
  filter(node == 130) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

c <-
DEC_all %>% 
  filter(node == 131) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

d <-
DEC_all %>% 
  filter(node == 132) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

e <- 
DEC_all %>% 
  filter(node == 133) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

f <- 
DEC_all %>% 
  filter(node == 186) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

modsDEC_all <- rbind(a,b,c,d,e,f)
#----------------------
# DECJ all
#----------------------
a <- 
  DECJ_all %>% 
  filter(node == 121) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

b <-
  DECJ_all %>% 
  filter(node == 130) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

c <-
  DECJ_all %>% 
  filter(node == 131) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

d <-
  DECJ_all %>% 
  filter(node == 132) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

e <- 
  DECJ_all %>% 
  filter(node == 133) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

f <- 
  DECJ_all %>% 
  filter(node == 186) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

modsDECJ_all <- rbind(a,b,c,d,e,f)

models_all <- rbind(modsDEC_all, modsDECJ_all)
#----------------------
# DEC fossil
#----------------------
a <- 
  DEC_fossil %>% 
  filter(node == 86) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

b <-
  DEC_fossil %>% 
  filter(node == 95) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

c <-
  DEC_fossil %>% 
  filter(node == 96) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

d <-
  DEC_fossil %>% 
  filter(node == 97) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

e <- 
  DEC_fossil %>% 
  filter(node == 98) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

f <- 
  DEC_fossil %>% 
  filter(node == 133) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

modsDEC_fossil <- rbind(a,b,c,d,e,f)
#----------------------
# DECJ fossil
#----------------------
a <- 
  DECJ_fossil %>% 
  filter(node == 86) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

b <-
  DECJ_fossil %>% 
  filter(node == 95) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

c <-
  DECJ_fossil %>% 
  filter(node == 96) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

d <-
  DECJ_fossil %>% 
  filter(node == 97) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

e <- 
  DECJ_fossil %>% 
  filter(node == 98) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

f <- 
  DECJ_fossil %>% 
  filter(node == 133) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

modsDECJ_fossil <- rbind(a,b,c,d,e,f)

#----------------------
# DEC extant
#----------------------
a <- 
  DEC_extant %>% 
  filter(node == 35) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_extant")

b <-
  DEC_extant %>% 
  filter(node == 36) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_extant")

c <-
  DEC_extant %>% 
  filter(node == 52) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_extant")

modsDEC_extant <- rbind(a,b,c)
#----------------------
# DECJ extant
#----------------------
a <- 
  DECJ_extant %>% 
  filter(node == 35) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_extant")

b <-
  DECJ_extant %>% 
  filter(node == 36) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_extant")

c <-
  DECJ_extant %>% 
  filter(node == 52) %>% 
  select(-ML, -colour) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_extant")

modsDECJ_extant <- rbind(a,b,c)

#-------------------
# Stick em all together
#-------------------

all_models <- rbind(modsDEC_all, modsDECJ_all, modsDEC_fossil, 
                    modsDECJ_fossil, modsDEC_extant, modsDECJ_extant)

write.csv(all_models, file = "biogeography/outputs/all-models-best-states.csv")
