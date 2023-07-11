# Extract data for plotting from BGB
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
load("biogeography/outputs/pinnipeds-fossil-DECJ_9areas_impossible.Rdata")
load("biogeography/outputs/pinnipeds-extant-DEC_9areas_impossible.Rdata")
load("biogeography/outputs/pinnipeds-extant-DECJ_9areas_impossible.Rdata")

#-------------------
# Read in the trees
#------------------
tree1 <- read.tree("biogeography/data/pinniped-tree-all_9areas.tre")
tree_fossil <- read.tree("biogeography/data/pinniped-tree-fossil_9areas.tre")
tree_extant <- read.tree("biogeography/data/pinniped-tree-extant_9areas.tre")

#----------------------------------------------
# Extract results
#----------------------------------------------
DEC_all <- extract_results(tree1, tipranges, resDEC)
DECJ_all <- extract_results(tree1, tipranges, resDECJ)
DEC_fossil <- extract_results(tree_fossil, tipranges_fossil, resDEC_fossil)
DECJ_fossil <- extract_results(tree_fossil, tipranges_fossil, resDECJ_fossil)
DEC_extant <- extract_results(tree_extant, tipranges_extant, resDEC_extant)
DECJ_extant <- extract_results(tree_extant, tipranges_extant, resDECJ_extant)

# Get node numbers. Check using
# ggtree(tree, branch.length = "none") + geom_text2(aes(subset=!isTip, label=node), size =2,  hjust=-.3) + geom_tiplab(size = 2)

#-----------------------
# Insets - all
#-----------------------
# Choose one taxon per group
# Needs to be living if possible
focal_species <- c("Arctocephalus australis", "Phoca vitulina", 
                   "Odobenus rosmarus", "Allodesmus demerei", 
                   "Potamotherium vallentoni", "Devinophoca claytoni", 
                   "Monachus monachus")

# Replace underscores with spaces in tree file
tree1$tip.label <- gsub("_", " ", tree1$tip.label)

# Remove all other taxa from the tree
tree_basic <- drop.tip(tree1, setdiff(tree1$tip.label, focal_species))

# Rename with group names
tree_basic$tip.label <- gsub("Arctocephalus australis", "Otariidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Phoca vitulina", "Phocinae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Odobenus rosmarus", "Odobenidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Allodesmus demerei", "Desmatophocidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Potamotherium vallentoni", "stem", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Devinophoca claytoni", "Devinophocinae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Monachus monachus", "Monachinae", tree_basic$tip.label) 

basic_tree <-
  ggtree(tree_basic, branch.length = "none") %>% ggtree::rotate(8) %>% ggtree::rotate(13) %>% ggtree::rotate(12) +
  xlim(-2,10) +
  geom_tiplab(geom = "text", size = 6) +
  geom_rootedge(rootedge = 1) 
#+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab(size = 2)

#------------------------------
# Add pies
# Get node numbers. Check using
# ggtree(tree, branch.length = "none") + geom_text2(aes(subset=!isTip, label=node), size =2,  hjust=-.3) + geom_tiplab(size = 2)

DEC_all2 <- filter(DEC_all, node == 121 | node == 130|  node == 131
              | node == 132 | node == 133 | node == 186)
# Change numbers to match tree node numbers
DEC_all2$node <- 8:13 # 7 taxa, so nodes are 8:13

# Identify states with probabilities > 0.1
colnames(DEC_all2)[which(DEC_all2[1,] > 0.1)]
colnames(DEC_all2)[which(DEC_all2[2,] > 0.1)]
colnames(DEC_all2)[which(DEC_all2[3,] > 0.1)]
colnames(DEC_all2)[which(DEC_all2[4,] > 0.1)]
colnames(DEC_all2)[which(DEC_all2[5,] > 0.1)]
colnames(DEC_all2)[which(DEC_all2[6,] > 0.1)]

# Create colour palette just for these colours
colours_pies <- data.frame(state = colnames(DEC_all), col = rep("#eeeeee", length(colnames(DEC_all))))

colours_pies2 <- 
  colours_pies %>%
  dplyr::slice(-c(257,258)) %>%
  mutate(col = case_when(state == "A" ~ "#D400D4",
                         state == "AD" ~ "#8EAA39",
                         state == "D" ~ "#FFED00",
                         state == "DH" ~ "#0000ff",
                         state == "DI" ~ "#000000",
                         state == "I" ~ "#74D7EE",
                         state == "AI" ~ "#CC8899",
                         state == "ADI" ~ "#88EE33",
                         state == "ADH" ~ "#ff4d00",
                         TRUE ~ as.character(col)))

xxx <- colours_pies2 %>% arrange(state)

# Create pies
# Colour argument gives strange warning but this is just related to change in base R
# Note weird colour set up is because the state - gets removed for some reason
# but needs to be the first in the list to match dd3
pies <- nodepie(DEC_all2, cols = 1:256, alpha = 1, color = c("#eeeeee", xxx$col))

# Add pies to plot
inset(basic_tree, pies, width = 0.4, height = 1)

# To save
inset_pies_DEC_all <- inset(basic_tree, pies, width = 0.3, height = 0.8)
#------------------
# Save plot
#------------------
ggsave(inset_pies_DEC_all, file = "biogeography/outputs/biogeography-inset-all-DEC.png", 
       width = 8, height = 6, dpi = 900)

#-----------------------
# Insets - fossils
#-----------------------
# Choose one taxon per group
# Needs to be living if possible
focal_speciesf <- c("Pelagiarctos sp", "Praepusa boeska", 
                   "Odobenidae indet", "Allodesmus demerei", 
                   "Potamotherium vallentoni", "Devinophoca claytoni", 
                   "Homiphoca sp")

# Replace underscores with spaces in tree file
tree_fossil$tip.label <- gsub("_", " ", tree_fossil$tip.label)

# Remove all other taxa from the tree
tree_basicf <- drop.tip(tree_fossil, setdiff(tree_fossil$tip.label, focal_speciesf))

# Rename with group names
tree_basicf$tip.label <- gsub("Pelagiarctos sp", "Otariidae", tree_basicf$tip.label) 
tree_basicf$tip.label <- gsub("Praepusa boeska", "Phocinae", tree_basicf$tip.label) 
tree_basicf$tip.label <- gsub("Odobenidae indet", "Odobenidae", tree_basicf$tip.label) 
tree_basicf$tip.label <- gsub("Allodesmus demerei", "Desmatophocidae", tree_basicf$tip.label) 
tree_basicf$tip.label <- gsub("Potamotherium vallentoni", "stem", tree_basicf$tip.label) 
tree_basicf$tip.label <- gsub("Devinophoca claytoni", "Devinophocinae", tree_basicf$tip.label) 
tree_basicf$tip.label <- gsub("Homiphoca sp", "Monachinae", tree_basicf$tip.label) 

basic_treef <-
  ggtree(tree_basicf, branch.length = "none") %>% ggtree::rotate(8) %>% ggtree::rotate(13) %>% ggtree::rotate(12) +
  xlim(-2,10) +
  geom_tiplab(geom = "text", size = 6) +
  geom_rootedge(rootedge = 1) 
#------------------------------
# Add pies
# Get node numbers. Check using
# ggtree(tree, branch.length = "none") + geom_text2(aes(subset=!isTip, label=node), size =2,  hjust=-.3) + geom_tiplab(size = 2)

DEC_fossil2 <- filter(DEC_fossil, node == 86 | node == 95|  node == 96
                   | node == 97 | node == 98 | node == 133)
# Change numbers to match tree node numbers
DEC_fossil2$node <- 8:13 # 7 taxa, so nodes are 8:13

# Identify states with probabilities > 0.1
colnames(DEC_fossil2)[which(DEC_fossil2[1,] > 0.1)]
colnames(DEC_fossil2)[which(DEC_fossil2[2,] > 0.1)]
colnames(DEC_fossil2)[which(DEC_fossil2[3,] > 0.1)]
colnames(DEC_fossil2)[which(DEC_fossil2[4,] > 0.1)]
colnames(DEC_fossil2)[which(DEC_fossil2[5,] > 0.1)]
colnames(DEC_fossil2)[which(DEC_fossil2[6,] > 0.1)]

# Create pies
# Colour argument gives strange warning but this is just related to change in base R
# Note weird colour set up is because the state - gets removed for some reason
# but needs to be the first in the list to match dd3
pies <- nodepie(DEC_fossil2, cols = 1:256, alpha = 1, color = c("#eeeeee", xxx$col))

# Add pies to plot
inset(basic_treef, pies, width = 0.4, height = 1)

# To save
inset_pies_DEC_fossil <- inset(basic_treef, pies, width = 0.3, height = 0.8)
#------------------
# Save plot
#------------------
ggsave(inset_pies_DEC_fossil, file = "biogeography/outputs/biogeography-inset-fossil-DEC.png", 
       width = 8, height = 6, dpi = 900)

#------------------------------
# DECJ
# Add pies
# Get node numbers. Check using
# ggtree(tree, branch.length = "none") + geom_text2(aes(subset=!isTip, label=node), size =2,  hjust=-.3) + geom_tiplab(size = 2)

DECJ_fossil2 <- filter(DECJ_fossil, node == 86 | node == 95|  node == 96
                      | node == 97 | node == 98 | node == 133)
# Change numbers to match tree node numbers
DECJ_fossil2$node <- 8:13 # 7 taxa, so nodes are 8:13

# Identify states with probabilities > 0.1
colnames(DECJ_fossil2)[which(DECJ_fossil2[1,] > 0.1)]
colnames(DECJ_fossil2)[which(DECJ_fossil2[2,] > 0.1)]
colnames(DECJ_fossil2)[which(DECJ_fossil2[3,] > 0.1)]
colnames(DECJ_fossil2)[which(DECJ_fossil2[4,] > 0.1)]
colnames(DECJ_fossil2)[which(DECJ_fossil2[5,] > 0.1)]
colnames(DECJ_fossil2)[which(DECJ_fossil2[6,] > 0.1)]

# Create pies
# Colour argument gives strange warning but this is just related to change in base R
# Note weird colour set up is because the state - gets removed for some reason
# but needs to be the first in the list to match dd3
pies <- nodepie(DECJ_fossil2, cols = 1:256, alpha = 1, color = c("#eeeeee", xxx$col))

# Add pies to plot
inset(basic_treef, pies, width = 0.4, height = 1)

# To save
inset_pies_DECJ_fossil <- inset(basic_treef, pies, width = 0.3, height = 0.8)
#------------------
# Save plot
#------------------
ggsave(inset_pies_DECJ_fossil, file = "biogeography/outputs/biogeography-inset-fossil-DECJ.png", 
       width = 8, height = 6, dpi = 900)

#-----------------------
# Insets - extant
#-----------------------
# Choose one taxon per group
# Needs to be living if possible
focal_speciese <- c("Arctocephalus townsendi", "Erignathus barbatus", 
                    "Odobenus rosmarus", 
                    "Monachus monachus")

# Replace underscores with spaces in tree file
tree_extant$tip.label <- gsub("_", " ", tree_extant$tip.label)

# Remove all other taxa from the tree
tree_basice <- drop.tip(tree_extant, setdiff(tree_extant$tip.label, focal_speciese))

# Rename with group names
tree_basice$tip.label <- gsub("Arctocephalus townsendi", "Otariidae", tree_basice$tip.label) 
tree_basice$tip.label <- gsub("Erignathus barbatus", "Phocinae", tree_basice$tip.label) 
tree_basice$tip.label <- gsub("Odobenus rosmarus", "Odobenidae", tree_basice$tip.label) 
tree_basice$tip.label <- gsub("Monachus monachus", "Monachinae", tree_basice$tip.label) 

basic_treee <-
  ggtree(tree_basice, branch.length = "none") +
  xlim(-2,6) +
  geom_tiplab(geom = "text", size = 6) +
  geom_rootedge(rootedge = 1) 
#------------------------------
# Add pies
# Get node numbers. Check using
# ggtree(tree, branch.length = "none") + geom_text2(aes(subset=!isTip, label=node), size =2,  hjust=-.3) + geom_tiplab(size = 2)

DEC_extant2 <- filter(DEC_extant, node == 35 | node == 36|  node == 52)
# Change numbers to match tree node numbers
DEC_extant2$node <- 5:7 # 4 taxa, so nodes are 5:7

# Identify states with probabilities > 0.1
colnames(DEC_extant2)[which(DEC_extant2[1,] > 0.1)]
colnames(DEC_extant2)[which(DEC_extant2[2,] > 0.1)]
colnames(DEC_extant2)[which(DEC_extant2[3,] > 0.1)]

# Create pies
# Colour argument gives strange warning but this is just related to change in base R
# Note weird colour set up is because the state - gets removed for some reason
# but needs to be the first in the list to match dd3
pies <- nodepie(DEC_extant2, cols = 1:256, alpha = 1, color = c("#eeeeee", xxx$col))

# Add pies to plot
inset(basic_treee, pies, width = 0.9, height = 1)

# To save
inset_pies_DEC_extant <- inset(basic_treee, pies, width = 0.5, height = 1)
#------------------
# Save plot
#------------------
ggsave(inset_pies_DEC_extant, file = "biogeography/outputs/biogeography-inset-extant-DEC.png", 
       width = 8, height = 6, dpi = 900)

#------------------------------
# DECJ
# Add pies
# Get node numbers. Check using
# ggtree(tree, branch.length = "none") + geom_text2(aes(subset=!isTip, label=node), size =2,  hjust=-.3) + geom_tiplab(size = 2)

DECJ_extant2 <- filter(DECJ_extant, node == 35 | node == 36|  node == 52)
# Change numbers to match tree node numbers
DECJ_extant2$node <- 5:7 # 7 taxa, so nodes are 8:13

# Identify states with probabilities > 0.1
colnames(DECJ_extant2)[which(DECJ_extant2[1,] > 0.1)]
colnames(DECJ_extant2)[which(DECJ_extant2[2,] > 0.1)]
colnames(DECJ_extant2)[which(DECJ_extant2[3,] > 0.1)]

# Create pies
# Colour argument gives strange warning but this is just related to change in base R
# Note weird colour set up is because the state - gets removed for some reason
# but needs to be the first in the list to match dd3
pies <- nodepie(DECJ_extant2, cols = 1:256, alpha = 1, color = c("#eeeeee", xxx$col))

# Add pies to plot
inset(basic_treee, pies, width = 0.6, height = 1)

# To save
inset_pies_DECJ_extant <- inset(basic_treee, pies, width = 0.5, height = 1)
#------------------
# Save plot
#------------------
ggsave(inset_pies_DECJ_extant, file = "biogeography/outputs/biogeography-inset-extant-DECJ.png", 
       width = 8, height = 6, dpi = 900)


#-----------------
# Make a big plot
#----------------
# Extract inset_pies from main plotting script

# Plot all together
library(patchwork)
(inset_pies_DEC_all + inset_pies) /
(inset_pies_DEC_fossil + inset_pies_DECJ_fossil) /
(inset_pies_DEC_extant + inset_pies_DECJ_extant) + plot_annotation(tag_levels = "A")


#----------------------
# DEC all
#----------------------
a <- 
  DEC_all %>% 
  filter(node == 121) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

b <-
  DEC_all %>% 
  filter(node == 130) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

c <-
DEC_all %>% 
  filter(node == 131) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

d <-
DEC_all %>% 
  filter(node == 132) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

e <- 
DEC_all %>% 
  filter(node == 133) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_all")

f <- 
DEC_all %>% 
  filter(node == 186) %>% 
  select(-ML) %>% 
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
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

b <-
  DECJ_all %>% 
  filter(node == 130) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

c <-
  DECJ_all %>% 
  filter(node == 131) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

d <-
  DECJ_all %>% 
  filter(node == 132) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

e <- 
  DECJ_all %>% 
  filter(node == 133) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_all")

f <- 
  DECJ_all %>% 
  filter(node == 186) %>% 
  select(-ML) %>% 
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
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

b <-
  DEC_fossil %>% 
  filter(node == 95) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

c <-
  DEC_fossil %>% 
  filter(node == 96) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

d <-
  DEC_fossil %>% 
  filter(node == 97) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

e <- 
  DEC_fossil %>% 
  filter(node == 98) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_fossil")

f <- 
  DEC_fossil %>% 
  filter(node == 133) %>% 
  select(-ML) %>% 
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
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

b <-
  DECJ_fossil %>% 
  filter(node == 95) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

c <-
  DECJ_fossil %>% 
  filter(node == 96) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

d <-
  DECJ_fossil %>% 
  filter(node == 97) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

e <- 
  DECJ_fossil %>% 
  filter(node == 98) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_fossil")

f <- 
  DECJ_fossil %>% 
  filter(node == 133) %>% 
  select(-ML) %>% 
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
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_extant")

b <-
  DEC_extant %>% 
  filter(node == 36) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DEC_extant")

c <-
  DEC_extant %>% 
  filter(node == 52) %>% 
  select(-ML) %>% 
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
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_extant")

b <-
  DECJ_extant %>% 
  filter(node == 36) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_extant")

c <-
  DECJ_extant %>% 
  filter(node == 52) %>% 
  select(-ML) %>% 
  pivot_longer(!node, names_to = "state", values_to = "proportion") %>%
  arrange(-proportion) %>%
  slice(1:5) %>%
  mutate(model = "DECJ_extant")

modsDECJ_extant <- rbind(a,b,c)

#-------------------
# Stick em all together
#-------------------

all_models <- cbind(modsDEC_all, modsDECJ_all, modsDEC_fossil, 
                    modsDECJ_fossil, modsDEC_extant, modsDECJ_extant)

write.csv(all_models, file = "biogeography/outputs/all-models-best-states.csv")

#------------------------

