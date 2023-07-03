# Pretty BGB results plots
# Make trees using ggtree.
#------------------------
# Load libraries
#------------------------
library(BioGeoBEARS)
library(geiger)
library(ape)
library(tidyverse)
library(ggtree)
library(ggimage)
library(patchwork)

#-------------------
# Read in the tree
#------------------
tree <- read.tree("biogeography/data/pinniped-tree-all_9areas.tre")
# Replace underscores with spaces in tree file
tree$tip.label <- gsub("_", " ", tree$tip.label)
# Reorder to match BGB code
tree <- reorder(tree, "pruningwise")

#-----------------
# Read in data
#-----------------
areas <- read_PHYLIP_data("biogeography/data/pinniped-all-geography_9areas.txt")
# Replace underscores with spaces
rownames(areas) <- gsub("_", " ", rownames(areas))

# To get different colours for each area 
# we need to replace 1 with numbers 2-9 IN B-I
areas2 <-  mutate(areas, B = case_when(B == "1" ~ "2", TRUE ~ as.character(B)),
                       C = case_when(C == "1" ~ "3", TRUE ~ as.character(C)),
                       D = case_when(D == "1" ~ "4", TRUE ~ as.character(D)),
                       E = case_when(E == "1" ~ "5", TRUE ~ as.character(E)),
                       F = case_when(F == "1" ~ "6", TRUE ~ as.character(F)),
                       G = case_when(G == "1" ~ "7", TRUE ~ as.character(G)),
                       H = case_when(H == "1" ~ "8", TRUE ~ as.character(H)),
                       I = case_when(I == "1" ~ "9", TRUE ~ as.character(I)))
#-----------------
# Define colours
#-----------------
basic_cols <- c("#eeeeee", "#D400D4", "#24408E", "#008026","#FFED00",
                "#FF8C00", "#E40303", "#613915", "#FFAFC8", "#74D7EE")

#----------------------------
# Identify MRCA of each group
# -----------------------------
# Read in file with taxonomy info
groups <- read.csv("data/taxa_groups.csv")
# Replace underscores with spaces
groups$Taxon <- gsub("_", " ", groups$Taxon)

# Match to tree
check <- name.check(tree, groups, data.names = groups$Taxon)

# Remove species not in the tree
matches <- match(groups$Taxon, check$data_not_tree, nomatch = 0)
groups <- subset(groups, matches == 0)

# Get node numbers
# Note that ggtree does something odd to node numbering meaning these don't all
# match 100%. Check using
# ggtree(tree) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()
phocid <- 116
otarid <- 177
walrus <- 155
desmo <- 201

# Make df of numbers and names
df <- data.frame(node = c(phocid, otarid, walrus, desmo),
                 name = c("Phocidae", "Otariidae", "Odobenidae", "Desmatophocidae"))

#-----------------------------
# Plot tree with areas at tips
#-----------------------------
# Make the tree base
base <- 
  ggtree(tree) +
  xlim(0,100) +
  geom_tiplab(geom = "text", fontface = "italic", size = 2)

#base
#------------------
# Add areas
# Ignore warning message about scales
area_plot <-
  gheatmap(base, areas2, offset = 10, width = 0.2,
         font.size = 2, colnames_position = "top", color = "black")+
         scale_fill_manual(values = basic_cols) +
  theme(legend.position = "none")
  
#area_plot  
#------------------
# Add clade labels
area_group_plot <- 
  area_plot + geom_cladelab(data = df, mapping = aes(node = node, label = name),
              offset = c(20,20,20,30), offset.text = 1)
  
#-------------------
# Add pies
# THIS DOESN'T LOOK GREAT BUT LEAVING HERE FOR POSSIBLE FUTURE USE
# To get dd2 need to source this...
source("biogeography/analyses/07-Extract-BGB-results-for-plotting.R")

# Fix node numbers (there are 104 taxa)
# dd2$node <- dd2$node + 104
# pies <- nodepie(dd2, cols = 1:256, alpha = 0.6, color = colour_list_impossible)
# This takes a while to run
i#nset(area_group_plot, pies, width = 0.1, height = 1)
#-------------------
# Add ML states
-------------------------------------------------
# Need to first fix the colours so they are in alphabetical order
# Take the legend code to make dff
dff <- data.frame(state = c("A", "B", "C", "D", "E", "F", "G", "H", "I", 
  "AB", "AC", "AD", "AH", "BC", "BE", "CF", "DG", "DH", "DI", "EI",
  "ABC", "ADH", "BCF", "BCG", "BDE", "CDI", "CEF", "EFG", "ABCG"), 
col = c("#D400D4", "#24408E", "#008026","#FFED00", "#FF8C00", "#E40303", "#613915", "#FFAFC8", "#74D7EE",
        "#732950","#4B0082", "#0000ff", "#00FFAA", "#FFED80", "#FF5500", "#FF3990", "#FFF333","#732982","#00AAFF","#AA718E",
        "#558080","#8EAA39", "#6A40AA","#2A6AEA", "#235347", "#6A6AAA", "#FFCC00","#FFA500","#800000" ))

# Now run this to see what states require colours
# base + geom_nodepoint(aes(colour = MLstates), size = 2) 

# Create list of colours, remember to remove non needed states
new_colours <- dff %>%
  arrange(state) %>%
  filter(state != "AH" & state != "BCG", state != "CEF", state != "CF", state != "EFG", state != "F") %>%
  pull(col)

# Plot 
area_group_ML_plot <- area_group_plot + geom_nodepoint(aes(colour = MLstates), size = 2) +
  scale_colour_manual(values = new_colours) +
  theme(legend.position = "none")

# replot withbase tree and legend to check colours match properly
# base + geom_nodepoint(aes(colour = MLstates), size = 2) +
#    scale_colour_manual(values = new_colours)

#-----------------------
# Inset
#-----------------------
# Choose one taxon per group
# Needs to be living if possible
focal_species <- c("Arctocephalus australis", "Phoca vitulina", 
                   "Odobenus rosmarus", "Allodesmus demerei", 
                   "Potamotherium vallentoni")

# Remove all other taxa from the tree
tree_basic <- drop.tip(tree, setdiff(tree$tip.label, focal_species))

# Rename with group names
tree_basic$tip.label <- gsub("Arctocephalus australis", "Otariidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Phoca vitulina", "Phocidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Odobenus rosmarus", "Odobenidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Allodesmus demerei", "Desmatophocidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Potamotherium vallentoni", "stem", tree_basic$tip.label) 

basic_tree <-
  ggtree(tree_basic) +
  xlim(-3,40) +
  geom_tiplab(geom = "text", size = 4) +
  geom_rootedge(rootedge = 3)


area_group_ML_plot + (plot_spacer()/basic_tree)

#------------------
# Save plot
#------------------
ggsave(area_group_ML_plot, file = "biogeography/outputs/biogeography-nice-figure.png", 
       width = 10, height = 10, dpi = 900)

