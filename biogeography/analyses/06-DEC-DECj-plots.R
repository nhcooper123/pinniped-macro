## BioGeoBEARS plots ###
## BTTW project July 2023

#------------------------
# Load libraries
#------------------------
library(BioGeoBEARS)
library(ape)
library(cladoRcpp)

#-----------------------------------------------------------------------
# Source script with modified BioGeoBEARs plotting function
#-----------------------------------------------------------------------
source("biogeography/analyses/modifiedBGBplot.R")

#-------------------
# Read in the trees
#------------------
tree <- read.tree("biogeography/data/pinniped-tree-all_9areas.tre")
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
load("biogeography/outputs/pinnipeds-all-DEC_9areas_neither.Rdata")
# resDEC1_fossil = fossil taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-fossil-DEC_9areas_neither.Rdata")
# resDEC1_extant = extant taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-extant-DEC_9areas_neither.Rdata")

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
load("biogeography/outputs/pinnipeds-all-DECj_9areas_neither.Rdata")
# resDECj1_fossil = fossil taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-fossil-DECj_9areas_neither.Rdata")
# resDECj1_extant = extant taxa, impossible and unlikely states removed
load("biogeography/outputs/pinnipeds-extant-DECj_9areas_neither.Rdata")

#-------------------
# Colours
#-------------------
source("biogeography/analyses/messing-about-with-colours.R")

#-----------------------------------------
# Raw trees
#-----------------------------------------
# All
#------------
png(file = "supplemental/figures/all-pinnipeds-tree.png",
    width = 5000, height = 5000, res = 300)

plot_BioGeoBEARS_mod(resDEC, titlecex = 0, plotwhat = "x", statecex = 1, 
                     tipcex = 0.8, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = TRUE)

dev.off()
#------------
# Fossil
#------------
png(file = "supplemental/figures/fossil-pinnipeds-tree.png",
    width = 5000, height = 5000, res = 300)

plot_BioGeoBEARS_mod(resDEC_fossil, titlecex = 0, plotwhat = "x", statecex = 1, 
                     tipcex = 0.8, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = TRUE)

dev.off()
#------------
# Extant
#------------
png(file = "supplemental/figures/extant-pinnipeds-tree.png",
    width = 5000, height = 5000, res = 300)

plot_BioGeoBEARS_mod(resDEC_extant, titlecex = 0, plotwhat = "x", statecex = 1, 
                     tipcex = 0.8, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = TRUE)

dev.off()

#-----------------------------------------
# Plots for appendices
#-----------------------------------------
# DEC - impossible
# 1. ML states at nodes
# 2. Pies at nodes
# DECj - impossible
# 1. ML states at nodes
# 2. Pies at nodes

# DEC - impossible and unlikely
# 1. ML states at nodes
# 2. Pies at nodes
# DECj - impossible and unlikely
# 1. ML states at nodes
# 2. Pies at nodes

# Repeat for fossils and extant only
#-----------------------------------------
# All pinnipeds
#-----------------------------------------
# DEC
# ML states at nodes
png(file = "supplemental/figures/all-pinnipeds-DEC-impossible-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC, 
                     # title info - cex = 0 suppresses title
                     titlecex = 0, 
                     # Plot ML states at nodes, and size them
                     plotwhat = "text", statecex = 1, 
                     # Size and offset of tip names
                     tipcex = 1, label.offset = 0.2, 
                     # Plot states at tips, resize tip states and offset of tip labels
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     # split splits or not?
                     plotsplits = FALSE,
                     # Include null ranges?
                     include_null_range = TRUE, 
                     # Input data - tree and geography
                     tr = tree, tipranges = tipranges,
                     # x label
                     xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     # colours
                     colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)
dev.off()
#-------------------
# DEC
# Pies at nodes
png(file = "supplemental/figures/all-pinnipeds-DEC-impossible-pies.png",
     width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#-------------------
# DEC+J
# ML states at nodes
png(file = "supplemental/figures/all-pinnipeds-DECj-impossible-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj,titlecex = 0, plotwhat = "text", statecex = 1, tipcex = 1, label.offset = 0.2, 
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, 
                     tr = tree, tipranges = tipranges, xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     colors_list_for_states = colour_list_impossible, show.tip.label = FALSE)
dev.off()
#-------------------
# DEC+J
# Pies at nodes
png(file = "supplemental/figures/all-pinnipeds-DECj-impossible-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#--------------------------
# Impossible and unlikely removed
# DEC
# ML states at nodes
png(file = "supplemental/figures/all-pinnipeds-DEC-unlikely-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC1, 
                     # title info - cex = 0 suppresses title
                     titlecex = 0, 
                     # Plot ML states at nodes, and size them
                     plotwhat = "text", statecex = 1, 
                     # Size and offset of tip names
                     tipcex = 1, label.offset = 0.2, 
                     # Plot states at tips, resize tip states and offset of tip labels
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     # split splits or not?
                     plotsplits = FALSE,
                     # Include null ranges?
                     include_null_range = TRUE, 
                     # Input data - tree and geography
                     tr = tree, tipranges = tipranges,
                     # x label
                     xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     # colours
                     colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)
dev.off()
#-------------------
# DEC
# Pies at nodes
png(file = "supplemental/figures/all-pinnipeds-DEC-unlikely-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC1, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#-------------------
# DEC+J
# ML states at nodes
png(file = "supplemental/figures/all-pinnipeds-DECj-unlikley-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj1,titlecex = 0, plotwhat = "text", statecex = 1, tipcex = 1, label.offset = 0.2, 
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, 
                     tr = tree, tipranges = tipranges, xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     colors_list_for_states = colour_list_impossible, show.tip.label = FALSE)
dev.off()
#-------------------
# DEC+J
# Pies at nodes
png(file = "supplemental/figures/all-pinnipeds-DECj-unlikely-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj1, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()

#--------------------------
# Fossil
#--------------------------
# DEC
# ML states at nodes
png(file = "supplemental/figures/fossil-pinnipeds-DEC-impossible-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC_fossil, 
                     # title info - cex = 0 suppresses title
                     titlecex = 0, 
                     # Plot ML states at nodes, and size them
                     plotwhat = "text", statecex = 1, 
                     # Size and offset of tip names
                     tipcex = 1, label.offset = 0.2, 
                     # Plot states at tips, resize tip states and offset of tip labels
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     # split splits or not?
                     plotsplits = FALSE,
                     # Include null ranges?
                     include_null_range = TRUE, 
                     # Input data - tree and geography
                     tr = tree, tipranges = tipranges,
                     # x label
                     xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     # colours
                     colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)
dev.off()
#-------------------
# DEC
# Pies at nodes
png(file = "supplemental/figures/fossil-pinnipeds-DEC-impossible-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC_fossil, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#-------------------
# DEC+J
# ML states at nodes
png(file = "supplemental/figures/fossil-pinnipeds-DECj-impossible-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj_fossil,titlecex = 0, plotwhat = "text", statecex = 1, tipcex = 1, label.offset = 0.2, 
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, 
                     tr = tree, tipranges = tipranges, xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     colors_list_for_states = colour_list_impossible, show.tip.label = FALSE)
dev.off()
#-------------------
# DEC+J
# Pies at nodes
png(file = "supplemental/figures/fossil-pinnipeds-DECj-impossible-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj_fossil, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#--------------------------
# Impossible and unlikely removed
# DEC
# ML states at nodes
png(file = "supplemental/figures/fossil-pinnipeds-DEC-unlikely-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC1_fossil, 
                     # title info - cex = 0 suppresses title
                     titlecex = 0, 
                     # Plot ML states at nodes, and size them
                     plotwhat = "text", statecex = 1, 
                     # Size and offset of tip names
                     tipcex = 1, label.offset = 0.2, 
                     # Plot states at tips, resize tip states and offset of tip labels
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     # split splits or not?
                     plotsplits = FALSE,
                     # Include null ranges?
                     include_null_range = TRUE, 
                     # Input data - tree and geography
                     tr = tree, tipranges = tipranges,
                     # x label
                     xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     # colours
                     colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)
dev.off()
#-------------------
# DEC
# Pies at nodes
png(file = "supplemental/figures/fossil-pinnipeds-DEC-unlikely-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC1_fossil, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#-------------------
# DEC+J
# ML states at nodes
png(file = "supplemental/figures/fossil-pinnipeds-DECj-unlikley-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj1_fossil,titlecex = 0, plotwhat = "text", statecex = 1, tipcex = 1, label.offset = 0.2, 
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, 
                     tr = tree, tipranges = tipranges, xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     colors_list_for_states = colour_list_impossible, show.tip.label = FALSE)
dev.off()
#-------------------
# DEC+J
# Pies at nodes
png(file = "supplemental/figures/fossil-pinnipeds-DECj-unlikely-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj1_fossil, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#------------------------
# Extant
#--------------------------
# DEC
# ML states at nodes
png(file = "supplemental/figures/extant-pinnipeds-DEC-impossible-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC_extant, 
                     # title info - cex = 0 suppresses title
                     titlecex = 0, 
                     # Plot ML states at nodes, and size them
                     plotwhat = "text", statecex = 1, 
                     # Size and offset of tip names
                     tipcex = 1, label.offset = 0.2, 
                     # Plot states at tips, resize tip states and offset of tip labels
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     # split splits or not?
                     plotsplits = FALSE,
                     # Include null ranges?
                     include_null_range = TRUE, 
                     # Input data - tree and geography
                     tr = tree, tipranges = tipranges,
                     # x label
                     xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     # colours
                     colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)
dev.off()
#-------------------
# DEC
# Pies at nodes
png(file = "supplemental/figures/extant-pinnipeds-DEC-impossible-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC_extant, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#-------------------
# DEC+J
# ML states at nodes
png(file = "supplemental/figures/extant-pinnipeds-DECj-impossible-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj_extant,titlecex = 0, plotwhat = "text", statecex = 1, tipcex = 1, label.offset = 0.2, 
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, 
                     tr = tree, tipranges = tipranges, xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     colors_list_for_states = colour_list_impossible, show.tip.label = FALSE)
dev.off()
#-------------------
# DEC+J
# Pies at nodes
png(file = "supplemental/figures/extant-pinnipeds-DECj-impossible-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj_extant, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#--------------------------
# Impossible and unlikely removed
# DEC
# ML states at nodes
png(file = "supplemental/figures/extant-pinnipeds-DEC-unlikely-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC1_extant, 
                     # title info - cex = 0 suppresses title
                     titlecex = 0, 
                     # Plot ML states at nodes, and size them
                     plotwhat = "text", statecex = 1, 
                     # Size and offset of tip names
                     tipcex = 1, label.offset = 0.2, 
                     # Plot states at tips, resize tip states and offset of tip labels
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     # split splits or not?
                     plotsplits = FALSE,
                     # Include null ranges?
                     include_null_range = TRUE, 
                     # Input data - tree and geography
                     tr = tree, tipranges = tipranges,
                     # x label
                     xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     # colours
                     colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)
dev.off()
#-------------------
# DEC
# Pies at nodes
png(file = "supplemental/figures/extant-pinnipeds-DEC-unlikely-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDEC1_extant, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
#-------------------
# DEC+J
# ML states at nodes
png(file = "supplemental/figures/extant-pinnipeds-DECj-unlikley-MLstates.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj1_extant,titlecex = 0, plotwhat = "text", statecex = 1, tipcex = 1, label.offset = 0.2, 
                     tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, 
                     tr = tree, tipranges = tipranges, xlab = "Time (Ma)", simplify_piecharts = TRUE, 
                     colors_list_for_states = colour_list_impossible, show.tip.label = FALSE)
dev.off()
#-------------------
# DEC+J
# Pies at nodes
png(file = "supplemental/figures/extant-pinnipeds-DECj-unlikely-pies.png",
    width = 6000, height = 6000, res = 300)

plot_BioGeoBEARS_mod(resDECj1_extant, titlecex = 0, plotwhat = "pie", statecex = 1, 
                     tipcex = 1, label.offset = 0.2, tipboxes_TF = TRUE, tip_stateshape = 22, tip_statecex = 3, tip_stateadj = 0.5,
                     plotsplits = FALSE, include_null_range = TRUE, tr = tree, tipranges = tipranges,
                     xlab = "Time (Ma)", simplify_piecharts = FALSE, colors_list_for_states = colour_list_impossible,
                     show.tip.label = FALSE)

dev.off()
