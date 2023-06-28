## BioGeoBEARS plots ###
## BTTW project July 2023

#------------------------
# Load libraries
#------------------------
library(BioGeoBEARS)
library(ape)
library(cladoRcpp)

# Plotting scripts
# scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))
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

#-----------------------------------------------------------------------
# Source script with modified BioGeoBEARs plotting function and colours
#-----------------------------------------------------------------------
source("biogeography/analyses/modifiedBGBplot.R")
source("biogeography/analyses/messing-about-with-colours.R")

#--------------------
# Main results plot
#--------------------
pdf("biogeography/outputs/faff.pdf", width = 20, height = 20)
plot_BioGeoBEARS_mod(resDEC, 
                     # title info - cex = 0 suppresses title
                     analysis_titletxt = "DEC all", titlecex = 0, 
                     # Plot pies at nodes, and size them
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

#--------------------------
# Make plots for appendices
#--------------------------
# All
pdf(file = "biogeography/outputs/all-pinnipeds-DEC-DECj-plots.pdf", height = 20)

dev.off()
#------------------------
# Fossil
pdf(file = "biogeography/outputs/fossil-pinnipeds-DEC-DECj-plots.pdf", height = 20)


dev.off()

#------------------------
# Extant
pdf(file = "biogeography/outputs/extant-pinnipeds-DEC-DECj-plots.pdf", height = 20)

dev.off()
