## BioGeoBEARS plots ###
## BTTW project July 2023

#------------------------
# Load libraries
#------------------------
library(BioGeoBEARS)
library(ape)

# Quick AIC function
AIC <- function(loglik, nparams){
  2*nparams - 2*loglik
}

# Plotting scripts
scriptdir <- np(system.file("extdata/a_scripts", package = "BioGeoBEARS"))
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

#------------------------
# Make plots
#------------------------
# All
pdf(file = "biogeography/outputs/all-pinnipeds-DEC-DECj-plots.pdf", height = 20)
plot_BioGeoBEARS_results(resDEC, "DEC all", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

plot_BioGeoBEARS_results(resDECj, "DECj all", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

plot_BioGeoBEARS_results(resDEC1, "DEC unlikely all", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

plot_BioGeoBEARS_results(resDECj1, "DECj unlikely all", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

dev.off()
#------------------------
# Fossil
pdf(file = "biogeography/outputs/fossil-pinnipeds-DEC-DECj-plots.pdf", height = 20)
plot_BioGeoBEARS_results(resDEC_fossil, "DEC fossil", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

plot_BioGeoBEARS_results(resDECj_fossil, "DECj fossil", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

plot_BioGeoBEARS_results(resDEC1_fossil, "DEC unlikely fossil", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

plot_BioGeoBEARS_results(resDECj1_fossil, "DECj unlikely fossil", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

dev.off()

#------------------------
# Extant
pdf(file = "biogeography/outputs/extant-pinnipeds-DEC-DECj-plots.pdf", height = 20)
plot_BioGeoBEARS_results(resDEC_extant, "DEC extant", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

plot_BioGeoBEARS_results(resDECj_extant, "DECj extant", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

plot_BioGeoBEARS_results(resDEC1_extant, "DEC unlikely extant", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

plot_BioGeoBEARS_results(resDECj1_fextant, "DECj unlikely extant", plotwhat = "pie", 
                         label.offset = 0.45, tipcex = 0.3, statecex = 0.7, splitcex = 0.6, 
                         titlecex = 0, plotsplits = TRUE, cornercoords_loc = scriptdir, 
                         include_null_range = TRUE, tr = tree, tipranges = tipranges,
                         tipboxes_TF = FALSE, xlab = "Time (Ma)")

dev.off()


#





#colors_list_for_states = c("#e3342f",
                           "#f6993f",
                           "#ffed4a",
                           "#38c172",
                           "#4dc0b5",
                           "#3490dc",
                           "#6574cd",
                           "#9561e2",
                           "#f66d9b")

# options
results_object, analysis_titletxt = NULL, addl_params = list(), 
plotwhat = "text", label.offset = NULL, tipcex = 0.8, statecex = 0.7, 
splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE, plotlegend = FALSE, 
legend_ncol = NULL, legend_cex = 1, cornercoords_loc = "auto", 
tr = NULL, tipranges = NULL, if_ties = "takefirst", pie_tip_statecex = 0.7, 
juststats = FALSE, xlab = "Millions of years ago", root.edge = TRUE, 
colors_list_for_states = NULL, skiptree = FALSE, show.tip.label = TRUE, 
tipcol = "black", dej_params_row = NULL, plot_max_age = NULL, 
skiplabels = FALSE, plot_stratum_lines = TRUE, include_null_range = NULL, 
plot_null_range = FALSE, simplify_piecharts = FALSE, tipboxes_TF = TRUE, 
tiplabel_adj = c(0.5), no.margin = FALSE, xlims = NULL, ylims = NULL) 