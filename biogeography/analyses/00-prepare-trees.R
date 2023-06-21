# Prep trees and data

# Load libraries
library(ape)

# Read in full tree
tree <- read.tree("biogeography/data/pinnipedia_mcc_clean.tre")

#-----------------------------------
# 9 areas tree
#-----------------------------------
# Read in 9 areas geographical data
geo <- read.delim("biogeography/data/pinniped-all-geography_9areas.txt")

# Drop species not in the geography file
tree2 <- drop.tip(tree, setdiff(tree$tip.label, geo$X105))

# Save tree
write.tree(tree2, "biogeography/data/pinniped-tree-all_9areas.tre")
