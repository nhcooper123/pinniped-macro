# Prep trees and data

# Load libraries
library(ape)

# Read in full tree
tree <- read.tree("biogeography/raw-data/pinnipedia_mcc_clean.tre")

#-----------------------------------
# 9 areas tree
#-----------------------------------
# Read in 9 areas geographical data
geo1 <- read.delim("biogeography/data/pinniped-all-geography_9areas.txt")

# Drop species not in the geography file
tree1 <- drop.tip(tree, setdiff(tree$tip.label, geo1$X105))

# Save tree
write.tree(tree1, "biogeography/data/pinniped-tree-all_9areas.tre")

#-----------------------------------
# 9 areas tree - extant only
#-----------------------------------
# Read in 9 areas geographical data
geo2 <- read.delim("biogeography/data/pinniped-extant-geography_9areas.txt")

# Drop species not in the geography file
tree2 <- drop.tip(tree, setdiff(tree$tip.label, geo2$X34))

# Save tree
write.tree(tree2, "biogeography/data/pinniped-tree-extant_9areas.tre")

#-----------------------------------
# 9 areas tree - fossil only
#-----------------------------------
# Read in 9 areas geographical data
geo3 <- read.delim("biogeography/data/pinniped-fossil-geography_9areas.txt")

# Drop species not in the geography file
tree3 <- drop.tip(tree, setdiff(tree$tip.label, geo3$X71))

# Save tree
write.tree(tree3, "biogeography/data/pinniped-tree-fossil_9areas.tre")
