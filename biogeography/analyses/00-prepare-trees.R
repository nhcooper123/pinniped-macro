# Prep trees and data

# Load libraries
library(ape)

# Read in full tree
tree <- read.tree("data/median_tree.tre")

# Replace 0 length branches with 0.001
for(i in 1:length(tree$edge.length)){
  if(tree$edge.length[i] == 0){
    tree$edge.length[i] <- 0.001
  }
}

# Change some taxon names to make the tree less complicated
tree$tip.label <- gsub("Pliophoca_etrusca_USNM_171221_et_181419_et_181504_et_187580_et_243697_et_250290_et_254327_et_34734", "Pliophoca_etrusca_USNM", tree$tip.label)
tree$tip.label <- gsub("Homiphoca_sp_SAM_PQL_30080", "Homiphoca_sp", tree$tip.label)
tree$tip.label <- gsub("Monachini_indet_NMVP160399", "Monachini_indet", tree$tip.label)
tree$tip.label <- gsub("cf_Miroungini_indet_CD35", "cf_Miroungini_indet", tree$tip.label)
tree$tip.label <- gsub("Phocidae_aff_Homiphoca_capensis_USNM", "aff_Homiphoca_capensis", tree$tip.label)
tree$tip.label <- gsub("Desmatophocidae_indet_USNM_335445", "Desmatophocidae_indet", tree$tip.label)
tree$tip.label <- gsub("Odobenidae_gen_et_spec_indet_LACM_135920", "Odobenidae_indet", tree$tip.label)
tree$tip.label <- gsub("Pelagiarctos_sp_SDNHM_131041", "Pelagiarctos_sp", tree$tip.label)

#-----------------------------------
# 9 areas tree
#-----------------------------------
# Read in 9 areas geographical data
geo1 <- read.delim("biogeography/data/pinniped-all-geography_9areas.txt")

# Drop species not in the geography file
tree1 <- drop.tip(tree, setdiff(tree$tip.label, geo1$X120))

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
tree3 <- drop.tip(tree, setdiff(tree$tip.label, geo3$X85))

# Save tree
write.tree(tree3, "biogeography/data/pinniped-tree-fossil_9areas.tre")

