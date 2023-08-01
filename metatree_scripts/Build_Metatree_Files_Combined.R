# SCRIPT TO BUILD METATREE FILES
rm(list=ls())
# Load metatree library:
library(metatree)
library(Claddis)

# Put your MRP files from TNT in a folder called MRP and your XML files in a folder called XML
# Set directory paths:
MRPDirectory <- "Pinniped_metatree/InputData/MRP/"
XMLDirectory <- "Pinniped_metatree/InputData/XML/"

# Standard exclusive data list (supertrees and the like):
# add names of any datasets you want to exclude here - e.g. "Koretsky_et_Rahmat_2015a"
setwd("Pinniped_metatree/InputData/XML/")
ExclusiveDataList <- c("Rule_etal_2020a", 
                       "Koretsky_2001a", 
                       "Koretsky_et_Grigorescu_2002a", 
                       "Koretsky_et_Holec_2002a",
                       "Koretsky_et_Rahmat_2013a",
                       "Koretsky_et_Rahmat_2015a",
                       "Koretsky_etal_2016a")

# Build safe pinniped metatree: for right now, this is just morphology!!
# the relative weights here correspond to 
			# 1) the input weights (the weights read in from the source MRP files), 
			# 2) the publication year weights (from equation 1 in the supplement of Lloyd et al. 2016), 
			# 3) the data set dependency weights (1 / the number of "sibling" data sets; see Lloyd et al. 2016), and 
			# 4) the within-matrix weights of individual clades (1 / number of conflicitng clades). 
#Zeroes exclude particular weighting types. These are the settings described in the whale paper

Pinnipedia_safe <- metatree::Metatree(MRPDirectory = MRPDirectory, XMLDirectory = XMLDirectory, 
	TargetClade = "Panpinnipedia", InclusiveDataList = c(), ExclusiveDataList = ExclusiveDataList, 
	MissingSpecies = "exclude", RelativeWeights = c(0, 1, 1, 1), 
	WeightCombination = "product", ReportContradictionsToScreen = FALSE)


## Write Morphology only output to file 

## Write out safe metatree files:
Claddis::write_nexus_matrix(Pinnipedia_safe$FullMRPMatrix, "Pinniped_metatree/metatree_files/morphology_only/Files/FULL.nex")
Claddis::write_nexus_matrix(Pinnipedia_safe$STRMRPMatrix, "Pinniped_metatree/metatree_files/morphology_only/Files/STR.nex")
Claddis::write_tnt_matrix(Pinnipedia_safe$FullMRPMatrix, "Pinniped_metatree/metatree_files/morphology_only/Files/FULL.tnt")
Claddis::write_tnt_matrix(Pinnipedia_safe$STRMRPMatrix, "Pinniped_metatree/metatree_files/morphology_only/Files/STR.tnt")
write.table(Pinnipedia_safe$SafelyRemovedTaxa, "Pinniped_metatree/metatree_files/morphology_only/Files/STR.txt", row.names = FALSE)
write.table(Pinnipedia_safe$CharacterWeights, "Pinniped_metatree/metatree_files/morphology_only/Files/CharacterWeights.txt", row.names = TRUE)
write.table(Pinnipedia_safe$DataSetWeights, "Pinniped_metatree/metatree_files/morphology_only/Files/DataSetWeights.txt", row.names = TRUE)
ape::write.tree(phy = Pinnipedia_safe$TaxonomyTree, file = "Pinniped_metatree/metatree_files/morphology_only/Files/TaxonomyTree.tre")


# Build Molecular only MRPs:
Pinnipedia_mol1<- metatree::Metatree(MRPDirectory = "Pinniped_metatree/InputData/Molecular_Data/Lopes/mrp", 
                    XMLDirectory = "Pinniped_metatree/InputData/Molecular_Data/Lopes/xml", 
                    TargetClade = "Panpinnipedia", InclusiveDataList = c(), ExclusiveDataList = c(), 
                    MissingSpecies = "exclude", RelativeWeights = c(1, 0, 0, 0), WeightCombination = "product", 
                    ReportContradictionsToScreen = FALSE, ExcludeTaxonomyMRP = TRUE)


# Prune any now duplicated characters in molecular data sets:
Pinnipedia_mol1$FullMRPMatrix <- metatree::PisaniMRPPrune(Pinnipedia_mol1$FullMRPMatrix)
# Place  character weights on 999 to 1000 scale:
Pinnipedia_mol1$FullMRPMatrix$matrix_1$character_weights <- plotrix::rescale(Pinnipedia_mol1$FullMRPMatrix$matrix_1$character_weights, c(999, 1000))

# Add in taxa missing from molecular as all NAs and add taxa missing from morphology to that matrix as NAs too:
missing_taxa <- setdiff(rownames(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix),rownames(Pinnipedia_mol1$FullMRPMatrix$matrix_1$matrix))
mol.dim<-dim(Pinnipedia_mol1$FullMRPMatrix$matrix_1$matrix) 
tmp <-rbind(Pinnipedia_mol1$FullMRPMatrix$matrix_1$matrix, matrix(data=rep(NA, length=(length(missing_taxa)*mol.dim[2])), nrow=length(missing_taxa), ncol=mol.dim[2]))
rownames(tmp) <- c(rownames(Pinnipedia_mol1$FullMRPMatrix$matrix_1$matrix), missing_taxa)

missing_taxa2<-setdiff(rownames(Pinnipedia_mol1$FullMRPMatrix$matrix_1$matrix), rownames(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix))
tmp2 <-rbind(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix, matrix(data=rep(NA, length=(length(missing_taxa2)*dim(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix)[2])), nrow=length(missing_taxa2), ncol=dim(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix)[2]))
rownames(tmp2) <- c(rownames(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix), missing_taxa2)

tmp2<-tmp2[match(rownames(tmp), rownames(tmp2)), ]

Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix <- tmp2
Pinnipedia_mol1$FullMRPMatrix$matrix_1$matrix <- tmp


# Duplicate molecular data twenty times to upweight it:
Pinnipedia_mol1$FullMRPMatrix <- metatree::EmbiggenMatrix(Pinnipedia_mol1$FullMRPMatrix, 100)



# Combine molecular and morphology MRPs:
Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix <- cbind(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix, Pinnipedia_mol1$FullMRPMatrix$matrix_1$matrix)
Pinnipedia_safe$FullMRPMatrix$matrix_1$ordering <- c(Pinnipedia_safe$FullMRPMatrix$matrix_1$ordering, Pinnipedia_mol1$FullMRPMatrix$matrix_1$ordering)
Pinnipedia_safe$FullMRPMatrix$matrix_1$character_weights <- c(Pinnipedia_safe$FullMRPMatrix$matrix_1$character_weights, Pinnipedia_mol1$FullMRPMatrix$matrix_1$character_weights)
Pinnipedia_safe$FullMRPMatrix$matrix_1$minimum_values <- c(Pinnipedia_safe$FullMRPMatrix$matrix_1$minimum_values, Pinnipedia_mol1$FullMRPMatrix$matrix_1$minimum_values)
Pinnipedia_safe$FullMRPMatrix$matrix_1$maximum_values <- c(Pinnipedia_safe$FullMRPMatrix$matrix_1$maximum_values, Pinnipedia_mol1$FullMRPMatrix$matrix_1$maximum_values)
Pinnipedia_safe$CharacterWeights <- rbind(Pinnipedia_safe$CharacterWeights, matrix(data = c((nrow(Pinnipedia_safe$CharacterWeights) + 1):ncol(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix), Pinnipedia_safe$FullMRPMatrix$matrix_1$character_weights[(nrow(Pinnipedia_safe$CharacterWeights) + 1):ncol(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix)]), ncol = 2, dimnames = list(rep("Pinnipedia_mol", times = ncol(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix) - nrow(Pinnipedia_safe$CharacterWeights)), c("Character", "Weight"))))

SafeSTR <- Claddis::safe_taxonomic_reduction(Pinnipedia_safe$FullMRPMatrix)

# Update STR matrices:
Pinnipedia_safe$STRMRPMatrix <- SafeSTR$reduced_matrix

# Update STR lists:
Pinnipedia_safe$SafelyRemovedTaxa <- SafeSTR$str_taxa


## now molecular 2

Pinnipedia_mol2<- metatree::Metatree(MRPDirectory = "Pinniped_metatree/InputData/Molecular_Data/fulton/mrp", 
                                     XMLDirectory = "Pinniped_metatree/InputData/Molecular_Data/fulton/xml", 
                                     TargetClade = "Panpinnipedia", InclusiveDataList = c(), ExclusiveDataList = c(), 
                                     MissingSpecies = "exclude", RelativeWeights = c(1, 0, 0, 0), WeightCombination = "product", 
                                     ReportContradictionsToScreen = FALSE, ExcludeTaxonomyMRP = TRUE)


# Prune any now duplicated characters in molecular data sets:
Pinnipedia_mol2$FullMRPMatrix <- metatree::PisaniMRPPrune(Pinnipedia_mol2$FullMRPMatrix)
# Place  character weights on 999 to 1000 scale:
Pinnipedia_mol2$FullMRPMatrix$matrix_1$character_weights <- plotrix::rescale(Pinnipedia_mol2$FullMRPMatrix$matrix_1$character_weights, c(999, 1000))

# Add in taxa missing from molecular as all NAs and add taxa missing from morphology to that matrix as NAs too:
missing_taxa <- setdiff(rownames(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix),rownames(Pinnipedia_mol2$FullMRPMatrix$matrix_1$matrix))
mol.dim<-dim(Pinnipedia_mol2$FullMRPMatrix$matrix_1$matrix) 
tmp <-rbind(Pinnipedia_mol2$FullMRPMatrix$matrix_1$matrix, matrix(data=rep(NA, length=(length(missing_taxa)*mol.dim[2])), nrow=length(missing_taxa), ncol=mol.dim[2]))
rownames(tmp) <- c(rownames(Pinnipedia_mol2$FullMRPMatrix$matrix_1$matrix), missing_taxa)

missing_taxa2<-setdiff(rownames(Pinnipedia_mol2$FullMRPMatrix$matrix_1$matrix), rownames(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix))
tmp2 <-rbind(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix, matrix(data=rep(NA, length=(length(missing_taxa2)*dim(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix)[2])), nrow=length(missing_taxa2), ncol=dim(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix)[2]))
rownames(tmp2) <- c(rownames(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix), missing_taxa2)

tmp2<-tmp2[match(rownames(tmp), rownames(tmp2)), ]

Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix <- tmp2
Pinnipedia_mol2$FullMRPMatrix$matrix_1$matrix <- tmp


# Duplicate molecular data twenty times to upweight it:
Pinnipedia_mol2$FullMRPMatrix <- metatree::EmbiggenMatrix(Pinnipedia_mol2$FullMRPMatrix, 100)



# Combine molecular and morphology MRPs:
Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix <- cbind(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix, Pinnipedia_mol2$FullMRPMatrix$matrix_1$matrix)
Pinnipedia_safe$FullMRPMatrix$matrix_1$ordering <- c(Pinnipedia_safe$FullMRPMatrix$matrix_1$ordering, Pinnipedia_mol2$FullMRPMatrix$matrix_1$ordering)
Pinnipedia_safe$FullMRPMatrix$matrix_1$character_weights <- c(Pinnipedia_safe$FullMRPMatrix$matrix_1$character_weights, Pinnipedia_mol2$FullMRPMatrix$matrix_1$character_weights)
Pinnipedia_safe$FullMRPMatrix$matrix_1$minimum_values <- c(Pinnipedia_safe$FullMRPMatrix$matrix_1$minimum_values, Pinnipedia_mol2$FullMRPMatrix$matrix_1$minimum_values)
Pinnipedia_safe$FullMRPMatrix$matrix_1$maximum_values <- c(Pinnipedia_safe$FullMRPMatrix$matrix_1$maximum_values, Pinnipedia_mol2$FullMRPMatrix$matrix_1$maximum_values)
Pinnipedia_safe$CharacterWeights <- rbind(Pinnipedia_safe$CharacterWeights, matrix(data = c((nrow(Pinnipedia_safe$CharacterWeights) + 1):ncol(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix), Pinnipedia_safe$FullMRPMatrix$matrix_1$character_weights[(nrow(Pinnipedia_safe$CharacterWeights) + 1):ncol(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix)]), ncol = 2, dimnames = list(rep("Pinnipedia_mol", times = ncol(Pinnipedia_safe$FullMRPMatrix$matrix_1$matrix) - nrow(Pinnipedia_safe$CharacterWeights)), c("Character", "Weight"))))

SafeSTR <- Claddis::safe_taxonomic_reduction(Pinnipedia_safe$FullMRPMatrix)

# Update STR matrices:
Pinnipedia_safe$STRMRPMatrix <- SafeSTR$reduced_matrix

# Update STR lists:
Pinnipedia_safe$SafelyRemovedTaxa <- SafeSTR$str_taxa


## Write out safe metatree files:
Claddis::write_nexus_matrix(Pinnipedia_safe$FullMRPMatrix, "Pinniped_metatree/metatree_files/Files/FULL.nex")
Claddis::write_nexus_matrix(Pinnipedia_safe$STRMRPMatrix, "Pinniped_metatree/metatree_files/Files/STR.nex")
Claddis::write_tnt_matrix(Pinnipedia_safe$FullMRPMatrix, "Pinniped_metatree/metatree_files/Files/FULL.tnt")
Claddis::write_tnt_matrix(Pinnipedia_safe$STRMRPMatrix, "Pinniped_metatree/metatree_files/Files/STR.tnt")
write.table(Pinnipedia_safe$SafelyRemovedTaxa, "Pinniped_metatree/metatree_files/Files/STR.txt", row.names = FALSE)
write.table(Pinnipedia_safe$CharacterWeights, "Pinniped_metatree/metatree_files/Files/CharacterWeights.txt", row.names = TRUE)
write.table(Pinnipedia_safe$DataSetWeights, "Pinniped_metatree/metatree_files/Files/DataSetWeights.txt", row.names = TRUE)
ape::write.tree(phy = Pinnipedia_safe$TaxonomyTree, file = "Pinniped_metatree/metatree_files/Files/TaxonomyTree.tre")

