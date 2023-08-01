rm(list=ls())
library(ape);
library(geiger)
library(paleotree)
source("Pinniped_metatree/Scripts/BuildBeastXML.R")

setwd("Pinniped_metatree/TimeTree/")
### read in constraint tree (MRP consensus, MRL tree etc )
### then read in stratigraphic data in csv format
### strat data should be csv file with name, fad, lad
### finally, read in molecular alignment for extant taxa

constraintTree <- read.tree("MRC.tre")
StratRanges <- read.csv("Taxa_FAD.csv", row.names = 1, stringsAsFactors = F)
StratRanges <- StratRanges[,-3]
colnames(StratRanges) <- c("fad", "lad")
rownames(StratRanges)<-gsub(pattern = " ", replacement = "_",rownames(StratRanges))
alignment <- read.nexus.data("pinnipeds_extant_beast.nex")
extant <- matrix(data=rep(0, 2*length(alignment)), nrow = length(alignment), ncol=2,dimnames = list(names(alignment), c("fad", "lad")))
StratRanges<-rbind(StratRanges, extant)
write.csv(as.matrix(setdiff(constraintTree$tip.label, rownames(StratRanges))), "missing_ages.csv")
# if taxon names in strat range contain spaces, uncomment and run line below
#rownames(StratRanges)<-gsub(" ", "_", rownames(StratRanges))

#if age ranges larger than some amount are not desirable run two lines below
max.age.range <- 5
StratRanges<-StratRanges[-which(StratRanges[,1]-StratRanges[,2]>max.age.range),] ## remove taxa with large strat ranges


write.csv(as.matrix(geiger:::name.check(constraintTree, StratRanges)$tree_not_data), "treenotdata.csv")
td <- treedata(constraintTree, StratRanges)


plot(constraintTree)
constraintTree <- td$phy
StratRanges <- td$data

fileStem = "Pinnipeds_BEAST"

PrepareBeastMetatree(constraintTree = constraintTree, 
                     StratRanges = StratRanges, 
                     alignment = alignment, 
                     myfileName = fileStem,
                     makeStartTree = TRUE, 
                     start.tree.method="mbl", 
                     vartime=0.5)

## output will be a series of files that can be used to create a BEAST xml through beauti
## i.e. nexus file and mean taxon ages
## as well as files that can be pasted into the xml (monophyly constraints, starting tree, age ranges).

