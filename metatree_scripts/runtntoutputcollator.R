# Load libraries:
rm(list=ls())
library(ape)
library(metatree)
# Set working directoyr to where tnt output resides:

setwd("Pinniped_metatree/metatree_files/Files/tnt/")
setwd("Pinniped_metatree/metatree_files/morphology_only/Files/tnt/")

# Get a list of the tree files found there:
treefiles <- list.files()[grep("trees_tnt", list.files())]
# Create empty lists to store results for each file:
AllTrees <- AllTreeLengths <- list()

# For every log file for which trees exist:

for(i in 1:length(treefiles)) {
    run_number<- sort(as.numeric(unlist(lapply(strsplit(treefiles, "trees_tnt|\\."), '[[', 3))))[i]
    # Build path to and read in current log file:
    currentlogfile <- readLines(paste("tnt_log.", run_number, ".txt", sep = ""))
    
    # Begin to whittle down to just the tree lengths by finding block that contains them:
    treelengths <- currentlogfile[(grep("Tree lengths", currentlogfile) + 1):(grep("COMMAND: LOG", currentlogfile) - 1)]
    
    # Collapse to only-non-empty lines:
    treelengths <- treelengths[lapply(lapply(strsplit(treelengths, " "), nchar), sum) > 0]
    
    # Remove top line:
    treelengths <- treelengths[2:length(treelengths)]
    
    # Trim ends away:
    treelengths <- gdata::trim(treelengths)
    
    # Remove leading count to just get tree lengths in order:
    treelengths <- as.numeric(unlist(lapply(strsplit(treelengths, " "), '[', -1)))
    
    # Store tree lengths in list:
    AllTreeLengths[[(i + 1)]] <- treelengths
    
}

# For each tree file:
for(i in 1:length(treefiles)) {
    run_number <- sort(as.numeric(unlist(lapply(strsplit(treefiles, "trees_tnt|\\."), '[[', 3))))[i]
    # Build path to tree file and read in:
    currenttreefile <- readLines(paste("trees_tnt.", run_number, ".tnt", sep = ""))
    
    # Get just the tree lines:
    treelines <- currenttreefile[unlist(lapply(strsplit(currenttreefile, ""), '[[', 1)) == "("]
    
    # Update ending of Newick string to semicolon:
    treelines <- gsub("*", ";", treelines, fixed = TRUE)
    
    # Update spaces to commas to move towards Newick format:
    treelines <- gsub(" ", ",", treelines, fixed = TRUE)
    
    # Add commas between clades to get to Newick format:
    treelines <- gsub(")(", "),(", treelines, fixed = TRUE)
    
    # Remove commas before close parantheses to get to Newick format:
    treelines <- gsub(",)", ")", treelines, fixed = TRUE)
    
    # Actually read in formatted data as tree files using ape:
    Trees <- read.tree(text = treelines)
    
    # Store tree files:
    AllTrees[[(i)]] <- Trees
    
}

# Convert all tree lengths to vector:
AllTreeLengths <- unlist(AllTreeLengths)

# Create empty new all trees to combine as a multiphylo object:
NewAllTrees <- list()

# For each tree file:
for(i in 1:length(AllTrees)) {
    if(class(AllTrees[[i]])=="phylo") {
      NewAllTrees[[length(NewAllTrees)+1]] <- ladderize(AllTrees[[i]])
    } else{
    # Store each indiviudal tree after ladderising it:
      for(j in 1:length(AllTrees[[i]])) NewAllTrees[[(length(NewAllTrees) + 1)]] <- ladderize(AllTrees[[i]][[j]])
    }
}


# Update to class multiPhylo:
class(NewAllTrees) <- "multiPhylo"

# Collapse to just minimum length trees:
NewAllTrees <- NewAllTrees[which(AllTreeLengths == min(AllTreeLengths))]

# Collapse to just unique trees:
NewAllTrees <- read.tree(text = unique(write.tree(NewAllTrees)))

# Get strict component consensus:
SCC <- ladderize(consensus(NewAllTrees))

# Get 50% MR consensus:
MRC <- ladderize(consensus(NewAllTrees,p=0.5))

# Write out MPTs:
write.tree(NewAllTrees, "../../MPTS/STR_MPTs.tre")
# Write out strict component consensus:
write.tree(SCC, "../../Consensus_Trees/STR_SCC.tre")
# Write out 50% MR consensus:
write.tree(MRC, "../../Consensus_Trees/STR_MRC.tre")


# Safely reinsert exclude taxa and write out to file:
 # Claddis::safe_taxonomic_reinsertion(input_filename = "../../MPTS/STR_MPTs.tre", output_filename  = "../../MPTS/MPTs.tre", 
 #                                     str_taxa = read.table("../STR.txt", header = TRUE, 
 #                                     stringsAsFactors = FALSE), multiple_placement_option = "random")
## fixed function by gjs to deal with similar names
str_fix(input_filename = "../../MPTS/STR_MPTs.tre", output_filename  = "../../MPTS/MPTs.tre", 
        str_taxa = read.table("../STR.txt", header = TRUE, 
                              stringsAsFactors = FALSE), multiple_placement_option = "random")
# Get exclude strict consensus and write to file:
MPTs <- ape::read.tree("../../MPTS/MPTs.tre")
phytools:::writeNexus(MPTs, "../../MPTS/MPTs.nex")
SCC <- ape::consensus(MPTs)
MR50<- ape::consensus(MPTs, p=0.5)
write.tree(root(MR50, outgroup = "allzero", resolve.root = T), "../../Consensus_Trees/MRC.tre")
write.tree(root(SCC, outgroup = "allzero", resolve.root = T), "../../Consensus_Trees/SCC.tre")

