# Open libraries:
library(Claddis)
library(metatree)

# Set working directory:
setwd("Pinniped_metatree/TNT_output/")

# Get list of mrp files:
mrp.list <- list.files()[grep("mrp.nex", list.files())]

# Get list of mrp files:
trees.list <- list.files()[grep("mpts_plus_strict.nex", list.files())]

# Make tree files:
for(i in 1:length(trees.list)) {
  
  # Read in TNT trees and split into mpts and strict consensus:
  mytrees <- Trees2MPTsAndStrict(trees.list[i])
  
  # If the tree limit of 100000 was hit (i.e., not all MPTs are guranteed to have been sampled) or no MRP could be created due to sheer number of trees:
  if(length(mytrees$mpts) == 100000 || sum(mrp.list == gsub("tntmpts_plus_strict.nex", "mrp.nex", trees.list[i])) == 0) {
    
    # Create MRP filename:
    mrp.filename <- gsub("tntmpts_plus_strict.nex", "mrp.nex", trees.list[i])
    
    # If MRP file was generated:
    if(length(which(mrp.list == mrp.filename)) > 0) {
      
      # Remove from MRP list:
      mrp.list <- mrp.list[-which(mrp.list == mrp.filename)]
      
      # Delete raw file as likely too big anyway:
      file.remove(mrp.filename)
      
    }
    
    # Create nexus file name:
    nexus.filename <- gsub("mrp\\.nex", ".nex", mrp.filename)
    
    # Read in original matrix:
    mymatrix <- read_nexus_matrix(paste("Pinniped_metatree/Nexus/",  nexus.filename, sep = ""))
    
    # Write out regular TNT file:
    write_tnt_matrix(mymatrix, gsub("\\.nex", ".tnt", nexus.filename))
    
    # Read in TNT lines:
    TNT.lines <- readLines(gsub("\\.nex", ".tnt", nexus.filename))
    
    # Get TNT data block:
    tnt.block <- TNT.lines[1:(grep("proc/;", TNT.lines) - 1)]
    
    # Create analysis block:
    anal.block <- paste(c("rseed*;\nhold 999;\nxmult=rss fuse 50 drift 50 ratchet 50;\nmult 50 =tbr drift;\ntsave scratch.tre;\nsave;\ntsave /;", rep("rseed*;\nhold 999;\nxmult=rss fuse 50 drift 50 ratchet 50;\nmult 50 =tbr drift;\ntsave scratch.tre +;\nsave;\ntsave /;", 19), "hold 5000;\nshortread scratch.tre;\nbbreak=tbr;"), collapse = "\n")
    
    # Cretae empty vector to store final block:
    full.block <- vector(mode = "character")
    
    # Fill out all blocks for analysis:
    for(j in 1:20) full.block <- c(full.block, paste(paste(tnt.block, collapse = "\n"), anal.block, "mrp;", paste("export ", gsub("\\.nex", "", nexus.filename), "mrp_", j, ".nex;", sep = ""), sep = "\n"))
    
    # Write out TNT file:
    write(paste(paste(full.block, collapse = "\n"), "\nproc/;\n", sep = ""), gsub("\\.nex", ".tnt", nexus.filename))
    
  }
  
  # Make file name:
  file.name <- gsub("tntmpts_plus_strict.nex", "", trees.list[i])
  
  # # Write out MPTs:
  write(mytrees$mpts, paste("Pinniped_metatree/MPT", "/", file.name, ".tre", sep = ""))
  # 
  # # Write out first MPT:
  write(mytrees$mpts[1], paste("Pinniped_metatree/firstMPT", "/", file.name, ".tre", sep = ""))
  # 
  # # Write out strict consensus:
  write(mytrees$strict, paste("Pinniped_metatree/SC", "/", file.name, ".tre", sep = ""))
  # 
  # # Delete trees file now no longer needed:
  file.remove(trees.list[i])
  # 
  # Spit out loop position:
  cat(i, " ")
  
}

# Make mrp files:
for(i in 1:length(mrp.list)) {
  
# Add assumptions block to MRP:
  x <- paste(c(readLines(mrp.list[i]), "BEGIN ASSUMPTIONS;", "OPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;", "END;"), collapse = "\n")
  
  # Write out MRP file with assumptions added (can then be read in with read_nexus_matrix):
  write(x = x, file = mrp.list[i])
  
  # Read in MRP file:
  mymrp <- read_nexus_matrix(mrp.list[i])
  
  # Remove root taxon:
  
      mymrp$matrix_1$matrix <- mymrp$matrix_1$matrix[-which(rownames(mymrp$matrix_1$matrix) == "ROOT"), ]
  
  # Collapse to just unique characters:
  mymrp <- compactify_matrix(mymrp)
  
  # Overwrite weights (set all to one):
  mymrp$matrix_1$character_weights <- rep(1, length(mymrp$matrix_1$character_weights))
  
  # Make file name:
  file.name <- gsub(".nex", "", mrp.list[i])
  
  # Isolate MPR taxon names:
  mrp.names <- rownames(mymrp$matrix_1$matrix)
  
  # Isolate full names:
  nexus.names <- rownames(read_nexus_matrix(paste("Pinniped_metatree/Nexus/", gsub("mrp", "", file.name), ".nex", sep = ""))$matrix_1$matrix)
  
  # Check to see if MRP names are contracted:
  if(length(setdiff(mrp.names, nexus.names)) > 0) {
    
    # List all contracted names:
    contracted.names <- setdiff(mrp.names, nexus.names)
    
    # For each contracted name:
    for(j in 1:length(contracted.names)) {
      
      # Get matching full name(s):
      full.name <- nexus.names[grep(contracted.names[j], nexus.names)]
      
      # Check that there are not multiple matches:
      if(length(full.name) > 1) stop("Multiple names match contracted form. Check manually.")
      
      # Overwrite contracted name with full name:
      rownames(mymrp$matrix_1$matrix)[which(rownames(mymrp$matrix_1$matrix) == contracted.names[j])] <- full.name
      
    }
    
  }
  
  # Write out MRP in #NEXUS format:
  write_nexus_matrix(mymrp, paste("Pinniped_metatree/MRP", file.name, ".nex", sep = ""))
  
  # Delete file once finished:
  file.remove(mrp.list[i])
  
  # Spit out loop position:
  cat(i, " ")
  
}
