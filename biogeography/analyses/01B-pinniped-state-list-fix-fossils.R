## Manual modification of states_list 
## (list of geographic ranges, a.k.a. the list of states in the state space)
## Code by Nick Matzke see http://phylo.wikidot.com/example-biogeobears-scripts#manual_modify_states
## This is to remove geographic ranges that are disjunct and therefore nonsensical states.

#-------------------
# Load libraries
#-------------------
library(BioGeoBEARS)
library(cladoRcpp)

#------------------------
# Load geographical data
#------------------------
geogfn <- "biogeography/data/pinniped-fossil-geography_9areas.txt"
tip.ranges <- getranges_from_LagrangePHYLIP(lgdata_fn = geogfn)

#----------------------------------------------
# Identify states
#----------------------------------------------
# Max range size
max_range_size <- 4
# Identify areas
areas <- getareas_from_tipranges_object(tip.ranges)

# Get the list of states/ranges, where each state/range is a list of areas, counting from 0
states_list_0based <- rcpp_areas_list_to_states_list(areas = areas, maxareas = max_range_size, 
                                                    include_null_range = TRUE)

# How many states/ranges? # n = 256
length(states_list_0based)

# Make the list of ranges
ranges_list <- NULL

for (i in 1:length(states_list_0based)){    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) ) {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

# Look at the ranges list
ranges_list

# How many states/ranges; n = 256
length(ranges_list)

#----------------------------------------------
# Remove disjunct ranges for each time bin
#----------------------------------------------
# TIME BIN 1 (youngest: 3.7 Ma)
#----------------------------------------------
# Identify non-adjacent ranges (impossible only)
nonadjacent1A <- c("AD","FI")
# Extract from range list
keep1A <- ranges_list %in% nonadjacent1A == FALSE
ranges_list_1A <- ranges_list[keep1A]
# Make the stratum-specific states list
states_list_0based_1Af<- states_list_0based[keep1A]

# Identify non-adjacent ranges (impossible and unlikely)
nonadjacent1B <- c("AD","FI", "AE", "AF", "AG", "AI", "BD", "BE", "BF", "BH", 
                   "BI", "CD", "CE", "CH", "CI", "DF", "DG", "EI", "FH", "GH", "GI")
# Extract from range list
keep1B <- ranges_list %in% nonadjacent1B == FALSE
ranges_list_1B <- ranges_list[keep1B]
# Make the stratum-specific states list
states_list_0based_1Bf <- states_list_0based[keep1B]

#----------------------------------------------
# TIME BIN 2 (5.32 Ma)
#----------------------------------------------
# Identify non-adjacent ranges (impossible only)
nonadjacent2A <- c("FI")
# Extract from range list
keep2A <- ranges_list %in% nonadjacent2A == FALSE
ranges_list_2A <- ranges_list[keep2A]
# Make the stratum-specific states list
states_list_0based_2Af <- states_list_0based[keep2A]

# Identify non-adjacent ranges (impossible and unlikely)
nonadjacent2B <- c("FI", "AE", "AF", "AG", "AI", "BD", "BE", "BF", "BH", 
                   "BI", "CD", "CE", "CH", "CI", "DF", "DG", "EI", "FH", "GH", "GI")
# Extract from range list
keep2B <- ranges_list %in% nonadjacent2B == FALSE
ranges_list_2B <- ranges_list[keep2B]
# Make the stratum-specific states list
states_list_0based_2Bf <- states_list_0based[keep2B]

#----------------------------------------------
# TIME BIN 3 (5.33 Ma)
#----------------------------------------------
# Identify non-adjacent ranges (impossible only)
nonadjacent3A <- c("AH", "FI")
# Extract from range list
keep3A <- ranges_list %in% nonadjacent3A == FALSE
ranges_list_3A <- ranges_list[keep3A]
# Make the stratum-specific states list
states_list_0based_3Af <- states_list_0based[keep3A]

# Identify non-adjacent ranges (impossible and unlikely)
nonadjacent3B <- c("AH", "FI", "AE", "AF", "AG", "AI", "BD", "BE", "BF", "BH", 
                   "BI", "CD", "CE", "CH", "CI", "DF", "DG", "EI", "FH", "GH", "GI")
# Extract from range list
keep3B <- ranges_list %in% nonadjacent3B == FALSE
ranges_list_3B <- ranges_list[keep3B]
# Make the stratum-specific states list
states_list_0based_3Bf <- states_list_0based[keep3B]

#----------------------------------------------
# TIME BIN 4 (5.97 Ma)
#----------------------------------------------
# Identify non-adjacent ranges (impossible only)
# Note that for models to fit for fossils Paratethys must be a possibility or
# the model breaks as the state is needed as an intermediate.
nonadjacent4A <- c("AH", "AI", "BI", "CI", "DI", "EI", "FI", "GI", "HI")
# Extract from range list
keep4A <- ranges_list %in% nonadjacent4A == FALSE
ranges_list_4A <- ranges_list[keep4A]
# Make the stratum-specific states list
states_list_0based_4Af <- states_list_0based[keep4A]

# Identify non-adjacent ranges (impossible and unlikely)
nonadjacent4B <- c("AH", "AI", "BI", "CI", "DI", "EI", "FI", "GI", "HI", "AE", "AF", "AG", 
                   "BD", "BE", "BF", "BH", "CD", "CE", "CH", "DF", "DG", "FH", "GH")
# Extract from range list
keep4B <- ranges_list %in% nonadjacent4B == FALSE
ranges_list_4B <- ranges_list[keep4B]
# Make the stratum-specific states list
states_list_0based_4Bf <- states_list_0based[keep4B]

#----------------------------------------------
# TIME BIN 5 (13.82 Ma)
#----------------------------------------------
# Same as time bin 3
states_list_0based_5Af <- states_list_0based_3Af
states_list_0based_5Bf <- states_list_0based_3Bf

#----------------------------------------------
# TIME BIN 6 (15.97 Ma)
#----------------------------------------------
# Identify non-adjacent ranges (impossible only)
nonadjacent6A <- c("AH")
# Extract from range list
keep6A <- ranges_list %in% nonadjacent6A == FALSE
ranges_list_6A <- ranges_list[keep6A]
# Make the stratum-specific states list
states_list_0based_6Af <- states_list_0based[keep6A]

# Identify non-adjacent ranges (impossible and unlikely)
nonadjacent6B <- c("AH", "AE", "AF", "AG", "AI", "BD", "BE", "BF", "BH", 
                   "BI", "CD", "CE", "CH", "CI", "DF", "DG", "EI", "FH", "GH", "GI")
# Extract from range list
keep6B <- ranges_list %in% nonadjacent6B == FALSE
ranges_list_6B <- ranges_list[keep6B]
# Make the stratum-specific states list
states_list_0based_6Bf <- states_list_0based[keep6B]

#----------------------------------------------
# TIME BIN 7 (20.44 Ma)
#----------------------------------------------
# Same as time bin 3
states_list_0based_7Af <- states_list_0based_3Af
states_list_0based_7Bf <- states_list_0based_3Bf

#----------------------------------------------
# TIME BIN 8 (37.8 Ma)
#----------------------------------------------
# Same as time bin 6
states_list_0based_8Af <- states_list_0based_6Af
states_list_0based_8Bf <- states_list_0based_6Bf