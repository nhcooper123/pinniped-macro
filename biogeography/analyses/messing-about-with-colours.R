#------------------------
# Colour matching
#------------------------
# A bit of a hack job but it gets there...

#-----------------------
library(BioGeoBEARS)
library(tidyverse)

#--------------------------
# Load results files
#--------------------------
load("biogeography/outputs/pinnipeds-all-DEC_9areas_impossible.Rdata")
#load("biogeography/outputs/pinnipeds-all-DEC_9areas_unlikely.Rdata")
#------------------------
# Load geographical data
#------------------------
geogfn <- "biogeography/data/pinniped-all-geography_9areas.txt"
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

#---------------------------------------------------
# Get colours for these ranges using biogeobears 
# defaults for impossible analyses
#---------------------------------------------------
colors_matrix <- get_colors_for_numareas(length(areas))
colors_list_for_states = mix_colors_for_states(colors_matrix, 
                                               states_list_0based, 
                                               plot_null_range = resDEC$inputs$include_null_range)
colors_list_for_states

# Combine
default_colours <- data.frame(ranges_list, colors_list_for_states)

#----------------------
# Make colour lists
#----------------------
# Add basic colours
default_colours <- mutate(default_colours, 
                          impossible = c("#eeeeee", # null
                                         "#D400D4", "#24408E", "#008026","#FFED00", "#FF8C00", "#E40303", "#613915", "#FFAFC8", "#74D7EE",# A-I
                                        rep("grey50", times = (length(default_colours$ranges_list))-10)))


# Add specific colours for tip combinations
default_colours <- mutate(default_colours,
                          impossible = case_when(#ranges_list == "AB" ~ "#732950",
                                                 #ranges_list == "AC" ~ "#4B0082",
                                                 ranges_list == "AD" ~  "#8EAA39",
                                                 #ranges_list == "AH" ~ "#00FFAA",
                                                 #ranges_list == "BC" ~ "#FFED80",
                                                 #ranges_list == "BE" ~ "#FF5500",
                                                 #ranges_list == "CF" ~ "#FF3990",
                                                 #ranges_list == "DG" ~ "#FFF333",
                                                 ranges_list == "DH" ~ "#0000ff",
                                                 ranges_list == "DI" ~ "#000000",
                                                 #ranges_list == "EI" ~ "#AA718E",
                                                 #ranges_list == "ABC" ~ "#558080",
                                                 #ranges_list == "ADH" ~ "#8EAA39",
                                                 #ranges_list == "BCF" ~ "#6A40AA",
                                                 #ranges_list == "BCG" ~ "#2A6AEA",
                                                 #ranges_list == "BDE" ~ "#235347",
                                                 #ranges_list == "CDI" ~ "#6A6AAA",
                                                 #ranges_list == "CEF" ~ "#FFCC00",
                                                 #ranges_list == "EFG" ~ "#FFA500",
                                                 ranges_list == "ABCG" ~  "#732982",
                                                 TRUE ~ as.character(impossible)))

colour_list_impossible <- pull(default_colours, impossible)


##### needs fixing
#-----------------------------------------
# Legend
#----------------------------------------
#png(file = "supplemental/figures/BGB-legend.png", width = 3500, height = 3100, res = 900)
#par(mar = c(1,1,1,1))
#plot(NULL, xaxt = 'n', yaxt = 'n',bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
#legend("topleft", legend = c("A", "B", "C", "D", "E", "F", "G", "H", "I", 
#                           "AB", "AC", "AD", "AH", "BC", "BE", "CF", "DG", "DH", "DI", "EI",
#                           "ABC", "ADH", "BCF", "BCG", "BDE", "CDI", "CEF", "EFG", "ABCG"), 
#      pch = 15, pt.cex = 2.4, cex = 1.1, bty = 'n', ncol = 4,
#      col = c("#D400D4", "#24408E", "#008026","#FFED00", "#FF8C00", "#E40303", "#613915", "#FFAFC8", "#74D7EE",
#               "#732950","#4B0082", "#0000ff", "#00FFAA", "#FFED80", "#FF5500", "#FF3990", "#FFF333","#732982","#00AAFF","#AA718E",
#               "#558080","#8EAA39", "#6A40AA","#2A6AEA", "#235347", "#6A6AAA", "#FFCC00","#FFA500","#800000" ))
#dev.off()
### SAVE!
