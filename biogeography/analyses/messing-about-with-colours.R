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
                          impossible = c("#000000", # null
                                         "#D400D4", "#24408E", "#008026","#FFED00", "#FF8C00", "#E40303", "#613915", "#FFAFC8", "#74D7EE",# A-I
                                        rep("grey50", times = (length(default_colours$ranges_list))-10)))


# Add specific colours for tip combinations
default_colours <- mutate(default_colours,
                          impossible = case_when(ranges_list == "AD" ~ "#732950",
                                                 ranges_list == "AH" ~ "#0000ff",
                                                 ranges_list == "BE" ~ "#00FFAA",
                                                 ranges_list == "CF" ~ "#FFED80",
                                                 ranges_list == "DH" ~ "#FF5500",
                                                 ranges_list == "DI" ~ "#FF3990",
                                                 ranges_list == "ADH" ~ "#732982",
                                                 ranges_list == "BCG" ~ "#00AAFF",
                                                 ranges_list == "CDI" ~ "#AA718E",
                                                 ranges_list == "CEF" ~ "#558080",
                                                 ranges_list == "EFG" ~ "#8EAA39",
                                                 TRUE ~ as.character(impossible)))

colour_list_impossible <- pull(default_colours, impossible)

#-----------------------------------------
# Legend
#----------------------------------------
plot(NULL, xaxt = 'n', yaxt = 'n',bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("topleft", legend = c("A", "B", "C", "D", "E", "F", "G", "H", "I", 
                            "AD", "AH", "BE", "CF", "DH", "DI", 
                            "ADH", "BCG", "CDI", "CEF", "EFG"), 
       pch = 15, pt.cex = 2.5, cex = 1.2, bty = 'n', ncol = 4,
       col = c("#D400D4", "#24408E", "#008026","#FFED00", "#FF8C00", "#E40303", "#613915", "#FFAFC8", "#74D7EE",
               "#732950","#0000ff","#00FFAA","#FFED80","#FF5500", "#FF3990",
               "#732982", "#00AAFF", "#AA718E","#558080","#8EAA39" ))
mtext("Areas", at = 0.35, cex = 1.5)

### SAVE!