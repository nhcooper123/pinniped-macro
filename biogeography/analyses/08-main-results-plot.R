# Pretty BGB results plots
# Make trees using ggtree.
#------------------------
# Load libraries
#------------------------
library(BioGeoBEARS)
library(geiger)
library(ape)
library(tidyverse)
library(ggtree)
library(ggimage)
library(patchwork)

#-------------------
# Read in the tree
#------------------
tree <- read.tree("biogeography/data/pinniped-tree-all_9areas.tre")
# Replace underscores with spaces in tree file
tree$tip.label <- gsub("_", " ", tree$tip.label)
# Reorder to match BGB code
tree <- reorder(tree, "pruningwise")


#-----------------
# Read in data
#-----------------
areas <- read_PHYLIP_data("biogeography/data/pinniped-all-geography_9areas.txt")
# Replace underscores with spaces
rownames(areas) <- gsub("_", " ", rownames(areas))

# To get different colours for each area 
# we need to replace 1 with numbers 2-9 IN B-I
areas2 <-  mutate(areas, B = case_when(B == "1" ~ "2", TRUE ~ as.character(B)),
                       C = case_when(C == "1" ~ "3", TRUE ~ as.character(C)),
                       D = case_when(D == "1" ~ "4", TRUE ~ as.character(D)),
                       E = case_when(E == "1" ~ "5", TRUE ~ as.character(E)),
                       F = case_when(F == "1" ~ "6", TRUE ~ as.character(F)),
                       G = case_when(G == "1" ~ "7", TRUE ~ as.character(G)),
                       H = case_when(H == "1" ~ "8", TRUE ~ as.character(H)),
                       I = case_when(I == "1" ~ "9", TRUE ~ as.character(I)))
#-----------------
# Define colours
#-----------------
basic_cols <- c("#eeeeee", "#D400D4", "#24408E", "#008026","#FFED00",
                "#FF8C00", "#E40303", "#613915", "#FFAFC8", "#74D7EE")

#----------------------------
# Identify MRCA of each group
# -----------------------------
# Read in file with taxonomy info
groups <- read.csv("data/taxa_groups.csv")
# Replace underscores with spaces
groups$Taxon <- gsub("_", " ", groups$Taxon)

# Match to tree
check <- name.check(tree, groups, data.names = groups$Taxon)

# Remove species not in the tree
matches <- match(groups$Taxon, check$data_not_tree, nomatch = 0)
groups <- subset(groups, matches == 0)

# Get node numbers
# Note that ggtree does something odd to node numbering meaning these don't all
# match 100%. Check using
# ggtree(tree, branch.length = "none") + geom_text2(aes(subset=!isTip, label=node), size =2,  hjust=-.3) + geom_tiplab(size = 2)
phocid <- 159
otarid <- 187
walrus <- 213
desmo <- 175
mona <- 134
dev <- 174

# Make df of numbers and names
df <- data.frame(node = c(phocid, otarid, walrus, desmo, mona, dev),
                 name = c("Phocinae", "Otariidae", "Odobenidae", "Desmatophocidae", "Monachinae", "Devinophocinae"))

#-----------------------------
# Plot tree with areas at tips
#-----------------------------
# Make the tree base

# Rotate branches to match order of bamm results:
#stem, odoben, otarid, desm, devi, phoci, monc
base <- 
  ggtree(tree) %>% rotate(121) +
  xlim(0,100) +
  geom_tiplab(geom = "text", fontface = "italic", size = 2)

#base
#------------------
# Add areas
# Ignore warning message about scales
area_plot <-
  gheatmap(base, areas2, offset = 15, width = 0.2,
         font.size = 2, colnames_position = "top", color = "black")+
         scale_fill_manual(values = basic_cols) +
  theme(legend.position = "none")
  
#area_plot  
#------------------
# Add clade labels
area_group_plot <- 
  area_plot + geom_cladelab(data = df, mapping = aes(node = node, label = name),
              offset = c(22,22,35.2,32.1,22,22), offset.text = 1)

#-------------------
# Add ML states
# Need to get results from BGB
source("biogeography/analyses/07-Extract-BGB-results-for-plotting.R")

#-------------------------------------------------
# Colours
#-------------------------------------------------
# Run this to see what states require colours
base + geom_nodepoint(aes(colour = MLstates), size = 2) 
# A-I (no F), ABCG, AD, DI

dff <- data.frame(state = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "ABCG", "AD", "DH", "DI"), 
                  col = c("#D400D4", "#24408E", "#008026","#FFED00", "#FF8C00", "#E40303", "#613915", "#FFAFC8", "#74D7EE",
                          "#732982", "#8EAA39" , "#0000ff", "#000000"))
        
# Create legend
# png(file = "biogeography/outputs/main-legend.png", width = 4000, height = 3100, res = 900)
plot(NULL, xaxt = 'n', yaxt = 'n',bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("topleft", legend = dff$state, 
       pch = 15, pt.cex = 2.4, cex = 1.1, bty = 'n', ncol = 3,
       col = dff$col)
#dev.off()

#------------------------------------------------------------
# Plot all areas
#------------------------------------------------------------
# Create list of colours, remove non needed states first
# Also need to arrange in alphabetical order
new_colours <- dff %>%
  arrange(state) %>%
  filter(state != "F" & state != "DI") %>%
  pull(col)

# Plot 
area_group_ML_plot <- area_group_plot + geom_nodepoint(aes(colour = MLstates), size = 2) +
  scale_colour_manual(values = new_colours) +
  theme(legend.position = "none")

#------------------
# Save plot
#------------------
ggsave(area_group_ML_plot, file = "biogeography/outputs/biogeography-nice-figure.png", 
       width = 9, height = 7, dpi = 900)

#-----------------------
# Inset
#-----------------------
# Choose one taxon per group
# Needs to be living if possible
focal_species <- c("Arctocephalus australis", "Phoca vitulina", 
                   "Odobenus rosmarus", "Allodesmus demerei", 
                   "Potamotherium vallentoni", "Devinophoca claytoni", 
                   "Monachus monachus")

# Replace underscores with spaces in tree file
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Remove all other taxa from the tree
tree_basic <- drop.tip(tree, setdiff(tree$tip.label, focal_species))

# Rename with group names
tree_basic$tip.label <- gsub("Arctocephalus australis", "Otariidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Phoca vitulina", "Phocinae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Odobenus rosmarus", "Odobenidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Allodesmus demerei", "Desmatophocidae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Potamotherium vallentoni", "stem", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Devinophoca claytoni", "Devinophocinae", tree_basic$tip.label) 
tree_basic$tip.label <- gsub("Monachus monachus", "Monachinae", tree_basic$tip.label) 

basic_tree <-
  ggtree(tree_basic, branch.length = "none") %>% rotate(8) %>% rotate(13) %>% rotate(12) +
  xlim(-2,10) +
  geom_tiplab(geom = "text", size = 6) +
  geom_rootedge(rootedge = 1) 
  #+ geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab(size = 2)

#------------------------------
# Add pies
# Get node numbers. Check using
# ggtree(tree, branch.length = "none") + geom_text2(aes(subset=!isTip, label=node), size =2,  hjust=-.3) + geom_tiplab(size = 2)

dd3 <- filter(dd2, node == 121 | node == 130|  node == 131
                  | node == 132 | node == 133 | node == 186)
# Change numbers to match tree node numbers
dd3$node <- 8:13 # 7 taxa, so nodes are 8:13

# Identify states with probabilities > 0.1
colnames(dd3)[which(dd3[1,] > 0.1)]
colnames(dd3)[which(dd3[2,] > 0.1)]
colnames(dd3)[which(dd3[3,] > 0.1)]
colnames(dd3)[which(dd3[4,] > 0.1)]
colnames(dd3)[which(dd3[5,] > 0.1)]
colnames(dd3)[which(dd3[6,] > 0.1)]


# Create colour palette just for these colours
colours_pies <- data.frame(state = colnames(dd2), col = rep("#eeeeee", length(colnames(dd2))))

colours_pies2 <- 
  colours_pies %>%
  dplyr::slice(-c(257,258,259)) %>%
  mutate(col = case_when(state == "A" ~ "#D400D4",
                         state == "AD" ~ "#8EAA39",
                         state == "D" ~ "#FFED00",
                         state == "DH" ~ "#0000ff",
                         state == "DI" ~ "#000000",
                         state == "I" ~ "#74D7EE",
                         TRUE ~ as.character(col)))

# Extract colours only
colours_pies3 <- 
  colours_pies2 %>%
  pull(col)

# Create pies
# Colour argument gives strange warning but this is just related to change in base R
# Note weird colour set up is because the state - gets removed for some reason
# but needs to be the first in the list to match dd3
pies <- nodepie(dd3, cols = 1:256, alpha = 1, color = c("#eeeeee", colours_pies3))

# Add pies to plot
inset(basic_tree, pies, width = 0.4, height = 1)

# To save
inset_pies <- inset(basic_tree, pies, width = 0.3, height = 0.8)
#------------------
# Save plot
#------------------
ggsave(inset_pies, file = "biogeography/outputs/biogeography-inset.png", 
       width = 8, height = 6, dpi = 900)
