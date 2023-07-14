# Make world map for legend
library(tidyverse)
world <- map_data("world")
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    colour = "lightgray", fill = "lightgrey") +
  theme_void() +          
  theme(rect = element_rect(fill = "transparent"))  # Make background transparent

ggsave(file = "biogeography/outputs/worldmap-for-legend.png", dpi = 900, 
       bg = "transparent")
