# alorenzetti 20211116

# description ####
# this script will load
# the required libs

# loading libs ####
library(pacman)

packs = c("tidyverse",
          "DESeq2",
          "tximport",
          "rtracklayer",
          "ggtext",
          "ggthemes",
          "ggrepel",
          "openxlsx",
          "eulerr",
          "topGO",
          "GO.db",
          "Rgraphviz")

p_load(char = packs)

# setting ggplot theme
theme_set(theme_bw())

# getting colors from tab10 color scheme
tab10 = list()
tab10$blue = ggthemes::ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[1]
tab10$red = ggthemes::ggthemes_data$tableau$`color-palettes`$regular$`Tableau 10`$value[3]
