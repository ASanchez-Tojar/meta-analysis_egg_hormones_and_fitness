################################################################################
# Authors: 
#
# Lucia Mentesana (luciamentesana@gmail.com) wrote most of the original code
#
# Alfredo Sanchez-Tojar (alfredo.tojar@gmail.com) revised & adjusted the code, 
# wrote some original code (nostly for the analytical part) and prepared the 
# code for publication

# Script first created in Feb 2023

################################################################################
# Description of script and Instructions
################################################################################

# Script used to generate figure with map for the following meta-analysis:

# Do egg hormones have fitness consequences in wild birds? A systematic review 
# and meta-analysis

# by Lucia Mentesana, Michaela Hau, Pietro B. D'Amelio, Nicolas M. Adreani, and
# Alfredo Sanchez-Tojar. 

# preprint available at: https://doi.org/10.1101/2024.10.29.620852

# Note that for all figures – they were exported as .svg format using R's export 
# function. Then, we used the free vector graphics editor Inkscape to enhance 
# their appearance. Specifically, we adjusted the size and font of the legends 
# and modified the colors of the data. Finally, we exported each figure in .png 
# format, which is the version we submitted to the journal.

################################################################################
# Packages needed
################################################################################

# install.packages("pacman") #if not already installed
# pacman::p_load(metafor, plyr, stringr, rotl, ape, orchaRd, dplyr, 
#                devtools, tidyverse, patchwork, R.rsp, emmeans, ggplot2, maps, 
#                phytools, flextable, wesanderson)
pacman::p_load(ggplot2,plyr,stringr)

# cleaning up
rm(list = ls())


################################################################################
# Functions:
################################################################################

# none

################################################################################
# Importing the raw data 
################################################################################

db_for_map <- read.csv("data/processed_data/meta-analysis_egg_hormones_processed_data_full.csv")

################################################################################
# Map showing from where we obtained the data
################################################################################

# I need to split the area into 
db_for_map[c('Location','Country')] <- str_split_fixed(db_for_map$Study_area, ';', 2)

# I will subset the data base so that I get one row per location. Like this, I 
# get one point per location independently of the number of effect sizes we have
# per place.
db_for_map_unique_rows <- db_for_map[-c(which(duplicated(db_for_map$Location))),]

db_for_map_unique_rows <- ddply(db_for_map_unique_rows, c("Location"), transform, size = count(Location)) 

db_for_map_unique_rows$size1 <- db_for_map_unique_rows$size.freq + 1


# create data for world coordinates using  
# map_data() function 
world_coordinates <- map_data("world") 

# create world map using ggplot() function 
# geom_map() function takes world coordinates as input to plot world map 
fig_map <- ggplot() + 
  geom_map( 
    data = world_coordinates, map = world_coordinates, 
    aes(long, lat, map_id = region), 
    color = "white", fill= "lightblue") + 
  geom_point(data = db_for_map_unique_rows, aes(x = Long_study_area, 
                                                y = Lat_study_area,
                                                size = size.freq),
             color="black", alpha = 0.7) +
  theme_classic()
