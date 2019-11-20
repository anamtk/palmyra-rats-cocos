library(tidyverse)
setwd("/Users/Ana/Dropbox/2017\ Palmyra\ Seed\ Predation\ Paper/Analyses/Draft_1_Mar_2019/forest_cover_percents")
canopy <- read.csv("forest_cover.csv")

#this is the total canopy area of Palmyra
canopy %>%
  summarize(sum = sum(Area, na.rm=T))

#this is the canopy area by canopy type, categorized as C. nucifera, C.nucifera+other species, no C. nucifera
canopy %>%
  group_by(Forest_Typ) %>%
  summarize(sum = sum(Area/1661031, na.rm=T))
