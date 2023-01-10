###----------------------------------------------------###
###   Advanced Nonparametric Statistics and Smoothing  ###
###     Density estimation over complicated domains    ###
###                      - Stops -                     ###
###----------------------------------------------------###


rm(list=ls())
graphics.off()

# ||||||||||||||
# Libraries ----
# ||||||||||||||

library(rjson)
library(geojsonio)
library(tidyverse)


# ||||||||||||||||
# Import data ----
# ||||||||||||||||

data <- fromJSON(file="data/output.json")


# ||||||||||
# Stops ----
# ||||||||||

haltes <- data[['haltes']]

coord = data.frame(matrix(ncol = 3, nrow = length(haltes)))
colnames(coord) <- c("lat", "long", "city_name")

for (i in 1:length(haltes)){
  coord[i,c(1:2)] = t(haltes[[i]][["geoCoordinaat"]])
  coord[i,3] = haltes[[i]][["omschrijvingGemeente"]]
}


good_cities <- c("Heverlee", "Kessel-Lo", "Leuven", "Wijgmaal", "Wilsele")

stops <- coord %>% filter(long >= 4.636 & long <= 4.77 & 
                                 lat >= 50.824 & lat <=50.944 &
                                 city_name %in% good_cities )

save(stops, file="data/stops.RData" )

