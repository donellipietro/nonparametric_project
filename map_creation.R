###----------------------------------------------------###
###   Advanced Nonparametric Statistics and Smoothing  ###
###     Density estimation over complicated domains    ###
###               - Density esyimation -               ###
###----------------------------------------------------###

rm(list=ls())
graphics.off()

# ||||||||||||||
# Libraries ----
# ||||||||||||||

# install.packages('rgdal')
# install.packages('mapproj')
# install.packages('rgeos')
# install.packages('sf')
# install.packages('concaveman')
# install.packages('geosphere')
# install.packages('rmapshaper')
# install.packages('FRK')

library(tidyverse)
library(rgdal)
library(rgeos)
library(ggplot2)
library(broom)
library(knitr)
library(mapproj)
library(maptools)

library(sf)
library(ggforce)
library(Rcpp)
library(tidyverse)

library(geosphere)
library(concaveman)

library(rmapshaper)
library(FRK)


# ||||||||||||||
# Functions ----
# ||||||||||||||

fill_gaps <- function(points, max_dist) {
  
  points <- data.frame(points)
  points <- points[order(points$order),]
  
  columns <- colnames(points)
  filled_points <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(filled_points) <- columns
  
  for(i in 1:(dim(points)[1]-1)){
    
    point1 <- points[i, c(1,2)]
    point2 <- points[i+1, c(1,2)]
    dist <- distm(point1, point2, fun = distHaversine)
    if(dist > max_dist && points[i,]$group == points[i+1,]$group){
      t = seq(0, 1, length = round(dist/max_dist)+1)
      newlong <- (1-t)*as.numeric(point1[1]) + t*as.numeric(point2[1])
      newlat <- (1-t)*as.numeric(point1[2]) + t*as.numeric(point2[2])
      newpoints <- data.frame(long = newlong, lat = newlat, 
                              order = points[i,]$order + t,
                              group = points[i,]$group)
      filled_points <- rbind(filled_points, newpoints)
    } else{
      filled_points <- rbind(filled_points, points[i, c("long", "lat", "order", "group")])
    }
  }
  
  return(filled_points)
  
}


# ||||||||||||||||
# Import data ----
# ||||||||||||||||

# Map
shapes_belgium_shp <- readOGR("data/sh_statbel_statistical_sectors_20200101.shp",
                              layer = "sh_statbel_statistical_sectors_20200101")

# Data aggregated by areas
data = read_csv("data/Data1_demographic_sector.csv", local = locale(encoding = "latin1"), na = c("", "NA", ".", "-"))

# ||||||||||||||||||
# Map of Leuven ----
# ||||||||||||||||||

theme_set(theme_bw())

# Belgium
proj4string(shapes_belgium_shp) = CRS("+init=EPSG:31370")
shapes_belgium_shp_trans = spTransform(shapes_belgium_shp, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
# plot(shapes_belgium_shp)
# shapes_belgium_shp_trans@data %>% filter(T_ARRD_NL == "Arrondissement Leuven")

# Leuven
shapes_leuven = shapes_belgium_shp_trans[shapes_belgium_shp_trans@data$T_ARRD_NL == "Arrondissement Leuven", ]
# plot(shapes_leuven)

# Conversion of shapes_leuven into a dataframe
shapes_leuven_fortify = fortify(shapes_leuven, region = "CS01012020") %>% as_tibble()

# Leuven boxes cordinates 
inner_leuven_coords = tibble(long = c(12.53, 12.58), lat = c(55.09, 55.125))
inner_leuven_coords_2 = tibble(long = c(4.63, 4.77), lat = c(50.82, 50.945))

# Indentification of inner Leuven areas
ids_inner_leuven = shapes_leuven_fortify %>%
  filter(between(long, inner_leuven_coords$long[1], inner_leuven_coords$long[2]),
         between(lat, inner_leuven_coords$lat[1], inner_leuven_coords$lat[2])) %>%
  pull(id) %>%
  unique()
names_inner_leuven = shapes_leuven@data %>%
  filter(CS01012020 %in% ids_inner_leuven) %>%
  pull(T_SEC_NL)

# Map of Leuven
shapes_leuven_fortify %>%
  left_join(
    data %>%
      select(id = ID) %>%
      mutate(info = "yes")) %>%
  mutate(info = ifelse(is.na(info), "no", "yes")) %>%
  ggplot() +
  geom_polygon(aes(long, lat, group = group, fill = info), color = "black") +
  coord_map()

# Map of inner Leuven
shapes_leuven_fortify %>%
  left_join(
    data %>%
      select(id = ID) %>%
      mutate(info = "yes")) %>%
  mutate(info = ifelse(is.na(info), "no", "yes")) %>%
  ggplot() +
  geom_polygon(aes(long, lat, group = group, fill=info), color = "black") +
  coord_map(xlim = c(inner_leuven_coords_2$long[1], inner_leuven_coords_2$long[2]),
            ylim = c(inner_leuven_coords_2$lat[1], inner_leuven_coords_2$lat[2]))


# |||||||||||||||||||||
# Map manipulation ----
# |||||||||||||||||||||

# Inner Leuven without external areas
only_inner_leuven <- shapes_leuven_fortify %>% 
  left_join(
    data %>% 
      select(id = ID) %>% 
      mutate(info = "yes")) %>% 
  drop_na()

# Plot of Inner Leuven without external areas
only_inner_leuven %>%
  ggplot() +
  geom_polygon(aes(long, lat, group = group, fill=info), color = "black", fill = "white") +
  coord_map(xlim = c(inner_leuven_coords_2$long[1], inner_leuven_coords_2$long[2]), 
            ylim = c(inner_leuven_coords_2$lat[1], inner_leuven_coords_2$lat[2]))


# Filling the gaps
filled_points <- fill_gaps(only_inner_leuven, max_dist = 10)

# Filled plot
ggplot(filled_points, aes(long, lat)) +
  geom_point(color = "black") +
  coord_map(xlim = c(inner_leuven_coords_2$long[1], inner_leuven_coords_2$long[2]), 
            ylim = c(inner_leuven_coords_2$lat[1], inner_leuven_coords_2$lat[2]))

# Concave hull
boundaries <- as.matrix(concaveman(as.matrix(filled_points[,c(2,1,3)]), concavity = 1, length_threshold = 0))
boundaries <- data.frame(boundaries)
colnames(boundaries) <- c("lat", "long", "order")

# Concave hull plot
ggplot(boundaries, aes(long, lat)) +
  geom_point(color = "black") +
  coord_map(xlim = c(inner_leuven_coords_2$long[1], inner_leuven_coords_2$long[2]), 
            ylim = c(inner_leuven_coords_2$lat[1], inner_leuven_coords_2$lat[2]))

# Removing extra points and duplicates
boundaries <- boundaries[boundaries$order %% 1 == 0,]
boundaries <- boundaries %>% distinct()

# Concave hull only original points plot
ggplot(boundaries, aes(long, lat)) +
  geom_point(color = "black") +
  coord_map(xlim = c(inner_leuven_coords_2$long[1], inner_leuven_coords_2$long[2]), 
            ylim = c(inner_leuven_coords_2$lat[1], inner_leuven_coords_2$lat[2]))

# Boundary Simplification
boundaries$id = 1
boundaries_polygon <- df_to_SpatialPolygons(boundaries,"id",c("long","lat"),CRS())
boundaries_gSimplify <- ms_simplify(boundaries_polygon, keep = 0.05)

boundaries_simplified = data.frame(boundaries_gSimplify@polygons[1][[1]]@Polygons[[1]]@coords)
colnames(boundaries_simplified) <- c("long", "lat")
boundaries_simplified <- boundaries_simplified %>% distinct()

# Simplified boundary plot
ggplot(boundaries_simplified, aes(long, lat)) +
  geom_point(color = "black") +
  coord_map(xlim = c(inner_leuven_coords_2$long[1], inner_leuven_coords_2$long[2]), 
            ylim = c(inner_leuven_coords_2$lat[1], inner_leuven_coords_2$lat[2]))

# Comparison

only_inner_leuven %>%
  ggplot() +
  geom_polygon(aes(long, lat, group = group, fill = info), color = "black", fill = "white") +
  coord_map(xlim = c(inner_leuven_coords_2$long[1], inner_leuven_coords_2$long[2]), 
            ylim = c(inner_leuven_coords_2$lat[1], inner_leuven_coords_2$lat[2]))

boundaries_simplified %>%
  ggplot() +
  geom_polygon(aes(long, lat),  color = "black", fill = "transparent") +
  coord_map(xlim = c(inner_leuven_coords_2$long[1], inner_leuven_coords_2$long[2]), 
            ylim = c(inner_leuven_coords_2$lat[1], inner_leuven_coords_2$lat[2]))

# ||||||||||||||||
# Data export ----
# ||||||||||||||||

boundary_nodes <- boundaries_simplified
leuven_box_coords <- inner_leuven_coords_2
leuven_complete <- only_inner_leuven
leuven_boundary <- boundaries
save(boundary_nodes, leuven_box_coords, leuven_complete, leuven_boundary, file = "data/map.RData")
