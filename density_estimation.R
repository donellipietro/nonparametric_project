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

# install.packages('mapproj')
# install.packages("deldir")

library(fdaPDE)
library(ggplot2)
require(rgeos)
require(sp)
library(tidyverse)
library(mapproj)
library(deldir)
library(pracma)


# ||||||||||||||
# Functions ----
# ||||||||||||||

filterByProximity <- function(xy, dist, mapUnits = F) {
  #xy can be either a SpatialPoints or SPDF object, or a matrix
  #dist is in km if mapUnits=F, in mapUnits otherwise
  if (!mapUnits) {
    d <- spDists(xy,longlat=T)
  }
  if (mapUnits) {
    d <- spDists(xy,longlat=F)
  }
  diag(d) <- NA
  close <- (d <= dist)
  diag(close) <- NA
  closePts <- which(close,arr.ind=T)
  discard <- matrix(nrow=2,ncol=2)
  if (nrow(closePts) > 0) {
    while (nrow(closePts) > 0) {
      if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
        discard <- rbind(discard, closePts[1,])
        closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
      }
    }
    discard <- discard[complete.cases(discard),]
    return(xy[-discard[,1],])
  }
  if (nrow(closePts) == 0) {
    return(xy)
  }
}


# |||||||||
# Data ----
# |||||||||

load("data/map.RData")
load("data/stops.RData")


# ||||||||||||||||||||
# ggplot settings ----
# ||||||||||||||||||||

theme_set(theme_bw())


# ||||||||||||||||||
# Standard Mesh ----
# ||||||||||||||||||

# Mesh Creation
# |||||||||||||

n <- dim(boundary_nodes)[1]
segments <- data.frame(1:n, c(2:n, 1))

nodes <- boundary_nodes[,c("long", "lat")]
mesh <- create.mesh.2D(nodes = nodes,  segments = segments)
mesh <- refine.mesh.2D(mesh, delaunay = TRUE, maximum_area = 0.000008, minimum_angle = 30)

# plot(mesh, xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
#      ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2]))

n.nodes <- dim(mesh$nodes)[1]
n.nodes

n.data <- dim(stops)[1]
n.data

# Plot mesh
# |||||||||

triangles <- data.frame(mesh$triangles)
colnames(triangles) <- c("x1", "x2", "x3")
triangles$id <- 1:dim(triangles)[1]
triangles <- reshape(triangles, 
                     idvar = "id", 
                     varying = c("x1", "x2", "x3"),
                     v.names = "index", 
                     timevar = "order",
                     times = c(1:3),
                     direction = "long")

nodes <- data.frame(mesh$nodes)
colnames(nodes) <- c("long", "lat")

triangles_cordinates <- nodes[triangles$index,]
triangles = cbind(triangles, triangles_cordinates)


# FEM basis
# |||||||||

FEMbasis <- create.FEM.basis(mesh)


# # lambda selection through k-CV
# # |||||||||||||||||||||||||||||
# 
# lambda <- 10^(seq(-3, -10, by = -0.5))
# 
# log.density <- DE.FEM(data = stops[,c(2,1)], FEMbasis = FEMbasis,
#                       lambda = lambda,
#                       preprocess_method = "RightCV", nfolds = 40,
#                       nsim = 10000, heatStep=0.1, heatIter=500,
#                       step_method="Fixed_Step", direction_method="BFGS", tol1=1e-3)
# 
# # CV error 
# 
# cv_error <- log.density$CV_err
# cv_error[cv_error == -Inf] = 0
# cv_error.min <- min(cv_error)
# cv_error.index_min <- which.min(cv_error)
# lambda.opt = lambda[cv_error.index_min]
# plot(log10(lambda), cv_error, ylab = "CV error", xlab = expression(log10(lambda)))
# points(log10(lambda.opt), cv_error.min, col = 'red', pch = 7)
# 
# lambda.opt


# Density estimation
# ||||||||||||||||||

log.density <- DE.FEM(data = stops[,c(2,1)], FEMbasis = FEMbasis,
                      lambda = 10^(-5.5),
                      # preprocess_method = "RightCV", nfolds = 5,
                      nsim = 10000, heatStep=0.1, heatIter=500,
                      step_method="Fixed_Step", direction_method="BFGS", tol1=1e-3)

FEMobj <- FEM(coeff=exp(log.density$g), FEMbasis = FEMbasis)



## Plots for presentation ----
## |||||||||||||||||||||||||||

# Leuven complete

leuven_complete %>%
  ggplot() +
  geom_polygon(aes(long, lat, group = group, fill = info), color = "black", fill = "white", linewidth = 0.2) +
  geom_polygon(data = leuven_boundary, aes(long, lat), color = "black", fill = "transparent", linewidth = 1) +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Leuven borders

boundary_nodes %>%
  ggplot() +
  geom_polygon(aes(long, lat), color = "black", fill = "transparent", linewidth = 1) +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Leuven + stops

leuven_complete %>%
  ggplot() +
  geom_polygon(aes(long, lat, group = group, fill = info), color = "black", fill = "white", linewidth = 0.5) +
  geom_polygon(data = leuven_boundary, aes(long, lat), color = "black", fill = "transparent", linewidth = 1) +
  geom_point(data = stops, aes(x = long, y = lat), color = "blue") +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Mesh

boundary_nodes %>%
  ggplot() +
  geom_polygon(data = triangles, aes(long, lat, group = id), color = "black", fill = "white", linewidth = 0.2) +
  geom_polygon(aes(long, lat), color = "black", fill = "transparent", linewidth = 1) +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Mesh + stops

boundary_nodes %>%
  ggplot() +
  geom_polygon(data = triangles, aes(long, lat, group = id), color = "black", fill = "white", linewidth = 0.2) +
  geom_polygon(aes(long, lat), color = "black", fill = "transparent", linewidth = 1) +
  geom_point(data = stops, aes(x = long, y = lat), color = "blue") +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# 3d plot

plot(FEMobj)

# Heatmap

X = seq(leuven_box_coords$long[1], leuven_box_coords$long[2], length = 100)
Y = seq(leuven_box_coords$lat[1], leuven_box_coords$lat[2], length = 300)

XY = meshgrid(X, Y)
XX = as.vector(XY[[1]])
YY = as.vector(XY[[2]])
grid = data.frame(long = XX, lat = YY)

Z = eval.FEM(FEMobj, locations = grid)

density <- data.frame(long = XX, lat = YY, Z = Z)
density <- density %>% drop_na()

density %>%
  ggplot(aes(long, lat, fill= Z)) + 
  geom_tile() +
  geom_polygon(data = leuven_complete, aes(long, lat, group = group, fill = info), color = "black", fill = "transparent", linewidth = 0.2) +
  scale_fill_gradient(low = "red", high = "yellow") +
  geom_polygon(data = boundary_nodes, aes(long, lat, fill = info), color = "black", fill = "transparent", linewidth = 1) +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# ||||||||||||||||||
# Adaptive Mesh ----
# ||||||||||||||||||

# Filter By Proximity
pts  <- as.matrix(stops[,c(1,2)])
stops_filtered <- data.frame(filterByProximity(pts,dist=0.002, mapUnits=T))

leuven_complete %>%
  ggplot() +
  geom_polygon(aes(long, lat, group = group), color = "black", fill = "white") +
  geom_point(data = stops, aes(x = long, y = lat), color = "red") +
  geom_point(data = stops_filtered, aes(x = long, y = lat), color = "green") +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]),
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2]))


# Calculate Voronoi Tesselation and tiles
tesselation <- deldir(stops_filtered$long, stops_filtered$lat,
                      rw = c(leuven_box_coords$long,leuven_box_coords$lat))
tiles <- tile.list(tesselation)

# Internal Voronoi Tesselation vertices
internal = point.in.polygon(tesselation$dirsgs[,1],
                                     tesselation$dirsgs[,2],
                                     leuven_boundary$long,
                                     leuven_boundary$lat)
internal_vertices <- tesselation$dirsgs[internal==1, 1:2]
colnames(internal_vertices) = c("long", "lat")

internal_vertices <- data.frame(filterByProximity(as.matrix(internal_vertices), dist=0.002, mapUnits=T))

# plot(tiles, pch = 19,
#      col.pts = "white",
#      border = "white",
#      fillcol = hcl.colors(50, "viridis"),
#      clipp = list(x = leuven_boundary$long, y = leuven_boundary$lat))
# points(internal_vertices, pch='x')


n <- dim(boundary_nodes)[1]
segments <- data.frame(1:n, c(2:n, 1))
nodes <- rbind(boundary_nodes[,c("long", "lat")], internal_vertices)
nodes <- nodes %>% distinct()

mesh <- create.mesh.2D(nodes = nodes,  segments = segments)
mesh <- refine.mesh.2D(mesh, delaunay = TRUE, minimum_angle = 30)

# plot(mesh, xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
#      ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2]))

n.nodes <- dim(mesh$nodes)[1]
n.nodes

n.data <- dim(stops)[1]
n.data

# Plot mesh
# |||||||||

triangles <- data.frame(mesh$triangles)
colnames(triangles) <- c("x1", "x2", "x3")
triangles$id <- 1:dim(triangles)[1]
triangles <- reshape(triangles, 
                     idvar = "id", 
                     varying = c("x1", "x2", "x3"),
                     v.names = "index", 
                     timevar = "order",
                     times = c(1:3),
                     direction = "long")

nodes <- data.frame(mesh$nodes)
colnames(nodes) <- c("long", "lat")

triangles_cordinates <- nodes[triangles$index,]
triangles = cbind(triangles, triangles_cordinates)


# FEM basis
# |||||||||

FEMbasis <- create.FEM.basis(mesh)


# lambda selection through k-CV
# |||||||||||||||||||||||||||||

# lambda <- 10^(seq(-3, -10, by = -0.2))
# 
# log.density <- DE.FEM(data = stops[,c(2,1)], FEMbasis = FEMbasis,
#                       lambda = lambda,
#                       preprocess_method = "RightCV", nfolds = 15,
#                       nsim = 10000, heatStep=0.1, heatIter=500,
#                       step_method="Fixed_Step", direction_method="BFGS", tol1=1e-3)

# CV error 

# cv_error <- log.density$CV_err
# cv_error[cv_error == -Inf] = 0
# cv_error.min <- min(cv_error)
# cv_error.index_min <- which.min(cv_error)
# plot(cv_error)
# points(cv_error.index_min, cv_error.min, col = 'red', pch = 7)
# 
# lambda.opt = lambda[cv_error.index_min]


# Density estimation
# ||||||||||||||||||

log.density <- DE.FEM(data = stops[,c(2,1)], FEMbasis = FEMbasis,
                      lambda = 10^(-5.86),
                      # preprocess_method = "RightCV", nfolds = 5,
                      nsim = 10000, heatStep=0.1, heatIter=500,
                      step_method="Fixed_Step", direction_method="BFGS", tol1=1e-3)

FEMobj <- FEM(coeff=exp(log.density$g), FEMbasis = FEMbasis)


## Plots for presentation ----
## |||||||||||||||||||||||||||

# Mesh

boundary_nodes %>%
  ggplot() +
  geom_polygon(data = triangles, aes(long, lat, group = id), color = "black", fill = "white", linewidth = 0.2) +
  geom_polygon(aes(long, lat), color = "black", fill = "transparent", linewidth = 1) +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Mesh + stops

boundary_nodes %>%
  ggplot() +
  geom_polygon(data = triangles, aes(long, lat, group = id), color = "black", fill = "white", linewidth = 0.2) +
  geom_polygon(aes(long, lat), color = "black", fill = "transparent", linewidth = 1) +
  geom_point(data = stops, aes(x = long, y = lat), color = "blue") +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# 3d plot

plot(FEMobj)

# Heatmap

X = seq(leuven_box_coords$long[1], leuven_box_coords$long[2], length = 100)
Y = seq(leuven_box_coords$lat[1], leuven_box_coords$lat[2], length = 300)

XY = meshgrid(X, Y)
XX = as.vector(XY[[1]])
YY = as.vector(XY[[2]])
grid = data.frame(long = XX, lat = YY)

Z = eval.FEM(FEMobj, locations = grid)

density <- data.frame(long = XX, lat = YY, Z = Z)
density <- density %>% drop_na()

density %>%
  ggplot(aes(long, lat, fill= Z)) + 
  geom_tile() +
  geom_polygon(data = leuven_complete, aes(long, lat, group = group, fill = info), color = "black", fill = "transparent", linewidth = 0.2) +
  scale_fill_gradient(low = "red", high = "yellow") +
  geom_polygon(data = boundary_nodes, aes(long, lat, fill = info), color = "black", fill = "transparent", linewidth = 1) +
  coord_map(xlim = c(leuven_box_coords$long[1], leuven_box_coords$long[2]), 
            ylim = c(leuven_box_coords$lat[1], leuven_box_coords$lat[2])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())