## making maps and preparing shapefiles

## load libraries
library(fields)
library(dplyr)
library(dismo)
library(maptools) 

## load data with taxa in different R objects
# import occurrence data and convert to the format required by maxent
Callisia.both <- read.csv(file="CallisiaCompletedData.csv") %>%
  select(Cytotype,Latitude,Longitude)
Callisia.both <- na.omit(Callisia.both)
diploid <- Callisia.both %>%
  filter(Cytotype=="2X")
diploid <- diploid[,c(3,2)]
tetraploid <- Callisia.both %>%
  filter(Cytotype=="4X")
tetraploid <- tetraploid[,c(3,2)]

## quick and dirty plot on map (could also plot points first and add map)
US(xlim=c(-85,-77), ylim=c(26,37))
points(diploid$lon, diploid$lat, col='orange', pch=20)
points(tetraploid$lon, tetraploid$lat, col='red', pch=20)

## a slightly more refined map (using built-in state outlines)
southeast <- c("florida", "georgia", "north carolina", "south carolina")
map(database="state", regions = southeast, interior=T, lwd=2)
points(diploid$lon, diploid$lat, col='orange', pch=20, cex=2)
points(tetraploid$lon, tetraploid$lat, col='red', pch=20, cex=2)

## using US census shapefiles, save custom (best) shapefile for modeling later
# download, unzip all state shapefiles to new directory
download.file("http://www2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_us_state_20m.zip", "cb_2015_us_state_20m.zip")
dir.create("shapefiles")
unzip("cb_2015_us_state_20m.zip", exdir="shapefiles")
# load shapefiles and set projection
state <- readShapePoly("shapefiles/cb_2015_us_state_20m.shp") 
projection(state) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
# extract shapefiles of interest and save to file 
southeastStatesCap <- c("Florida", "Georgia", "North Carolina", "South Carolina")
SEstates <- state[as.character(state@data$NAME) %in% southeastStatesCap, ]
writeSpatialShape(SEstates, "shapefiles/SEstates")

# map using custom shapefile and save to file
dir.create("figures")
pdf(file="figures/mapping.pdf")
map(SEstates)
points(diploid$lon, diploid$lat, col='orange', pch=20, cex=2)
points(tetraploid$lon, tetraploid$lat, col='red', pch=20, cex=2)
dev.off()
