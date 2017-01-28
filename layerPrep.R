## prepping layers for modeling

## load libraries (not all used directly, but required as dependencies)
library(dismo)
library(fields)
library(maps)
library(rgdal)
library(raster)
library(maptools)

# load previously created shapefile
SEstates <- readShapePoly("shapefiles/SEstates.shp") 

# if shapefile hasn't been created, use following code to download US Census states
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

###Past 1930
## load PRISM1930 ppt layers
ppt_4kmM2_0 <- raster("~data/dataLayers/PRISM1930/PRISM_ppt_stable_4kmM2_1930_all_bil/PRISM_ppt_stable_4kmM2_1930_bil.bil")

##change projection of past data layers
projection(ppt_4kmM2_0) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
ppt_0 <- mask(ppt_4kmM2_0, SEstates)
ppt.0 <- crop(ppt_0, extent(SEstates))
writeRaster(ppt.0, "~data/PastLayers/1930/ppt0.asc", format="ascii", overwrite=TRUE)

## load PRISM1930 tmax layers
tmax_4kmM2_0 <- raster("~data/dataLayers/PRISM1930/PRISM_tmax_stable_4kmM2_1930_all_bil/PRISM_tmax_stable_4kmM2_1930_bil.bil")

##change projection of past data layers
projection(tmax_4kmM2_0) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmax_0 <- mask(tmax_4kmM2_0, SEstates)
tmax.0 <- crop(tmax_0, extent(SEstates))
writeRaster(tmax.0, "/Users/hertwecklab/Documents/Thesis/data/PastLayers/1930/tmax0.asc", format="ascii", overwrite=TRUE)

## load PRISM1930 tmean layers
tmean_4kmM2_0 <- raster("/Users/hertwecklab/Documents/Thesis/data/dataLayers/PRISM1930/PRISM_tmean_stable_4kmM2_1930_all_bil/PRISM_tmean_stable_4kmM2_1930_bil.bil")

##change projection of past data layers
projection(tmean_4kmM2_0) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmean_0 <- mask(tmean_4kmM2_0, SEstates)
tmean.0 <- crop(tmean_0, extent(SEstates))
writeRaster(tmean.0, "/Users/hertwecklab/Documents/Thesis/data/PastLayers/1930/tmean0.asc", format="ascii", overwrite=TRUE)

## load PRISM1930 tmin layers
tmin_4kmM2_0 <- raster("~data/dataLayers/PRISM1930/PRISM_tmin_stable_4kmM2_1930_all_bil/PRISM_tmin_stable_4kmM2_1930_bil.bil")

##change projection of past data layers
projection(tmin_4kmM2_0) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmin_0 <- mask(tmin_4kmM2_0, SEstates)
tmin.0 <- crop(tmin_0, extent(SEstates))
writeRaster(tmin.0, "~data/PastLayers/1930/tmin0.asc", format="ascii", overwrite=TRUE)

## if layers have already been clipped, masked and saved and you need to reload them:
ppt0 <- raster("layers/avg_1930_ppt0.asc")
tmax0 <- raster("layers/avg_1930_tmax0.asc")
tmean0 <- raster("layers/avg_1930_tmean0.asc")
tmin0 <- raster("layers/avg_1930_tmin0.asc")

## correlation analysis
stack0 <- stack(tmin0, tmean0, tmax0, ppt0) 
corr0 <- layerStats(stack0, 'pearson', na.rm=TRUE)
c0 <- corr0$`pearson correlation coefficient`
write.csv(c0, "correlation1930.csv")
# inspect output for correlations between layers
#   0.7 and above (or -0.7 and below) are correlated
#   for this analysis, retain ppt0 

###Past 2014
## load PRISM2014 ppt layers
ppt_4kmM3_1 <- raster("~data/dataLayers/PRISM2014/PRISM_ppt_stable_4kmM3_2014_all_bil/PRISM_ppt_stable_4kmM3_2014_bil.bil")

##change projection of past data layers
projection(ppt_4kmM3_1) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
ppt_1 <- mask(ppt_4kmM3_1, SEstates)
ppt.1 <- crop(ppt_1, extent(SEstates))
writeRaster(ppt.1, "~data/PastLayers/2014/ppt1.asc", format="ascii", overwrite=TRUE)

## load PRISM2014 tmax layers
tmax_4kmM2_1 <- raster("~data/dataLayers/PRISM2014/PRISM_tmax_stable_4kmM2_2014_all_bil/PRISM_tmax_stable_4kmM2_2014_bil.bil")

##change projection of past data layers
projection(tmax_4kmM2_1) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmax_1 <- mask(tmax_4kmM2_1, SEstates)
tmax.1 <- crop(tmax_1, extent(SEstates))
writeRaster(tmax.1, "~data/PastLayers/2014/tmax1.asc", format="ascii", overwrite=TRUE)

## load PRISM2014 tmean layers
tmean_4kmM2_1 <- raster("~data/dataLayers/PRISM2014/PRISM_tmean_stable_4kmM2_2014_all_bil/PRISM_tmean_stable_4kmM2_2014_bil.bil")

##change projection of past data layers
projection(tmean_4kmM2_1) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmean_1 <- mask(tmean_4kmM2_1, SEstates)
tmean.1 <- crop(tmean_1, extent(SEstates))
writeRaster(tmean.1, "~data/PastLayers/2014/tmean1.asc", format="ascii", overwrite=TRUE)

## load PRISM2014 tmin layers
tmin_4kmM2_1 <- raster("~data/dataLayers/PRISM2014/PRISM_tmin_stable_4kmM2_2014_all_bil/PRISM_tmin_stable_4kmM2_2014_bil.bil")

##change projection of past data layers
projection(tmin_4kmM2_1) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmin_1 <- mask(tmin_4kmM2_1, SEstates)
tmin.1 <- crop(tmin_1, extent(SEstates))
writeRaster(tmin.1, "~data/PastLayers/2014/tmin1.asc", format="ascii", overwrite=TRUE)

## if layers have already been clipped, masked and saved and you need to reload them:
ppt1 <- raster("layers/avg_2014_ppt0.asc")
tmax1 <- raster("layers/avg_2014_tmax0.asc")
tmean1 <- raster("layers/avg_2014_tmean0.asc")
tmin1 <- raster("layers/avg_2014_tmin0.asc")

## correlation analysis
stack1 <- stack(tmin1, tmean1, tmax1, ppt1) 
corr1 <- layerStats(stack1, 'pearson', na.rm=TRUE)
c1 <- corr1$`pearson correlation coefficient`
write.csv(c1, "correlation2014.csv")
# inspect output for correlations between layers
#   0.7 and above (or -0.7 and below) are correlated
#   for this analysis, retain tmax1, ppt1
