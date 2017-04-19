## prepping layers for modeling on Ubuntu computer

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

###Past 1929
## load PRISM1929 ppt layers
ppt_4kmM2_9 <- raster("~/Documents/Thesis/data/dataLayers/PRISM1929/PRISM_ppt_stable_4kmM2_1929_all_bil/PRISM_ppt_stable_4kmM2_1929_bil.bil")

##change projection of past data layers
projection(ppt_4kmM2_9) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
ppt_9 <- mask(ppt_4kmM2_9, SEstates)
ppt.9 <- crop(ppt_9, extent(SEstates))
writeRaster(ppt.9, "~/Desktop/Callisia_nicheModeling/layers/ppt9.asc", format="ascii", overwrite=TRUE)

## load PRISM1929 tmax layers
tmax_4kmM2_9 <- raster("~/Documents/Thesis/data/dataLayers/PRISM1929/PRISM_tmax_stable_4kmM2_1929_all_bil/PRISM_tmax_stable_4kmM2_1929_bil.bil")

##change projection of past data layers
projection(tmax_4kmM2_9) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmax_9 <- mask(tmax_4kmM2_9, SEstates)
tmax.9 <- crop(tmax_9, extent(SEstates))
writeRaster(tmax.9, "~/Desktop/Callisia_nicheModeling/layers/tmax9.asc", format="ascii", overwrite=TRUE)

## load PRISM1929 tmean layers
tmean_4kmM2_9 <- raster("~/Documents/Thesis/data/dataLayers/PRISM1929/PRISM_tmean_stable_4kmM2_1929_all_bil/PRISM_tmean_stable_4kmM2_1929_bil.bil")

##change projection of past data layers
projection(tmean_4kmM2_9) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmean_9 <- mask(tmean_4kmM2_9, SEstates)
tmean.9 <- crop(tmean_9, extent(SEstates))
writeRaster(tmean.9, "~/Desktop/Callisia_nicheModeling/layers/tmean9.asc", format="ascii", overwrite=TRUE)

## load PRISM1929 tmin layers
tmin_4kmM2_9 <- raster("~/Documents/Thesis/data/dataLayers/PRISM1929/PRISM_tmin_stable_4kmM2_1929_all_bil/PRISM_tmin_stable_4kmM2_1929_bil.bil")

##change projection of past data layers
projection(tmin_4kmM2_9) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmin_9 <- mask(tmin_4kmM2_9, SEstates)
tmin.9 <- crop(tmin_9, extent(SEstates))
writeRaster(tmin.9, "~/Desktop/Callisia_nicheModeling/layers/tmin9.asc", format="ascii", overwrite=TRUE)

## load PRISM1929 vpdmax layers
vpdmax_4kmM1_9 <- raster("~/Documents/Thesis/data/dataLayers/PRISM1929/PRISM_vpdmax_stable_4kmM1_1929_all_bil/PRISM_vpdmax_stable_4kmM1_1929_bil.bil")

##change projection of past data layers
projection(vpdmax_4kmM1_9) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
vpdmax_9 <- mask(vpdmax_4kmM1_9, SEstates)
vpdmax.9 <- crop(vpdmax_9, extent(SEstates))
writeRaster(vpdmax.9, "~/Desktop/Callisia_nicheModeling/layers/vpdmax9.asc", format="ascii", overwrite=TRUE)

## load PRISM1929 vpdmin layers
vpdmin_4kmM1_9 <- raster("~/Documents/Thesis/data/dataLayers/PRISM1929/PRISM_vpdmin_stable_4kmM1_1929_all_bil/PRISM_vpdmin_stable_4kmM1_1929_bil.bil")

##change projection of past data layers
projection(vpdmin_4kmM1_9) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
vpdmin_9 <- mask(vpdmin_4kmM1_9, SEstates)
vpdmin.9 <- crop(vpdmin_9, extent(SEstates))
writeRaster(vpdmin.9, "~/Desktop/Callisia_nicheModeling/layers/vpdmin9.asc", format="ascii", overwrite=TRUE)

## load PRISM1929 tdmean layers
tdmean_4kmM1_9 <- raster("~/Documents/Thesis/data/dataLayers/PRISM1929/PRISM_tdmean_stable_4kmM1_1929_all_bil/PRISM_tdmean_stable_4kmM1_1929_bil.bil")

##change projection of past data layers
projection(tdmean_4kmM1_9) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tdmean_9 <- mask(tdmean_4kmM1_9, SEstates)
tdmean.9 <- crop(tdmean_9, extent(SEstates))
writeRaster(tdmean.9, "~/Desktop/Callisia_nicheModeling/layers/tdmean9.asc", format="ascii", overwrite=TRUE)

## if layers have already been clipped, masked and saved and you need to reload them:
ppt9 <- raster("layers/ppt9.asc")
tmax9 <- raster("layers/tmax9.asc")
tmean9 <- raster("layers/tmean9.asc")
tmin9 <- raster("layers/tmin9.asc")
vpdmax9 <- raster("layers/vpdmax9.asc")
vpdmin9 <- raster("layers/vpdmin9.asc")
tdmean9 <- raster("layers/tdmean9.asc")

## correlation analysis
stack9 <- stack(tmin9, tmean9, tmax9, ppt9, vpdmax9, vpdmin9, tdmean9) 
corr9 <- layerStats(stack9, 'pearson', na.rm=TRUE)
c9 <- corr9$`pearson correlation coefficient`
write.csv(c9, "correlation1929.csv")
# inspect output for correlations between layers
# absolute value of 0.7 or greater are correlated
# for this analysis, retain tmean9, ppt9, vpdmin9

###Past 2011
## load PRISM2011 ppt layers
ppt_4kmM3_11 <- raster("~/Documents/Thesis/data/dataLayers/PRISM2011/PRISM_ppt_stable_4kmM3_2011_all_bil/PRISM_ppt_stable_4kmM3_2011_bil.bil")

##change projection of past data layers
projection(ppt_4kmM3_11) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
ppt_11 <- mask(ppt_4kmM3_11, SEstates)
ppt.11 <- crop(ppt_11, extent(SEstates))
writeRaster(ppt.11, "~/Documents/Thesis/data/PastLayers_clipped/2011/ppt11.asc", format="ascii", overwrite=TRUE)

## load PRISM2011 tmax layers
tmax_4kmM2_11 <- raster("~/Documents/Thesis/data/dataLayers/PRISM2011/PRISM_tmax_stable_4kmM2_2011_all_bil/PRISM_tmax_stable_4kmM2_2011_bil.bil")

##change projection of past data layers
projection(tmax_4kmM2_11) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmax_11 <- mask(tmax_4kmM2_11, SEstates)
tmax.11 <- crop(tmax_11, extent(SEstates))
writeRaster(tmax.11, "~/Documents/Thesis/data/PastLayers_clipped/2011/tmax11.asc", format="ascii", overwrite=TRUE)

## load PRISM2011 tmean layers
tmean_4kmM2_11 <- raster("~/Documents/Thesis/data/dataLayers/PRISM2011/PRISM_tmean_stable_4kmM2_2011_all_bil/PRISM_tmean_stable_4kmM2_2011_bil.bil")

##change projection of past data layers
projection(tmean_4kmM2_11) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmean_11 <- mask(tmean_4kmM2_11, SEstates)
tmean.11 <- crop(tmean_11, extent(SEstates))
writeRaster(tmean.11, "~/Documents/Thesis/data/PastLayers_clipped/2011/tmean11.asc", format="ascii", overwrite=TRUE)

## load PRISM2011 tmin layers
tmin_4kmM2_11 <- raster("~/Documents/Thesis/data/dataLayers/PRISM2011/PRISM_tmin_stable_4kmM2_2011_all_bil/PRISM_tmin_stable_4kmM2_2011_bil.bil")

##change projection of past data layers
projection(tmin_4kmM2_11) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tmin_11 <- mask(tmin_4kmM2_11, SEstates)
tmin.11 <- crop(tmin_11, extent(SEstates))
writeRaster(tmin.11, "~/Documents/Thesis/data/PastLayers_clipped/2011/tmin11.asc", format="ascii", overwrite=TRUE)

## load PRISM2011 vpdmax layers
vpdmax_4kmM1_11 <- raster("~/Documents/Thesis/data/dataLayers/PRISM2011/PRISM_vpdmax_stable_4kmM1_2011_all_bil/PRISM_vpdmax_stable_4kmM1_2011_bil.bil")

##change projection of past data layers
projection(vpdmax_4kmM1_11) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
vpdmax_11 <- mask(vpdmax_4kmM1_11, SEstates)
vpdmax.11 <- crop(vpdmax_11, extent(SEstates))
writeRaster(vpdmax.11, "~/Desktop/Callisia_nicheModeling/layers/vpdmax11.asc", format="ascii", overwrite=TRUE)

## load PRISM2011 vpdmin layers
vpdmin_4kmM1_11 <- raster("~/Documents/Thesis/data/dataLayers/PRISM2011/PRISM_vpdmin_stable_4kmM1_2011_all_bil/PRISM_vpdmin_stable_4kmM1_2011_bil.bil")

##change projection of past data layers
projection(vpdmin_4kmM1_11) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
vpdmin_11 <- mask(vpdmin_4kmM1_11, SEstates)
vpdmin.11 <- crop(vpdmin_11, extent(SEstates))
writeRaster(vpdmin.11, "~/Desktop/Callisia_nicheModeling/layers/vpdmin11.asc", format="ascii", overwrite=TRUE)

## load PRISM2011 tdmean layers
tdmean_4kmM1_11 <- raster("~/Documents/Thesis/data/dataLayers/PRISM2011/PRISM_tdmean_stable_4kmM1_2011_all_bil/PRISM_tdmean_stable_4kmM1_2011_bil.bil")

##change projection of past data layers
projection(tdmean_4kmM1_11) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project

## clip data layers
tdmean_11 <- mask(tdmean_4kmM1_11, SEstates)
tdmean.11 <- crop(tdmean_11, extent(SEstates))
writeRaster(tdmean.11, "~/Desktop/Callisia_nicheModeling/layers/tdmean11.asc", format="ascii", overwrite=TRUE)

## if layers have already been clipped, masked and saved and you need to reload them:
ppt11 <- raster("layers/ppt11.asc")
tmax11 <- raster("layers/tmax11.asc")
tmean11 <- raster("layers/tmean11.asc")
tmin11 <- raster("layers/tmin11.asc")
vpdmax11 <- raster("layers/vpdmax11.asc")
vpdmin11 <- raster("layers/vpdmin11.asc")
tdmean11 <- raster("layers/tdmean11.asc")

## correlation analysis
stack11 <- stack(tmin11, tmean11, tmax11, ppt11, vpdmax11, vpdmin11, tdmean11) 
corr11 <- layerStats(stack11, 'pearson', na.rm=TRUE)
c11 <- corr11$`pearson correlation coefficient`
write.csv(c11, "correlation2011.csv")
# inspect output for correlations between layers
# absolute value of 0.7 or greater are correlated
# for this analysis, retain tmean11, ppt11, vpdmax11, vpdmin11