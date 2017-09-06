## running maxent in R
# use ?maxent to check documentation for correct installation of maxent.jar

# load libraries
library(dismo)
library(fields)
library(maps)
library(rgdal)
library(raster)
library(maptools)
library(dplyr)
library(rJava)

# create directory for saving models later
dir.create("models")

## import occurrence data and convert to format required by maxent
Callisia.both <- read.csv(file="CallisiaCompletedData.csv") %>%
  dplyr::select(Cytotype,Latitude,Longitude)
Callisia.both <- na.omit(Callisia.both)
diploid <- Callisia.both %>%
  filter(Cytotype=="2X") 
diploid <- diploid[,c(3,2)]
tetraploid <- Callisia.both %>%
  filter(Cytotype=="4X")
tetraploid <- tetraploid[,c(3,2)]
#deleting rows whose points are outside of SEstates object
tetraploid<- tetraploid[-c(10,46,56), ]
#occurrence data for both Callisia cytotypes
Callisia.both <- read.csv(file="CallisiaCompletedData.csv") %>%
  dplyr::select(Cytotype,Latitude,Longitude)
both <- Callisia.both[,c(3,2)]
both <- na.omit(both)

#layers ending in 9 are for PRISM1929
#layers ending in 11 are for PRISM2011
# import layers with CRS specified
CRS <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ppt9 <- raster("layers/ppt9.asc", crs=CRS)
tmax9 <- raster("layers/tmax9.asc", crs=CRS)
tmean9 <- raster("layers/tmean9.asc", crs=CRS)
tmin9 <- raster("layers/tmin9.asc", crs=CRS)
vpdmax9 <- raster("layers/vpdmax9.asc", crs=CRS)
vpdmin9 <- raster("layers/vpdmin9.asc", crs=CRS)
tdmean9 <- raster("layers/tdmean9.asc", crs=CRS)
ppt11 <- raster("layers/ppt11.asc", crs=CRS)
tmax11 <- raster("layers/tmax11.asc", crs=CRS)
tmean11 <- raster("layers/tmean11.asc", crs=CRS)
tmin11 <- raster("layers/tmin11.asc", crs=CRS)
vpdmax11 <- raster("layers/vpdmax11.asc", crs=CRS)
vpdmin11 <- raster("layers/vpdmin11.asc", crs=CRS)
tdmean11 <- raster("layers/tdmean11.asc", crs=CRS)

## create stack of non-correlated layers (as determined by layerPrep.R)
predictors9<- stack(tmean9, ppt9, vpdmin9)
predictors11<- stack(tmean11, ppt11, vpdmax11, vpdmin11)

## Default maxent modeling
# run maxent for diploid (default parameters for dismo)
maxDip9 <- maxent(predictors9, diploid)
maxDip9 # views results in browser window
response(maxDip9) # show response curves for each layer
rDip9 <- predict(maxDip9, predictors9) # create model
plot(rDip9) # plot predictive model
#dev.copy2pdf(file="figures/MaxDiploid1929.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(diploid) # add points to predictive model
writeRaster(rDip9, "models/diploid1929.grd")

# run maxent for tetraploid (default parameters for dismo)
maxTetra9 <- maxent(predictors9, tetraploid) 
maxTetra9 
response(maxTetra9) 
rTetra9 <- predict(maxTetra9, predictors9) 
plot(rTetra9)
#dev.copy2pdf(file="figures/MaxTetraploid1929.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(tetraploid)
writeRaster(rTetra9, "models/tetraploid1929.grd")

# run maxent for both cytotypes (default parameters for dismo)
maxBoth9 <- maxent(predictors9, both)
maxBoth9 # views results in browser window
response(maxBoth9) # show response curves for each layer
rBoth9 <- predict(maxBoth9, predictors9) # create model
plot(rBoth9) # plot predictive model
#dev.copy2pdf(file="figures/MaxBoth1929.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(both) # add points to predictive model
writeRaster(rBoth9, "models/both1929.grd")

## Advanced modeling
# maxent with jackknife, random seed, and response curves, followed by cross
maxDipAdv9 <- maxent(
  x=predictors9,
  p=diploid,
  removeDuplicates=TRUE,
  nbg=10000,
  args=c(
    'randomseed=true', #default=false
    'threads=2', #default=1
    'responsecurves=true', #default=false
    'replicates=10', #default=1
    'replicatetype=crossvalidate',
    'maximumiterations=1000' #default=500
  )
)
maxDipAdv9 #view output as html
response(maxDipAdv9) # show response curves for each layer
rDipAdv9 <- predict(maxDipAdv9, predictors9) # create model
plot(rDipAdv9) # plot predictive model
#dev.copy2pdf(file="figures/AdvDiploid1929.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(diploid) # add points to predictive model
writeRaster(rDipAdv9, "models/diploidAdv1929.grd")

# maxent with jackknife, random seed, and response curves, followed by cross-validation
maxTetraAdv9 <- maxent(
  x=predictors9,
  p=tetraploid,
  removeDuplicates=TRUE,
  nbg=10000,
  args=c(
    'randomseed=true', #default=false
    'threads=2', #default=1
    'responsecurves=true', #default=false
    'replicates=10', #default=1
    'replicatetype=crossvalidate',
    'maximumiterations=1000' #default=500
  )
)
maxTetraAdv9 #view output as html
response(maxTetraAdv9) # show response curves for each layer
rTetraAdv9 <- predict(maxTetraAdv9, predictors9) # create model
plot(rTetraAdv9) # plot predictive model
#dev.copy2pdf(file="figures/AdvTetraploid1929.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(tetraploid) # add points to predictive model
writeRaster(rTetraAdv9, "models/tetraploidAdv1929.grd")

# maxent with jackknife, random seed, and response curves, followed by cross
maxBothAdv9 <- maxent(
  x=predictors9,
  p=both,
  removeDuplicates=TRUE,
  nbg=10000,
  args=c(
    'randomseed=true', #default=false
    'threads=2', #default=1
    'responsecurves=true', #default=false
    'replicates=10', #default=1
    'replicatetype=crossvalidate',
    'maximumiterations=1000' #default=500
  )
)
maxBothAdv9 #view output as html
response(maxBothAdv9) # show response curves for each layer
rBothAdv9 <- predict(maxBothAdv9, predictors9) # create model
plot(rBothAdv9) # plot predictive model
#dev.copy2pdf(file="figures/AdvBoth1929.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(tetraploid) # add points to predictive model
writeRaster(rBothAdv9, "models/bothAdv1929.grd")

## Default maxent modeling
# run maxent for diploid (default parameters for dismo)
maxDip11 <- maxent(predictors11, diploid)
maxDip11 # views results in browser window
response(maxDip11) # show response curves for each layer
rDip11 <- predict(maxDip11, predictors11) # create model
plot(rDip11) # plot predictive model
#dev.copy2pdf(file="figures/MaxDiploid2011.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(diploid) # add points to predictive model
writeRaster(rDip11, "models/diploid2011.grd")

# run maxent for tetraploid (default parameters for dismo)
maxTetra11 <- maxent(predictors11, tetraploid) 
maxTetra11 
response(maxTetra11) 
rTetra11 <- predict(maxTetra11, predictors11) 
plot(rTetra11)
#dev.copy2pdf(file="figures/MaxTetraploid2011.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(tetraploid)
writeRaster(rTetra11, "models/tetraploid2011.grd")

# run maxent for both cytotypes (default parameters for dismo)
maxBoth11 <- maxent(predictors11, both)
maxBoth11 # views results in browser window
response(maxBoth11) # show response curves for each layer
rBoth11 <- predict(maxBoth11, predictors11) # create model
plot(rBoth11) # plot predictive model
#dev.copy2pdf(file="figures/MaxBoth2011.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(both) # add points to predictive model
writeRaster(rBoth11, "models/both2011.grd")

## Advanced modeling
# maxent with jackknife, random seed, and response curves, followed by cross
maxDipAdv11 <- maxent(
  x=predictors11,
  p=diploid,
  removeDuplicates=TRUE,
  nbg=10000,
  args=c(
    'randomseed=true', #default=false
    'threads=2', #default=1
    'responsecurves=true', #default=false
    'replicates=10', #default=1
    'replicatetype=crossvalidate',
    'maximumiterations=1000' #default=500
  )
)
maxDipAdv11 #view output as html
response(maxDipAdv11) # show response curves for each layer
rDipAdv11 <- predict(maxDipAdv11, predictors11) # create model
plot(rDipAdv11) # plot predictive model
#dev.copy2pdf(file="figures/AdvDiploid2011.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(diploid) # add points to predictive model
writeRaster(rDipAdv11, "models/diploidAdv2011.grd")

# maxent with jackknife, random seed, and response curves, followed by cross-validation
maxTetraAdv11 <- maxent(
  x=predictors11,
  p=tetraploid,
  removeDuplicates=TRUE,
  nbg=10000,
  args=c(
    'randomseed=true', #default=false
    'threads=2', #default=1
    'responsecurves=true', #default=false
    'replicates=10', #default=1
    'replicatetype=crossvalidate',
    'maximumiterations=1000' #default=500
  )
)
maxTetraAdv11 #view output as html
response(maxTetraAdv11) # show response curves for each layer
rTetraAdv11 <- predict(maxTetraAdv11, predictors11) # create model
plot(rTetraAdv11) # plot predictive model
#dev.copy2pdf(file="figures/AdvTetraploid2011.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(tetraploid) # add points to predictive model
writeRaster(rTetraAdv11, "models/tetraploidAdv2011.grd")

# maxent with jackknife, random seed, and response curves, followed by cross
maxBothAdv11 <- maxent(
  x=predictors11,
  p=both,
  removeDuplicates=TRUE,
  nbg=10000,
  args=c(
    'randomseed=true', #default=false
    'threads=2', #default=1
    'responsecurves=true', #default=false
    'replicates=10', #default=1
    'replicatetype=crossvalidate',
    'maximumiterations=1000' #default=500
  )
)
maxBothAdv11 #view output as html
response(maxBothAdv11) # show response curves for each layer
rBothAdv11 <- predict(maxBothAdv11, predictors11) # create model
plot(rBothAdv11) # plot predictive model
#dev.copy2pdf(file="figures/AdvBoth2011.pdf", width = 7, height = 5) #save plot as pdf to the figures directory
points(tetraploid) # add points to predictive model
writeRaster(rBothAdv11, "models/bothAdv2011.grd")
