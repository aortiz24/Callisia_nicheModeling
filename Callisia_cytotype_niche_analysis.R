##running maxent and creating models with Callisia graminea cytotype data

# load libraries
library(dismo)
library(fields)
library(maps)
library(rgdal)
library(raster)
library(maptools)
library(dplyr)
library(rJava)

# import occurrence data and convert to format required by maxent
Callisia.both <- read.csv(file="taxaData/CallisiaCompletedData.csv") %>%
  select(Cytotype,Latitude,Longitude)
Callisia.both <- na.omit(Callisia.both)
diploid <- Callisia.both %>%
  filter(Cytotype=="2X")
diploid <- diploid[,c(3,2)]
tetraploid <- Callisia.both %>%
  filter(Cytotype=="4X")
tetraploid <- tetraploid[,c(3,2)]

# create stack of non-correlated layers (as determined by layerPrep.R)
predictors <- stack(alt, bio2, bio3, bio5, bio6, bio8, bio9, bio12, bio13, bio14, bio19) 
plot(predictors)

# run maxent for diploid (default parameters)
maxDip <- maxent(predictors, diploid)
maxDip # views results in browser window
response(maxDip) # show response curves for each layer
rDip <- predict(maxDip, predictors) # create model
plot(rDip)
points(diploid)
writeRaster(rDip, "models/diploid_Callisia.grd")

# run maxent for tetraploid (default parameters)
maxTetra <- maxent(predictors, tetraploid) # run with default parameters
maxTetra # views results in browser window
response(maxTetra) # show response curves for each layer
rTetra <- predict(maxTetra, predictors) # create model 
plot(rTetra)
points(tetraploid)
writeRaster(rTetra, "models/tetraploid_Callisia.grd")

# more complicated maxent modeling for diploids
maxAdvanced1 <- maxent(predictors, diploid, args=c("randomseed=true", "replicatetype=crossvalidate", "replicates=640", "-J")) # takes much longer, but includes cross validation, random seed, runs jacknife
maxAdvanced1 # views results in browser window
response(maxAdvanced1) # show response curves for each layer
rDip <- predict(maxAdvanced1, predictors) # create model
plot(rDip)
points(diploid)
writeRaster(rDip, "models/diploid_Callisia1.grd")


# more complicated maxent modeling for tetraploids
maxAdvanced2 <- maxent(predictors, tetraploid, args=c("randomseed=true", "replicatetype=crossvalidate", "replicates=640", "-J")) # takes much longer, but includes cross validation, random seed, runs jacknife
maxAdvanced2 # views results in browser window
response(maxAdvanced2) # show response curves for each layer
rTetra <- predict(maxAdvanced2, predictors) # create model 
plot(rTetra)
points(tetraploid)
writeRaster(rTetra, "models/tetraploid_Callisia2.grd")


##nicheOverlap analysis
# read in models (raster)
rDip <- raster("models/diploid.grd")
rTetra <- raster("models/tetraploid.grd")
#comparing niches
nicheOverlap(rDip, rTetra, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip, rTetra, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic
