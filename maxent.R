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
both <- read.csv(file="taxaData/bothTaxa.csv")
diploid <- both %>%
  filter(cytotype=="2X")
diploid <- diploid[,c(3,2)]
tetraploid <- both %>%
  filter(cytotype=="4X")
tetraploid <- tetraploid[,c(3,2)]

# import layers with CRS specified
CRS <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
alt <- raster("layers/alt.asc", crs=CRS)
bio2 <- raster("layers/bio2.asc", crs=CRS)
bio3 <- raster("layers/bio3.asc", crs=CRS)
bio5 <- raster("layers/bio5.asc", crs=CRS)
bio6 <- raster("layers/bio6.asc", crs=CRS)
bio8 <- raster("layers/bio8.asc", crs=CRS)
bio9 <- raster("layers/bio9.asc", crs=CRS)
bio12 <- raster("layers/bio12.asc", crs=CRS)
bio13 <- raster("layers/bio13.asc", crs=CRS)
bio14 <- raster("layers/bio14.asc", crs=CRS)
bio19 <- raster("layers/bio19.asc", crs=CRS)

## create stack of non-correlated layers (as determined by layerPrep.R)
predictors <- stack(alt, bio2, bio3, bio5, bio6, bio8, bio9, bio12, bio13, bio14, bio19) 
# plot each layer individually
plot(predictors)

## basic bioclim modeling with diploids then tetraploids
# extract layer data for each point
dipPts <- extract(predictors, diploid)
# create bioclim model
dipBC <- bioclim(dipPts)
# predict bioclim model
dipBCpredict <- predict(predictors, dipBC)
# plot bioclim model
plot(dipBCpredict)

# extract layer data for each point
tetraPts <- extract(predictors, tetraploid)
# create bioclim model
tetraBC <- bioclim(tetraPts)
# predict bioclim model
tetraBCpredict <- predict(predictors, tetraBC)
# plot bioclim model
plot(tetraBCpredict)

## Default maxent modeling
# run maxent for diploid (default parameters for dismo)
maxDip <- maxent(predictors, diploid)
maxDip # views results in browser window
response(maxDip) # show response curves for each layer
rDip <- predict(maxDip, predictors) # create model
plot(rDip) # plot predictive model
points(diploid) # add points to predictive model
writeRaster(rDip, "models/diploid.grd")

# run maxent for tetraploid (default parameters for dismo)
maxTetra <- maxent(predictors, tetraploid) 
maxTetra 
response(maxTetra) 
rTetra <- predict(maxTetra, predictors) 
plot(rTetra)
points(tetraploid)
writeRaster(rTetra, "models/tetraploid.grd")

## Advanced modeling
# develop testing and training sets for diploid
fold <- kfold(diploid, k=5) #split occurence points into 5 sets
dipTest <- diploid[fold == 1, ] #take 20% (1/5) for testing
dipTrain <- diploid[fold != 1, ] #leave 40% for training
# fit training model for diploid
maxDipTrain <- maxent(predictors, dipTrain) #fit maxent model
maxDipTrain #view results in html
rDipTrain <- predict(maxDipTrain, predictors) #predict full model
plot(rDipTrain) #visualize full model
points(diploid) #add points to plot
# testing model for diploid
# extract background points
bg <- randomPoints(predictors, 1000)
# cross-validate model for diploid
maxDipTest <- evaluate(maxDipTrain, p=dipTest, a=bg, x=predictors)
maxDipTest #print results
threshold(maxDipTest) #identify threshold for presence or absence
plot(maxDipTest, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for diploid
pvtest <- data.frame(extract(predictors, dipTest))
avtest <- data.frame(extract(predictors, bg))
# cross-validate model
maxDipTest2 <- evaluate(maxDipTrain, p=pvtest, a=avtest)
maxDipTest2
# Alternative 2: predict to testing points for diploid
testp <- predict(maxDipTrain, pvtest)
testa <- predict(maxDipTrain, avtest)
maxDipTest3 <- evaluate(p=testp, a=testa)
maxDipTest3
# maxent with jackknife, random seed, and response curves, followed by cross-validation
maxDipAdv <- maxent(
  x=predictors,
  p=diploid,
  removeDuplicates=TRUE,
  nbg=10000,
  args=c(
    'randomseed=true', #default=false
    'threads=2', #default=1
    'responsecurves=true', #default=false
    'jackknife=true', #default=false
    'replicates=10', #default=1
    'replicatetype=crossvalidate',
    'maximumiterations=1000' #default=500
  )
)
maxDipAdv #view output as html

# develop testing and training sets for tetraploid
fold <- kfold(tetraploid, k=5) #split occurence points into 5 sets
tetraTest <- tetraploid[fold == 1, ] #take 20% (1/5) for testing
tetraTrain <- tetraploid[fold != 1, ] #leave 40% for training
# fit training model for tetraploid
maxTetraTrain <- maxent(predictors, tetraTrain) #fit maxent model
maxTetraTrain #view results in html
rTetraTrain <- predict(maxTetraTrain, predictors) #predict full model
plot(rTetraTrain) #visualize full model
points(tetraploid) #add points to plot
# testing model for tetraploid
# extract background points
bg <- randomPoints(predictors, 1000)
# cross-validate model for tetraploid
maxTetraTest <- evaluate(maxTetraTrain, p=tetraTest, a=bg, x=predictors)
maxTetraTest #print results
threshold(maxTetraTest) #identify threshold for presence or absence
plot(maxTetraTest, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for tetraploid
pvtest <- data.frame(extract(predictors, tetraTest))
avtest <- data.frame(extract(predictors, bg))
# cross-validate model
maxTetraTest2 <- evaluate(maxTetraTrain, p=pvtest, a=avtest)
maxTetraTest2
# Alternative 2: predict to testing points for tetraploid
testp <- predict(maxTetraTrain, pvtest)
testa <- predict(maxTetraTrain, avtest)
maxTetraTest3 <- evaluate(p=testp, a=testa)
maxTetraTest3
# maxent with jackknife, random seed, and response curves, followed by cross-validation
maxTetraAdv <- maxent(
  x=predictors,
  p=tetraploid,
  removeDuplicates=TRUE,
  nbg=10000,
  args=c(
    'randomseed=true', #default=false
    'threads=2', #default=1
    'responsecurves=true', #default=false
    'jackknife=true', #default=false
    'replicates=10', #default=1
    'replicatetype=crossvalidate',
    'maximumiterations=1000' #default=500
  )
)
maxTetraAdv #view output as html