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
  select(Cytotype,Latitude,Longitude)
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
  select(Cytotype,Latitude,Longitude)
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

# plot each layer individually
plot(predictors9)
plot(predictors11)

### basic bioclim modeling with PRISM 1929 layers for diploids then tetraploids
# extract layer data for each point
dipPts9 <- extract(predictors9, diploid)
# create bioclim model
dipBC9 <- bioclim(dipPts9)
# predict bioclim model
dipBCpredict9 <- predict(predictors9, dipBC9)
# plot bioclim model
plot(dipBCpredict9)

# extract layer data for each point
tetraPts9 <- extract(predictors9, tetraploid)
# create bioclim model
tetraBC9 <- bioclim(tetraPts9)
# predict bioclim model
tetraBCpredict9 <- predict(predictors9, tetraBC9)
# plot bioclim model
plot(tetraBCpredict9)

# extract layer data for each point for both cytotypes
bothPts9 <- extract(predictors9, both)
# create bioclim model
bothBC9 <- bioclim(bothPts9)
# predict bioclim model
bothBCpredict9 <- predict(predictors9, bothBC9)
# plot bioclim model
plot(bothBCpredict9)

## Default maxent modeling
# run maxent for diploid (default parameters for dismo)
maxDip9 <- maxent(predictors9, diploid)
maxDip9 # views results in browser window
response(maxDip9) # show response curves for each layer
rDip9 <- predict(maxDip9, predictors9) # create model
plot(rDip9) # plot predictive model
points(diploid) # add points to predictive model
writeRaster(rDip9, "models/diploid1929.grd")

# run maxent for tetraploid (default parameters for dismo)
maxTetra9 <- maxent(predictors9, tetraploid) 
maxTetra9 
response(maxTetra9) 
rTetra9 <- predict(maxTetra9, predictors9) 
plot(rTetra9)
points(tetraploid)
writeRaster(rTetra9, "models/tetraploid1929.grd")

# run maxent for both cytotypes (default parameters for dismo)
maxBoth9 <- maxent(predictors9, both)
maxBoth9 # views results in browser window
response(maxBoth9) # show response curves for each layer
rBoth9 <- predict(maxBoth9, predictors9) # create model
plot(rBoth9) # plot predictive model
points(both) # add points to predictive model
writeRaster(rBoth9, "models/both1929.grd")

## Advanced modeling
# develop testing and training sets for diploid
fold <- kfold(diploid, k=5) #split occurence points into 5 sets
dipTest9 <- diploid[fold == 1, ] #take 20% (1/5) for testing
dipTrain9 <- diploid[fold != 1, ] #leave 40% for training
# fit training model for diploid
maxDipTrain9 <- maxent(predictors9, dipTrain9) #fit maxent model
maxDipTrain9 #view results in html
rDipTrain9 <- predict(maxDipTrain9, predictors9) #predict full model
plot(rDipTrain9) #visualize full model
points(diploid) #add points to plot
# testing model for diploid
# extract background points
bg9 <- randomPoints(predictors9, 1000)
# cross-validate model for diploid
maxDipTest9 <- evaluate(maxDipTrain9, p=dipTest9, a=bg9, x=predictors9)
maxDipTest9 #print results
threshold(maxDipTest9) #identify threshold for presence or absence
plot(maxDipTest9, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for diploid
pvtest9 <- data.frame(extract(predictors9, dipTest9))
avtest9 <- data.frame(extract(predictors9, bg9))
# cross-validate model
maxDipTest29 <- evaluate(maxDipTrain9, p=pvtest9, a=avtest9)
maxDipTest29
# Alternative 2: predict to testing points for diploid
testp9 <- predict(maxDipTrain9, pvtest9)
testa9 <- predict(maxDipTrain9, avtest9)
maxDipTest39 <- evaluate(p=testp9, a=testa9)
maxDipTest39
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
points(diploid) # add points to predictive model
writeRaster(rDipAdv9, "models/diploidAdv1929.grd")

# develop testing and training sets for tetraploid
fold <- kfold(tetraploid, k=5) #split occurence points into 5 sets
tetraTest9 <- tetraploid[fold == 1, ] #take 20% (1/5) for testing
tetraTrain9 <- tetraploid[fold != 1, ] #leave 40% for training
# fit training model for tetraploid
maxTetraTrain9 <- maxent(predictors9, tetraTrain9) #fit maxent model
maxTetraTrain9 #view results in html
rTetraTrain9 <- predict(maxTetraTrain9, predictors9) #predict full model
plot(rTetraTrain9) #visualize full model
points(tetraploid) #add points to plot
# testing model for tetraploid
# extract background points
bg9 <- randomPoints(predictors9, 1000)
# cross-validate model for tetraploid
maxTetraTest9 <- evaluate(maxTetraTrain9, p=tetraTest9, a=bg9, x=predictors9)
maxTetraTest9 #print results
threshold(maxTetraTest9) #identify threshold for presence or absence
plot(maxTetraTest9, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for tetraploid
pvtest9 <- data.frame(extract(predictors9, tetraTest9))
avtest9 <- data.frame(extract(predictors9, bg9))
# cross-validate model
maxTetraTest29 <- evaluate(maxTetraTrain9, p=pvtest9, a=avtest9)
maxTetraTest29
# Alternative 2: predict to testing points for tetraploid
testp9 <- predict(maxTetraTrain9, pvtest9)
testa9 <- predict(maxTetraTrain9, avtest9)
maxTetraTest39 <- evaluate(p=testp9, a=testa9)
maxTetraTest39
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
points(tetraploid) # add points to predictive model
writeRaster(rTetraAdv9, "models/tetraploidAdv1929.grd")

# develop testing and training sets for both cytotypes
fold <- kfold(both, k=5) #split occurence points into 5 sets
bothTest9 <- both[fold == 1, ] #take 20% (1/5) for testing
bothTrain9 <- both[fold != 1, ] #leave 40% for training
# fit training model for both cytotypes
maxBothTrain9 <- maxent(predictors9, BothTrain9) #fit maxent model
maxBothTrain9 #view results in html
rBothTrain9 <- predict(maxBothTrain9, predictors9) #predict full model
plot(rBothTrain9) #visualize full model
points(both) #add points to plot
# testing model for both cytotypes
# extract background points
bg9 <- randomPoints(predictors9, 1000)
# cross-validate model for both cytotypes
maxBothTest9 <- evaluate(maxBothTrain9, p=bothTest9, a=bg9, x=predictors9)
maxBothTest9 #print results
threshold(maxBothTest9) #identify threshold for presence or absence
plot(maxBothTest9, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for both cytotypes
pvtest9 <- data.frame(extract(predictors9, bothTest9))
avtest9 <- data.frame(extract(predictors9, bg9))
# cross-validate model
maxBothTest29 <- evaluate(maxBothTrain9, p=pvtest9, a=avtest9)
maxBothTest29
# Alternative 2: predict to testing points for both cytotypes
testp9 <- predict(maxBothTrain9, pvtest9)
testa9 <- predict(maxBothTrain9, avtest9)
maxBothTest39 <- evaluate(p=testp9, a=testa9)
maxBothTest39
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
points(tetraploid) # add points to predictive model
writeRaster(rBothAdv9, "models/bothAdv1929.grd")

### basic bioclim modeling with PRISM 2011 layers for diploids then tetraploids
# extract layer data for each point
dipPts11 <- extract(predictors11, diploid)
# create bioclim model
dipBC11 <- bioclim(dipPts11)
# predict bioclim model
dipBCpredict11 <- predict(predictors11, dipBC11)
# plot bioclim model
plot(dipBCpredict11)

# extract layer data for each point
tetraPts11 <- extract(predictors11, tetraploid)
# create bioclim model
tetraBC11 <- bioclim(tetraPts11)
# predict bioclim model
tetraBCpredict11 <- predict(predictors11, tetraBC11)
# plot bioclim model
plot(tetraBCpredict11)

# extract layer data for each point for both cytotypes
bothPts11 <- extract(predictors11, both)
# create bioclim model
bothBC11 <- bioclim(bothPts11)
# predict bioclim model
bothBCpredict11 <- predict(predictors11, bothBC11)
# plot bioclim model
plot(bothBCpredict11)

## Default maxent modeling
# run maxent for diploid (default parameters for dismo)
maxDip11 <- maxent(predictors11, diploid)
maxDip11 # views results in browser window
response(maxDip11) # show response curves for each layer
rDip11 <- predict(maxDip11, predictors11) # create model
plot(rDip11) # plot predictive model
points(diploid) # add points to predictive model
writeRaster(rDip11, "models/diploid2011.grd")

# run maxent for tetraploid (default parameters for dismo)
maxTetra11 <- maxent(predictors11, tetraploid) 
maxTetra11 
response(maxTetra11) 
rTetra11 <- predict(maxTetra11, predictors11) 
plot(rTetra11)
points(tetraploid)
writeRaster(rTetra11, "models/tetraploid2011.grd")

# run maxent for both cytotypes (default parameters for dismo)
maxBoth11 <- maxent(predictors11, both)
maxBoth11 # views results in browser window
response(maxBoth11) # show response curves for each layer
rBoth11 <- predict(maxBoth11, predictors11) # create model
plot(rBoth11) # plot predictive model
points(both) # add points to predictive model
writeRaster(rBoth11, "models/both2011.grd")

## Advanced modeling
# develop testing and training sets for diploid
fold <- kfold(diploid, k=5) #split occurence points into 5 sets
dipTest11 <- diploid[fold == 1, ] #take 20% (1/5) for testing
dipTrain11 <- diploid[fold != 1, ] #leave 40% for training
# fit training model for diploid
maxDipTrain11 <- maxent(predictors11, dipTrain11) #fit maxent model
maxDipTrain11 #view results in html
rDipTrain11 <- predict(maxDipTrain11, predictors11) #predict full model
plot(rDipTrain11) #visualize full model
points(diploid) #add points to plot
# testing model for diploid
# extract background points
bg11 <- randomPoints(predictors11, 1000)
# cross-validate model for diploid
maxDipTest11 <- evaluate(maxDipTrain11, p=dipTest11, a=bg11, x=predictors11)
maxDipTest11 #print results
threshold(maxDipTest11) #identify threshold for presence or absence
plot(maxDipTest11, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for diploid
pvtest11 <- data.frame(extract(predictors11, dipTest11))
avtest11 <- data.frame(extract(predictors11, bg11))
# cross-validate model
maxDipTest211 <- evaluate(maxDipTrain11, p=pvtest11, a=avtest11)
maxDipTest211
# Alternative 2: predict to testing points for diploid
testp11 <- predict(maxDipTrain11, pvtest11)
testa11 <- predict(maxDipTrain11, avtest11)
maxDipTest311 <- evaluate(p=testp11, a=testa11)
maxDipTest311
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
points(diploid) # add points to predictive model
writeRaster(rDipAdv11, "models/diploidAdv2011.grd")

# develop testing and training sets for tetraploid
fold <- kfold(tetraploid, k=5) #split occurence points into 5 sets
tetraTest11 <- tetraploid[fold == 1, ] #take 20% (1/5) for testing
tetraTrain11 <- tetraploid[fold != 1, ] #leave 40% for training
# fit training model for tetraploid
maxTetraTrain11 <- maxent(predictors11, tetraTrain11) #fit maxent model
maxTetraTrain11 #view results in html
rTetraTrain11 <- predict(maxTetraTrain11, predictors11) #predict full model
plot(rTetraTrain11) #visualize full model
points(tetraploid) #add points to plot
# testing model for tetraploid
# extract background points
bg11 <- randomPoints(predictors11, 1000)
# cross-validate model for tetraploid
maxTetraTest11 <- evaluate(maxTetraTrain11, p=tetraTest11, a=bg11, x=predictors11)
maxTetraTest11 #print results
threshold(maxTetraTest11) #identify threshold for presence or absence
plot(maxTetraTest11, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for tetraploid
pvtest11 <- data.frame(extract(predictors11, tetraTest11))
avtest11 <- data.frame(extract(predictors11, bg11))
# cross-validate model
maxTetraTest211 <- evaluate(maxTetraTrain11, p=pvtest11, a=avtest11)
maxTetraTest211
# Alternative 2: predict to testing points for tetraploid
testp11 <- predict(maxTetraTrain11, pvtest11)
testa11 <- predict(maxTetraTrain11, avtest11)
maxTetraTest311 <- evaluate(p=testp11, a=testa11)
maxTetraTest311
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
points(tetraploid) # add points to predictive model
writeRaster(rTetraAdv11, "models/tetraploidAdv2011.grd")

# develop testing and training sets for both cytotypes
fold <- kfold(both, k=5) #split occurence points into 5 sets
bothTest11 <- both[fold == 1, ] #take 20% (1/5) for testing
bothTrain11 <- both[fold != 1, ] #leave 40% for training
# fit training model for both cytotypes
maxBothTrain11 <- maxent(predictors11, BothTrain11) #fit maxent model
maxBothTrain11 #view results in html
rBothTrain11 <- predict(maxBothTrain11, predictors11) #predict full model
plot(rBothTrain11) #visualize full model
points(both) #add points to plot
# testing model for both cytotypes
# extract background points
bg11 <- randomPoints(predictors11, 1000)
# cross-validate model for both cytotypes
maxBothTest11 <- evaluate(maxBothTrain11, p=bothTest11, a=bg11, x=predictors11)
maxBothTest11 #print results
threshold(maxBothTest11) #identify threshold for presence or absence
plot(maxBothTest11, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for both cytotypes
pvtest11 <- data.frame(extract(predictors11, bothTest11))
avtest11 <- data.frame(extract(predictors11, bg11))
# cross-validate model
maxBothTest211 <- evaluate(maxBothTrain11, p=pvtest11, a=avtest11)
maxBothTest211
# Alternative 2: predict to testing points for both cytotypes
testp11 <- predict(maxBothTrain11, pvtest11)
testa11 <- predict(maxBothTrain11, avtest11)
maxBothTest311 <- evaluate(p=testp11, a=testa11)
maxBothTest311
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
points(tetraploid) # add points to predictive model
writeRaster(rBothAdv11, "models/bothAdv2011.grd")
