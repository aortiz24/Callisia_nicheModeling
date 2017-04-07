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

#layers ending in 0 are for PRISM1930
#layers ending in 1 are for PRISM2014
# import layers with CRS specified
CRS <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"
ppt0 <- raster("layers/avg_1930_ppt0.asc", crs=CRS)
tmax0 <- raster("layers/avg_1930_tmax0.asc", crs=CRS)
tmean0 <- raster("layers/avg_1930_tmean0.asc", crs=CRS)
tmin0 <- raster("layers/avg_1930_tmin0.asc", crs=CRS)
ppt1 <- raster("layers/avg_2014_ppt0.asc", crs=CRS)
tmax1 <- raster("layers/avg_2014_tmax0.asc", crs=CRS)
tmean1 <- raster("layers/avg_2014_tmean0.asc", crs=CRS)
tmin1 <- raster("layers/avg_2014_tmin0.asc", crs=CRS)
ppt9 <- raster("layers/ppt9.asc")
tmax9 <- raster("layers/tmax9.asc")
tmean9 <- raster("layers/tmean9.asc")
tmin9 <- raster("layers/tmin9.asc")
vpdmax9 <- raster("layers/vpdmax9.asc")
vpdmin9 <- raster("layers/vpdmin9.asc")
tdmean9 <- raster("layers/tdmean9.asc")
ppt11 <- raster("layers/ppt11.asc")
tmax11 <- raster("layers/tmax11.asc")
tmean11 <- raster("layers/tmean11.asc")
tmin11 <- raster("layers/tmin11.asc")
vpdmax11 <- raster("layers/vpdmax11.asc")
vpdmin11 <- raster("layers/vpdmin11.asc")
tdmean11 <- raster("layers/tdmean11.asc")

## create stack of non-correlated layers (as determined by layerPrep.R)
predictors0<- stack(ppt0)
predictors1<- stack(ppt1,tmax1)
predictors9<- stack(tmean9, vpdmin9)
predictors11<- stack(tmin11, ppt11, vpdmin11)

# plot each layer individually
plot(predictors0)
plot(predictors1)
plot(predictors9)
plot(predictors11)

### basic bioclim modeling with PRISM 1930 layers for diploids then tetraploids
# extract layer data for each point
dipPts0 <- extract(predictors0, diploid)
# create bioclim model
dipBC0 <- bioclim(dipPts0)
# predict bioclim model
dipBCpredict0 <- predict(predictors0, dipBC0)
# plot bioclim model
plot(dipBCpredict0)

# extract layer data for each point
tetraPts0 <- extract(predictors0, tetraploid)
# create bioclim model
tetraBC0 <- bioclim(tetraPts0)
# predict bioclim model
tetraBCpredict0 <- predict(predictors0, tetraBC0)
# plot bioclim model
plot(tetraBCpredict0)

# extract layer data for each point for both cytotypes
bothPts0 <- extract(predictors0, both)
# create bioclim model
bothBC0 <- bioclim(bothPts0)
# predict bioclim model
bothBCpredict0 <- predict(predictors0, bothBC0)
# plot bioclim model
plot(bothBCpredict0)

## Default maxent modeling
# run maxent for diploid (default parameters for dismo)
maxDip0 <- maxent(predictors0, diploid)
maxDip0 # views results in browser window
response(maxDip0) # show response curves for each layer
rDip0 <- predict(maxDip0, predictors0) # create model
plot(rDip0) # plot predictive model
points(diploid) # add points to predictive model
writeRaster(rDip0, "models/diploid1930.grd")

# run maxent for tetraploid (default parameters for dismo)
maxTetra0 <- maxent(predictors0, tetraploid) 
maxTetra0 
response(maxTetra0) 
rTetra0 <- predict(maxTetra0, predictors0) 
plot(rTetra0)
points(tetraploid)
writeRaster(rTetra0, "models/tetraploid1930.grd")

# run maxent for both cytotypes (default parameters for dismo)
maxBoth0 <- maxent(predictors0, both)
maxBoth0 # views results in browser window
response(maxBoth0) # show response curves for each layer
rBoth0 <- predict(maxBoth0, predictors0) # create model
plot(rBoth0) # plot predictive model
points(both) # add points to predictive model
writeRaster(rBoth0, "models/both1930.grd")

## Advanced modeling
# develop testing and training sets for diploid
fold <- kfold(diploid, k=5) #split occurence points into 5 sets
dipTest0 <- diploid[fold == 1, ] #take 20% (1/5) for testing
dipTrain0 <- diploid[fold != 1, ] #leave 40% for training
# fit training model for diploid
maxDipTrain0 <- maxent(predictors0, dipTrain0) #fit maxent model
maxDipTrain0 #view results in html
rDipTrain0 <- predict(maxDipTrain0, predictors0) #predict full model
plot(rDipTrain0) #visualize full model
points(diploid) #add points to plot
# testing model for diploid
# extract background points
bg0 <- randomPoints(predictors0, 1000)
# cross-validate model for diploid
maxDipTest0 <- evaluate(maxDipTrain0, p=dipTest0, a=bg0, x=predictors0)
maxDipTest0 #print results
threshold(maxDipTest0) #identify threshold for presence or absence
plot(maxDipTest0, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for diploid
pvtest0 <- data.frame(extract(predictors0, dipTest0))
avtest0 <- data.frame(extract(predictors0, bg0))
# cross-validate model
maxDipTest20 <- evaluate(maxDipTrain0, p=pvtest0, a=avtest0)
maxDipTest20
# Alternative 2: predict to testing points for diploid
testp0 <- predict(maxDipTrain0, pvtest0)
testa0 <- predict(maxDipTrain0, avtest0)
maxDipTest30 <- evaluate(p=testp0, a=testa0)
maxDipTest30
# maxent with jackknife, random seed, and response curves, followed by cross
#jackknife not run because only one layer in predictor0
maxDipAdv0 <- maxent(
  x=predictors0,
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
maxDipAdv0 #view output as html
response(maxDipAdv0) # show response curves for each layer
rDipAdv0 <- predict(maxDipAdv0, predictors0) # create model
plot(rDipAdv0) # plot predictive model
points(diploid) # add points to predictive model
writeRaster(rDipAdv0, "models/diploidAdv1930.grd")

# develop testing and training sets for tetraploid
fold <- kfold(tetraploid, k=5) #split occurence points into 5 sets
tetraTest0 <- tetraploid[fold == 1, ] #take 20% (1/5) for testing
tetraTrain0 <- tetraploid[fold != 1, ] #leave 40% for training
# fit training model for tetraploid
maxTetraTrain0 <- maxent(predictors0, tetraTrain0) #fit maxent model
maxTetraTrain0 #view results in html
rTetraTrain0 <- predict(maxTetraTrain0, predictors0) #predict full model
plot(rTetraTrain0) #visualize full model
points(tetraploid) #add points to plot
# testing model for tetraploid
# extract background points
bg0 <- randomPoints(predictors0, 1000)
# cross-validate model for tetraploid
maxTetraTest0 <- evaluate(maxTetraTrain0, p=tetraTest0, a=bg0, x=predictors0)
maxTetraTest0 #print results
threshold(maxTetraTest0) #identify threshold for presence or absence
plot(maxTetraTest0, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for tetraploid
pvtest0 <- data.frame(extract(predictors0, tetraTest0))
avtest0 <- data.frame(extract(predictors0, bg0))
# cross-validate model
maxTetraTest20 <- evaluate(maxTetraTrain0, p=pvtest0, a=avtest0)
maxTetraTest20
# Alternative 2: predict to testing points for tetraploid
testp0 <- predict(maxTetraTrain0, pvtest0)
testa0 <- predict(maxTetraTrain0, avtest0)
maxTetraTest30 <- evaluate(p=testp0, a=testa0)
maxTetraTest30
# maxent with jackknife, random seed, and response curves, followed by cross-validation
#jackknife not run because only one layer in predictor0
maxTetraAdv0 <- maxent(
  x=predictors0,
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
maxTetraAdv0 #view output as html
response(maxTetraAdv0) # show response curves for each layer
rTetraAdv0 <- predict(maxTetraAdv0, predictors0) # create model
plot(rTetraAdv0) # plot predictive model
points(tetraploid) # add points to predictive model
writeRaster(rTetraAdv0, "models/tetraploidAdv1930.grd")

# develop testing and training sets for both cytotypes
fold <- kfold(both, k=5) #split occurence points into 5 sets
bothTest0 <- both[fold == 1, ] #take 20% (1/5) for testing
bothTrain0 <- both[fold != 1, ] #leave 40% for training
# fit training model for both cytotypes
maxBothTrain0 <- maxent(predictors0, BothTrain0) #fit maxent model
maxBothTrain0 #view results in html
rBothTrain0 <- predict(maxBothTrain0, predictors0) #predict full model
plot(rBothTrain0) #visualize full model
points(both) #add points to plot
# testing model for both cytotypes
# extract background points
bg0 <- randomPoints(predictors0, 1000)
# cross-validate model for both cytotypes
maxBothTest0 <- evaluate(maxBothTrain0, p=bothTest0, a=bg0, x=predictors0)
maxBothTest0 #print results
threshold(maxBothTest0) #identify threshold for presence or absence
plot(maxBothTest0, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for both cytotypes
pvtest0 <- data.frame(extract(predictors0, bothTest0))
avtest0 <- data.frame(extract(predictors0, bg0))
# cross-validate model
maxBothTest20 <- evaluate(maxBothTrain0, p=pvtest0, a=avtest0)
maxBothTest20
# Alternative 2: predict to testing points for both cytotypes
testp0 <- predict(maxBothTrain0, pvtest0)
testa0 <- predict(maxBothTrain0, avtest0)
maxBothTest30 <- evaluate(p=testp0, a=testa0)
maxBothTest30
# maxent with jackknife, random seed, and response curves, followed by cross
#jackknife not run because only one layer in predictor0
maxBothAdv0 <- maxent(
  x=predictors0,
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
maxBothAdv0 #view output as html
response(maxBothAdv0) # show response curves for each layer
rBothAdv0 <- predict(maxBothAdv0, predictors0) # create model
plot(rBothAdv0) # plot predictive model
points(tetraploid) # add points to predictive model
writeRaster(rBothAdv0, "models/bothAdv1930.grd")

### basic bioclim modeling with PRISM 2014 layers for diploids then tetraploids
# extract layer data for each point
dipPts1 <- extract(predictors1, diploid)
# create bioclim model
dipBC1 <- bioclim(dipPts1)
# predict bioclim model
dipBCpredict1 <- predict(predictors1, dipBC1)
# plot bioclim model
plot(dipBCpredict1)

# extract layer data for each point
tetraPts1 <- extract(predictors1, tetraploid)
# create bioclim model
tetraBC1 <- bioclim(tetraPts1)
# predict bioclim model
tetraBCpredict1 <- predict(predictors1, tetraBC1)
# plot bioclim model
plot(tetraBCpredict1)

# extract layer data for each point of both cytotypes
bothPts1 <- extract(predictors1, both)
# create bioclim model
bothBC1 <- bioclim(bothPts1)
# predict bioclim model
bothBCpredict1 <- predict(predictors1, bothBC1)
# plot bioclim model
plot(bothBCpredict1)

## Default maxent modeling
# run maxent for diploid (default parameters for dismo)
maxDip1 <- maxent(predictors1, diploid)
maxDip1 # views results in browser window
response(maxDip1) # show response curves for each layer
rDip1 <- predict(maxDip1, predictors1) # create model
plot(rDip1) # plot predictive model
points(diploid) # add points to predictive model
writeRaster(rDip1, "models/diploid2014.grd")

# run maxent for tetraploid (default parameters for dismo)
maxTetra1 <- maxent(predictors1, tetraploid) 
maxTetra1 
response(maxTetra1) 
rTetra1 <- predict(maxTetra1, predictors1) 
plot(rTetra1)
points(tetraploid)
writeRaster(rTetra1, "models/tetraploid2014.grd")

# run maxent for both cytotypes (default parameters for dismo)
maxBoth1 <- maxent(predictors1, both)
maxBoth1 # views results in browser window
response(maxBoth1) # show response curves for each layer
rBoth1 <- predict(maxBoth1, predictors1) # create model
plot(rBoth1) # plot predictive model
points(both) # add points to predictive model
writeRaster(rBoth1, "models/both2014.grd")

## Advanced modeling
# develop testing and training sets for diploid
fold <- kfold(diploid, k=5) #split occurence points into 5 sets
dipTest1 <- diploid[fold == 1, ] #take 20% (1/5) for testing
dipTrain1 <- diploid[fold != 1, ] #leave 40% for training
# fit training model for diploid
maxDipTrain1 <- maxent(predictors1, dipTrain1) #fit maxent model
maxDipTrain1 #view results in html
rDipTrain1 <- predict(maxDipTrain1, predictors1) #predict full model
plot(rDipTrain1) #visualize full model
points(diploid) #add points to plot
# testing model for diploid
# extract background points
bg1 <- randomPoints(predictors1, 1000)
# cross-validate model for diploid
maxDipTest1 <- evaluate(maxDipTrain1, p=dipTest1, a=bg1, x=predictors1)
maxDipTest1 #print results
threshold(maxDipTest1) #identify threshold for presence or absence
plot(maxDipTest1, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for diploid
pvtest1 <- data.frame(extract(predictors1, dipTest1))
avtest1 <- data.frame(extract(predictors1, bg1))
# cross-validate model
maxDipTest21 <- evaluate(maxDipTrain1, p=pvtest1, a=avtest1)
maxDipTest21
# Alternative 2: predict to testing points for diploid
testp1 <- predict(maxDipTrain1, pvtest1)
testa1 <- predict(maxDipTrain1, avtest1)
maxDipTest31 <- evaluate(p=testp1, a=testa1)
maxDipTest31
# maxent with jackknife, random seed, and response curves, followed by cross-validation
maxDipAdv1 <- maxent(
  x=predictors1,
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
maxDipAdv1 #view output as html
response(maxDipAdv1) # show response curves for each layer
rDipAdv1 <- predict(maxDipAdv1, predictors1) # create model
plot(rDipAdv1) # plot predictive model
points(diploid) # add points to predictive model
writeRaster(rDipAdv1, "models/diploidAdv2014.grd")

# develop testing and training sets for tetraploid
fold <- kfold(tetraploid, k=5) #split occurence points into 5 sets
tetraTest1 <- tetraploid[fold == 1, ] #take 20% (1/5) for testing
tetraTrain1 <- tetraploid[fold != 1, ] #leave 40% for training
# fit training model for tetraploid
maxTetraTrain1 <- maxent(predictors1, tetraTrain1) #fit maxent model
maxTetraTrain1 #view results in html
rTetraTrain1 <- predict(maxTetraTrain1, predictors1) #predict full model
plot(rTetraTrain1) #visualize full model
points(tetraploid) #add points to plot
# testing model for tetraploid
# extract background points
bg1 <- randomPoints(predictors1, 1000)
# cross-validate model for tetraploid
maxTetraTest1 <- evaluate(maxTetraTrain1, p=tetraTest1, a=bg1, x=predictors1)
maxTetraTest1 #print results
threshold(maxTetraTest1) #identify threshold for presence or absence
plot(maxTetraTest1, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for tetraploid
pvtest1 <- data.frame(extract(predictors1, tetraTest1))
avtest1 <- data.frame(extract(predictors1, bg1))
# cross-validate model
maxTetraTest21 <- evaluate(maxTetraTrain1, p=pvtest1, a=avtest1)
maxTetraTest21
# Alternative 2: predict to testing points for tetraploid
testp1 <- predict(maxTetraTrain1, pvtest1)
testa1 <- predict(maxTetraTrain1, avtest1)
maxTetraTest31 <- evaluate(p=testp1, a=testa1)
maxTetraTest31
# maxent with jackknife, random seed, and response curves, followed by cross-validation
maxTetraAdv1 <- maxent(
  x=predictors1,
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
maxTetraAdv1 #view output as html
response(maxTetraAdv1) # show response curves for each layer
rTetraAdv1 <- predict(maxTetraAdv1, predictors1) # create model
plot(rTetraAdv1) # plot predictive model
points(tetraploid) # add points to predictive model
writeRaster(rTetraAdv1, "models/tetraploidAdv2014.grd")

# develop testing and training sets for both cytotypes
fold <- kfold(both, k=5) #split occurence points into 5 sets
bothTest1 <- both[fold == 1, ] #take 20% (1/5) for testing
bothTrain1 <- both[fold != 1, ] #leave 40% for training
# fit training model for both cytotypes
maxBothTrain1 <- maxent(predictors1, bothTrain1) #fit maxent model
maxBothTrain1 #view results in html
rBothTrain1 <- predict(maxBothTrain1, predictors1) #predict full model
plot(rBothTrain1) #visualize full model
points(both) #add points to plot
# testing model for both cytotypes
# extract background points
bg1 <- randomPoints(predictors1, 1000)
# cross-validate model for both cytotypes
maxBothTest1 <- evaluate(maxBothTrain1, p=bothTest1, a=bg1, x=predictors1)
maxBothTest1 #print results
threshold(maxBothTest1) #identify threshold for presence or absence
plot(maxBothTest1, 'ROC') #plot AUC
# alternative methods for testing models (should give same answers)
# Alternative 1: another way to test model for both cytotypes
pvtest1 <- data.frame(extract(predictors1, bothTest1))
avtest1 <- data.frame(extract(predictors1, bg1))
# cross-validate model
maxBothTest21 <- evaluate(maxBothTrain1, p=pvtest1, a=avtest1)
maxBothTest21
# Alternative 2: predict to testing points for both cytotypes
testp1 <- predict(maxBothTrain1, pvtest1)
testa1 <- predict(maxBothTrain1, avtest1)
maxBothTest31 <- evaluate(p=testp1, a=testa1)
maxBothTest31
# maxent with jackknife, random seed, and response curves, followed by cross-validation
maxBothAdv1 <- maxent(
  x=predictors1,
  p=both,
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
maxBothAdv1 #view output as html
response(maxBothAdv1) # show response curves for each layer
rBothAdv1 <- predict(maxBothAdv1, predictors1) # create model
plot(rBothAdv1) # plot predictive model
points(tetraploid) # add points to predictive model
writeRaster(rBothAdv1, "models/bothAdv2014.grd")
