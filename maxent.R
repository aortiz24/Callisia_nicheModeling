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
ppt0 <- raster("layers/avg_1930_ppt0.asc", crs=CRS)
tmax0 <- raster("layers/avg_1930_tmax0.asc", crs=CRS)
tmean0 <- raster("layers/avg_1930_tmean0.asc", crs=CRS)
tmin0 <- raster("layers/avg_1930_tmin0.asc", crs=CRS)
ppt1 <- raster("layers/avg_2014_ppt0.asc", crs=CRS)
tmax1 <- raster("layers/avg_2014_tmax0.asc", crs=CRS)
tmean1 <- raster("layers/avg_2014_tmean0.asc", crs=CRS)
tmin1 <- raster("layers/avg_2014_tmin0.asc", crs=CRS)

## create stack of non-correlated layers (as determined by layerPrep.R)
predictors <- stack(alt, bio2, bio3, bio5, bio6, bio8, bio9, bio12, bio13, bio14, bio19)
predictors0<- stack(ppt0)
predictors1<- stack(ppt1,tmax1)

# plot each layer individually
plot(predictors)
plot(predictors0)
plot(predictors1)

### basic bioclim modeling with Bioclim layers for diploids then tetraploids
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

## Advanced modeling
# develop testing and training sets for diploid
fold0 <- kfold(diploid, k=5) #split occurence points into 5 sets
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
# maxent with jackknife, random seed, and response curves, followed by cross-validation
maxDipAdv0 <- maxent(
  x=predictors0,
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
maxDipAdv0 #view output as html

# develop testing and training sets for tetraploid
fold0 <- kfold(tetraploid, k=5) #split occurence points into 5 sets
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
maxTetraAdv0 <- maxent(
  x=predictors0,
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
maxTetraAdv0 #view output as html

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

## Advanced modeling
# develop testing and training sets for diploid
fold1 <- kfold(diploid, k=5) #split occurence points into 5 sets
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

# develop testing and training sets for tetraploid
fold1 <- kfold(tetraploid, k=5) #split occurence points into 5 sets
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