## evaluating niche overlap between cytotypes

library(dplyr)
library(raster)
library(dismo)
library(ecospat)
library(ENMeval)
library(fossil)
library(vegan)
library(rms)

## multi-variate climate space comparisons (standard statistical testing, non-model based)
# import occurrence data and convert to format required by maxent
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

## create stack of non-correlated layers (as determined by layerPrep.R)
predictors0<- stack(ppt0)
predictors1<- stack(ppt1,tmax1)

# extract layer data for each point and add label
dipPts0 <- raster::extract(predictors0, diploid)
dipPts0 <- cbind.data.frame(species="diploid", dipPts0) #add column for diploid
dipPts1 <- raster::extract(predictors1, diploid)
dipPts1 <- cbind.data.frame(species="diploid", dipPts1) #add column for diploid
tetraPts0 <- raster::extract(predictors0, tetraploid)
tetraPts0 <- cbind.data.frame(species="tetraploid", tetraPts0)
tetraPts0<-na.omit(tetraPts0)#removing NA values
tetraPts1 <- raster::extract(predictors1, tetraploid)
tetraPts1 <- cbind.data.frame(species="tetraploid", tetraPts1)
tetraPts1<-na.omit(tetraPts1)#removing NA values

# combine diploid and tetraploid
bothPts0 <- as.data.frame(rbind(dipPts0, tetraPts0))
bothPts1 <- as.data.frame(rbind(dipPts1, tetraPts1))

###Using PRISM 1930 weather data
##Instead of using an ANOVA, I will use a Logistic Regression
#independent continuous variables:weather layers
#dependent categoric variable:species/ploidy column

#import occurrence data for both Callisia cytotypes
both <- read.csv(file="both.csv")

#to make the dataframe with columns: cytotype labels, latitude,longitude, environmental variables containing extracted data 
#1930
master0<- cbind(both,bothPts0)
master0<- master0[,-4] 

#2014
master1<- cbind(both,bothPts1)
master1<- master1[,-4] 

#enter dataframe in logistic regression code
x0 <-master0
colnames(x0) <- tolower(colnames(x0))
x0[,2] <- as.numeric(x0[,2])
x0[,3] <- as.numeric(x0[,3]) 

#set up matrix to correct for spatial autocorrelation in logistic regression for 1930 
x.matrix0 <- earth.dist(x0[,c("latitude", "longitude")])
x.matrix0 <- as.matrix(x.matrix0)
pcnm.x.matrix0 <- pcnm(x.matrix0)
x.env.pcnm0 <- cbind(x0, pcnm.x.matrix0$vectors[,1:round((length(which(pcnm.x.matrix0$values > 0)))/2)])

for(i in 4:ncol(x0)){
  
  x.env.pcnm0 <- cbind(x.env.pcnm0, scale(x.env.pcnm0[,i], scale = sd(x.env.pcnm0[,i], na.rm = TRUE)))
  colnames(x.env.pcnm0)[ncol(x.env.pcnm0)] <- paste("std", colnames(x0)[i], sep = ".")
  
}

colnames(x.env.pcnm0)

#1930 - Look at the column names in the data frame "x.env.pcnm0". 
#Put all the columns named "PCNM..." into the model below, 
#followed by all the environmental variable columns starting with "std..."
model.lrm0 <- lrm(x.env.pcnm0[,"cytotype"] ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + std.avg_1930_ppt0, y = TRUE, data=x.env.pcnm0, penalty = 0.001)

#view results, variables with a p-value[Pr(>|z|)] of <0.05 are significant
#the PCNM variables are adjusting for spatial autocorrelation, 
#and allow me to make bolder statements about the variables that are influencing cytotype distribution
model.lrm0

###Using PRISM 2014 weather data
##Instead of using an ANOVA, I will use a Logistic Regression
#independent continuous variables:weather layers
#dependent categoric variable:species/ploidy column
x1 <-master1
colnames(x1) <- tolower(colnames(x1))
x1[,2] <- as.numeric(x1[,2])
x1[,3] <- as.numeric(x1[,3]) 

#set up matrix to correct for spatial autocorrelation in logistic regression for 2014 
x.matrix1 <- earth.dist(x1[,c("latitude", "longitude")])
x.matrix1 <- as.matrix(x.matrix1)
pcnm.x.matrix1 <- pcnm(x.matrix1)
x.env.pcnm1 <- cbind(x1, pcnm.x.matrix1$vectors[,1:round((length(which(pcnm.x.matrix1$values > 0)))/2)])

for(i in 4:ncol(x1)){
  
  x.env.pcnm1 <- cbind(x.env.pcnm1, scale(x.env.pcnm1[,i], scale = sd(x.env.pcnm1[,i], na.rm = TRUE)))
  colnames(x.env.pcnm1)[ncol(x.env.pcnm1)] <- paste("std", colnames(x1)[i], sep = ".")
  
}

colnames(x.env.pcnm1)

#2014 - Look at the column names in the data frame "x.env.pcnm1". 
#Put all the columns named "PCNM..." into the model below, 
#followed by all the environmental variable columns starting with "std..."
model.lrm1 <- lrm(x.env.pcnm1[,"cytotype"] ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + std.avg_2014_ppt0 + std.avg_2014_tmax0, y = TRUE, data=x.env.pcnm1, penalty = 0.001)

#view results, variables with a p-value[Pr(>|z|)] of <0.05 are significant
#the PCNM variables are adjusting for spatial autocorrelation, 
#and allow me to make bolder statements about the variables that are influencing cytotype distribution
model.lrm1

## principle component analysis(PCA)
bothNum0 <- bothPts0[ ,-1] #remove species names
pca_both0 <- prcomp(bothNum0, center = TRUE, scale. = TRUE) #PCA=Error because only has one weather variable
print(pca_both0) #print deviations and rotations=Error because only has one weather variable
summary(pca_both0) #print importance of components=Error because only has one weather variable
plot(pca_both0, type="l") #plot variances=Error because only has one weather variable
ncomp <- 1#specify number of components to load (representing 99% of variation)=Error because only has one weather variable

## model-based approaches
# read in default maxent models
rDip0 <- raster("models/diploid1930.grd")
rTetra0 <- raster("models/tetraploid1930.grd")
rBoth0<- raster("models/both1930.grd")

# assessing niche overlap by comparing diploids and tetraploids in 1930
nicheOverlap(rDip0, rTetra0, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip0, rTetra0, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

# read in advanced maxent models
rDipAdv0 <- raster("models/diploidAdv1930.grd")
rTetraAdv0 <- raster("models/tetraploidAdv1930.grd")
rBothAdv0 <- raster("models/bothAdv1930.grd")

# assessing niche overlap by comparing diploids and tetraploids in 1930
nicheOverlap(rDipAdv0, rTetraAdv0, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDipAdv0, rTetraAdv0, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

# principle component analysis(PCA)
bothNum1 <- bothPts1[ ,-1] #remove species names
pca_both1 <- prcomp(bothNum1, center = TRUE, scale. = TRUE) #PCA
print(pca_both1) #print deviations and rotations
summary(pca_both1) #print importance of components
plot(pca_both1, type="l") #plot variances
ncomp <- 2 #specify number of components to load (representing 99% of variation)

### model-based approaches
## read in default maxent models
rDip1 <- raster("models/diploid2014.grd")
rTetra1 <- raster("models/tetraploid2014.grd")
rBoth1<- raster("models/both2014.grd")

# assessing niche overlap by comparing diploids and tetraploids in 2014
nicheOverlap(rDip1, rTetra1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip1, rTetra1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in diploid niche from 1930 to 2014
nicheOverlap(rDip0, rDip1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip0, rDip1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in tetraploid niche from 1930 to 2014
nicheOverlap(rTetra0, rTetra1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rTetra0, rTetra1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in both cytotypes niche from 1930 to 2014
nicheOverlap(rBoth0, rBoth1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rBoth0, rBoth1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

## read in advanced maxent models
rDipAdv1 <- raster("models/diploidAdv2014.grd")
rTetraAdv1 <- raster("models/tetraploidAdv2014.grd")
rBothAdv1 <- raster("models/bothAdv2014.grd")

# assessing niche overlap by comparing diploids and tetraploids in 2014
nicheOverlap(rDipAdv1, rTetraAdv1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDipAdv1, rTetraAdv1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in diploid niche from 1930 to 2014
nicheOverlap(rDipAdv0, rDipAdv1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDipAdv0, rDipAdv1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in tetraploid niche from 1930 to 2014
nicheOverlap(rTetraAdv0, rTetraAdv1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rTetraAdv0, rTetraAdv1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in both cytotypes niche from 1930 to 2014
nicheOverlap(rBothAdv0, rBothAdv1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rBothAdv0, rBothAdv1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

##Permutation test which determines if the two cytotpes distribution are truly different not just different by chance
#copy and paste the _sorted.csv file in the rresult folder within ENMTools to the working directory
# load results of permuted model replications csv
y <- read.csv("IDENTITY_tetraploid_vs_diploid_sorted.csv")

#This data is permuted, meaning that the data points of the two species were
#scambled for each run, so that there really is no niche difference among 
#the two "species." Then for each of those runs, the I and D stats were 
#calculated. This gives you an idea of the range of I and D scores you would
#expect to see when in fact the underlying points DO NOT differ in their 
#niches. If you want to claim that your two species differ in their niches,
#the I or D score you get should be LOWER than the I or D 
#values that were calculated using the fake data (lower than 95% of the I or
#D scores from the fake data or more).
quantile(y[,1], 0.05)

#This gives you the 5% threshold of permuted I values
#Change this to y[,2] if you want to get this stat for D instead of I.
quantile(y[,2], 0.05)

#Now compare this 5% threshold to the observed I point estimate (you should
#have calculated this point estimate earlier). If the point estimate is 
#lower than the 5% threshold, then the niches are significantly different.
=======
  ## evaluating niche overlap between taxa
  
  library(dplyr)
library(raster)
library(dismo)
library(ecospat)
library(ENMeval)
library(fossil)
library(vegan)
library(rms)

## multi-variate climate space comparisons (standard statistical testing, non-model based)
# import occurrence data and convert to format required by maxent
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

## create stack of non-correlated layers (as determined by layerPrep.R)
predictors0<- stack(ppt0)
predictors1<- stack(ppt1,tmax1)

# extract layer data for each point and add label
dipPts0 <- raster::extract(predictors0, diploid)
dipPts0 <- cbind.data.frame(species="diploid", dipPts0) #add column for diploid
dipPts1 <- raster::extract(predictors1, diploid)
dipPts1 <- cbind.data.frame(species="diploid", dipPts1) #add column for diploid
tetraPts0 <- raster::extract(predictors0, tetraploid)
tetraPts0 <- cbind.data.frame(species="tetraploid", tetraPts0)
tetraPts0<-na.omit(tetraPts0)#removing NA values
tetraPts1 <- raster::extract(predictors1, tetraploid)
tetraPts1 <- cbind.data.frame(species="tetraploid", tetraPts1)
tetraPts1<-na.omit(tetraPts1)#removing NA values

# combine diploid and tetraploid
bothPts0 <- as.data.frame(rbind(dipPts0, tetraPts0))
bothPts1 <- as.data.frame(rbind(dipPts1, tetraPts1))

###Using PRISM 1930 weather data
##Instead of using an ANOVA, I will use a Logistic Regression
#independent continuous variables:weather layers
#dependent categoric variable:species/ploidy column

#import occurrence data for both Callisia cytotypes
both <- read.csv(file="both.csv")

#to make the dataframe with columns: cytotype labels, latitude,longitude, environmental variables containing extracted data 
#1930
master0<- cbind(both,bothPts0)
master0<- master0[,-4] 

#2014
master1<- cbind(both,bothPts1)
master1<- master1[,-4] 

#enter dataframe in logistic regression code
x0 <-master0
colnames(x0) <- tolower(colnames(x0))
x0[,2] <- as.numeric(x0[,2])
x0[,3] <- as.numeric(x0[,3]) 

#set up matrix to correct for spatial autocorrelation in logistic regression for 1930 
x.matrix0 <- earth.dist(x0[,c("latitude", "longitude")])
x.matrix0 <- as.matrix(x.matrix0)
pcnm.x.matrix0 <- pcnm(x.matrix0)
x.env.pcnm0 <- cbind(x0, pcnm.x.matrix0$vectors[,1:round((length(which(pcnm.x.matrix0$values > 0)))/2)])

for(i in 4:ncol(x0)){
  
  x.env.pcnm0 <- cbind(x.env.pcnm0, scale(x.env.pcnm0[,i], scale = sd(x.env.pcnm0[,i], na.rm = TRUE)))
  colnames(x.env.pcnm0)[ncol(x.env.pcnm0)] <- paste("std", colnames(x0)[i], sep = ".")
  
}

colnames(x.env.pcnm0)

#1930 - Look at the column names in the data frame "x.env.pcnm0". 
#Put all the columns named "PCNM..." into the model below, 
#followed by all the environmental variable columns starting with "std..."
model.lrm0 <- lrm(x.env.pcnm0[,"cytotype"] ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + std.avg_1930_ppt0, y = TRUE, data=x.env.pcnm0, penalty = 0.001)

#view results, variables with a p-value[Pr(>|z|)] of <0.05 are significant
#the PCNM variables are adjusting for spatial autocorrelation, 
#and allow me to make bolder statements about the variables that are influencing cytotype distribution
model.lrm0

###Using PRISM 2014 weather data
##Instead of using an ANOVA, I will use a Logistic Regression
#independent continuous variables:weather layers
#dependent categoric variable:species/ploidy column
x1 <-master1
colnames(x1) <- tolower(colnames(x1))
x1[,2] <- as.numeric(x1[,2])
x1[,3] <- as.numeric(x1[,3]) 

#set up matrix to correct for spatial autocorrelation in logistic regression for 2014 
x.matrix1 <- earth.dist(x1[,c("latitude", "longitude")])
x.matrix1 <- as.matrix(x.matrix1)
pcnm.x.matrix1 <- pcnm(x.matrix1)
x.env.pcnm1 <- cbind(x1, pcnm.x.matrix1$vectors[,1:round((length(which(pcnm.x.matrix1$values > 0)))/2)])

for(i in 4:ncol(x1)){
  
  x.env.pcnm1 <- cbind(x.env.pcnm1, scale(x.env.pcnm1[,i], scale = sd(x.env.pcnm1[,i], na.rm = TRUE)))
  colnames(x.env.pcnm1)[ncol(x.env.pcnm1)] <- paste("std", colnames(x1)[i], sep = ".")
  
}

colnames(x.env.pcnm1)

#2014 - Look at the column names in the data frame "x.env.pcnm1". 
#Put all the columns named "PCNM..." into the model below, 
#followed by all the environmental variable columns starting with "std..."
model.lrm1 <- lrm(x.env.pcnm1[,"cytotype"] ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + std.avg_2014_ppt0 + std.avg_2014_tmax0, y = TRUE, data=x.env.pcnm1, penalty = 0.001)

#view results, variables with a p-value[Pr(>|z|)] of <0.05 are significant
#the PCNM variables are adjusting for spatial autocorrelation, 
#and allow me to make bolder statements about the variables that are influencing cytotype distribution
model.lrm1

## principle component analysis(PCA)
bothNum0 <- bothPts0[ ,-1] #remove species names
pca_both0 <- prcomp(bothNum0, center = TRUE, scale. = TRUE) #PCA=Error because only has one weather variable
print(pca_both0) #print deviations and rotations=Error because only has one weather variable
summary(pca_both0) #print importance of components=Error because only has one weather variable
plot(pca_both0, type="l") #plot variances=Error because only has one weather variable
ncomp <- 1#specify number of components to load (representing 99% of variation)=Error because only has one weather variable

## model-based approaches
# read in default maxent models
rDip0 <- raster("models/diploid1930.grd")
rTetra0 <- raster("models/tetraploid1930.grd")
rBoth0<- raster("models/both1930.grd")

# assessing niche overlap by comparing diploids and tetraploids in 1930
nicheOverlap(rDip0, rTetra0, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip0, rTetra0, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

# read in advanced maxent models
rDipAdv0 <- raster("models/diploidAdv1930.grd")
rTetraAdv0 <- raster("models/tetraploidAdv1930.grd")
rBothAdv0 <- raster("models/bothAdv1930.grd")

# assessing niche overlap by comparing diploids and tetraploids in 1930
nicheOverlap(rDipAdv0, rTetraAdv0, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDipAdv0, rTetraAdv0, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

# principle component analysis(PCA)
bothNum1 <- bothPts1[ ,-1] #remove species names
pca_both1 <- prcomp(bothNum1, center = TRUE, scale. = TRUE) #PCA
print(pca_both1) #print deviations and rotations
summary(pca_both1) #print importance of components
plot(pca_both1, type="l") #plot variances
ncomp <- 2 #specify number of components to load (representing 99% of variation)

### model-based approaches
## read in default maxent models
rDip1 <- raster("models/diploid2014.grd")
rTetra1 <- raster("models/tetraploid2014.grd")
rBoth1<- raster("models/both2014.grd")

# assessing niche overlap by comparing diploids and tetraploids in 2014
nicheOverlap(rDip1, rTetra1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip1, rTetra1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in diploid niche from 1930 to 2014
nicheOverlap(rDip0, rDip1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip0, rDip1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in tetraploid niche from 1930 to 2014
nicheOverlap(rTetra0, rTetra1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rTetra0, rTetra1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in both cytotypes niche from 1930 to 2014
nicheOverlap(rBoth0, rBoth1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rBoth0, rBoth1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

## read in advanced maxent models
rDipAdv1 <- raster("models/diploidAdv2014.grd")
rTetraAdv1 <- raster("models/tetraploidAdv2014.grd")
rBothAdv1 <- raster("models/bothAdv2014.grd")

# assessing niche overlap by comparing diploids and tetraploids in 2014
nicheOverlap(rDipAdv1, rTetraAdv1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDipAdv1, rTetraAdv1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in diploid niche from 1930 to 2014
nicheOverlap(rDipAdv0, rDipAdv1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDipAdv0, rDipAdv1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in tetraploid niche from 1930 to 2014
nicheOverlap(rTetraAdv0, rTetraAdv1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rTetraAdv0, rTetraAdv1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in both cytotypes niche from 1930 to 2014
nicheOverlap(rBothAdv0, rBothAdv1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rBothAdv0, rBothAdv1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

##Permutation test which determines if the two cytotpes distribution are truly different not just different by chance
#copy and paste the _sorted.csv file in the rresult folder within ENMTools to the working directory
# load results of permuted model replications csv
y <- read.csv("IDENTITY_tetraploid_vs_diploid_sorted.csv")

#This data is permuted, meaning that the data points of the two species were
#scambled for each run, so that there really is no niche difference among 
#the two "species." Then for each of those runs, the I and D stats were 
#calculated. This gives you an idea of the range of I and D scores you would
#expect to see when in fact the underlying points DO NOT differ in their 
#niches. If you want to claim that your two species differ in their niches,
#the I or D score you get should be LOWER than the I or D 
#values that were calculated using the fake data (lower than 95% of the I or
#D scores from the fake data or more).
quantile(y[,1], 0.05)

#This gives you the 5% threshold of permuted I values
#Change this to y[,2] if you want to get this stat for D instead of I.
quantile(y[,2], 0.05)

#Now compare this 5% threshold to the observed I point estimate (you should
#have calculated this point estimate earlier). If the point estimate is 
#lower than the 5% threshold, then the niches are significantly different.