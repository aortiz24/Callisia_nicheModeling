## evaluating niche overlap between cytotypes
#load libraries
library(dplyr)
library(raster)
library(dismo)
library(ecospat)
library(ENMeval)
library(fossil)
library(vegan)
library(rms)
library(devtools)
library(ggbiplot)

## multi-variate climate space comparisons (standard statistical testing, non-model based)
#import occurrence data and convert to format required by maxent
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

#layers ending in 9 are for PRISM11929
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
predictors9<- stack(tmean9, ppt9, vpdmax9)
predictors11<- stack(tmean11, ppt11, vpdmin11, tdmean11)

# extract layer data for each point and add label
dipPts9 <- raster::extract(predictors9, diploid)
dipPts9 <- cbind.data.frame(species="diploid", dipPts9) #add column for diploid
dipPts11 <- raster::extract(predictors11, diploid)
dipPts11 <- cbind.data.frame(species="diploid", dipPts11) #add column for diploid
tetraPts9 <- raster::extract(predictors9, tetraploid)
tetraPts9 <- cbind.data.frame(species="tetraploid", tetraPts9)
tetraPts9<-na.omit(tetraPts9)#removing NA values
tetraPts11 <- raster::extract(predictors11, tetraploid)
tetraPts11 <- cbind.data.frame(species="tetraploid", tetraPts11)
tetraPts11<-na.omit(tetraPts11)#removing NA values

# combine diploid and tetraploid
bothPts9 <- as.data.frame(rbind(dipPts9, tetraPts9))
bothPts11 <- as.data.frame(rbind(dipPts11, tetraPts11))

#import occurrence data for both Callisia cytotypes
both <- read.csv(file="both.csv")

###Using PRISM 1929 weather data
##Instead of using an ANOVA, I will use a Logistic Regression
#independent continuous variables:weather layers
#dependent categoric variable:species/ploidy column

#to make the dataframe with columns: cytotype labels, latitude,longitude, environmental variables containing extracted data 
#1929
master9<- cbind(both,bothPts9)
master9<- master9[,-4] 

#enter dataframe in logistic regression code
x9 <-master9
colnames(x9) <- tolower(colnames(x9))
x9[,2] <- as.numeric(x9[,2])
x9[,3] <- as.numeric(x9[,3]) 

#set up matrix to correct for spatial autocorrelation in logistic regression for 1929 
x.matrix9 <- earth.dist(x9[,c("latitude", "longitude")])
x.matrix9 <- as.matrix(x.matrix9)
pcnm.x.matrix9 <- pcnm(x.matrix9)
x.env.pcnm9 <- cbind(x9, pcnm.x.matrix9$vectors[,1:round((length(which(pcnm.x.matrix9$values > 0)))/2)])

for(i in 4:ncol(x9)){
  
  x.env.pcnm9 <- cbind(x.env.pcnm9, scale(x.env.pcnm9[,i], scale = sd(x.env.pcnm9[,i], na.rm = TRUE)))
  colnames(x.env.pcnm9)[ncol(x.env.pcnm9)] <- paste("std", colnames(x9)[i], sep = ".")
  
}

colnames(x.env.pcnm9)

#1929 - Look at the column names in the data frame "x.env.pcnm9". 
#Put all the columns named "PCNM..." into the model below, 
#followed by all the environmental variable columns starting with "std..."
model.lrm9 <- lrm(x.env.pcnm9[,"cytotype"] ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + std.tmean9 + std.ppt9 + std.vpdmax9, y = TRUE, data=x.env.pcnm9, penalty = 0.001)

#view results, variables with a p-value[Pr(>|z|)] of <0.05 are significant
#the PCNM variables are adjusting for spatial autocorrelation, 
#and allow me to make bolder statements about the variables that are influencing cytotype distribution
summary(model.lrm9)

#1929 logistic regression model results
sink("logistic_regression_results/logistic-regression9.txt")#creates a text file called logistic-regression9.txt in your logistic_regression_results directory
print(model.lrm9)
print(summary(model.lrm9))
sink()

###Using PRISM 2011 weather data
##Instead of using an ANOVA, I will use a Logistic Regression
#independent continuous variables:weather layers
#dependent categoric variable:species/ploidy column

#to make the dataframe with columns: cytotype labels, latitude,longitude, environmental variables containing extracted data
#2011
master11<- cbind(both,bothPts11)
master11<- master11[,-4] 

#enter dataframe in logistic regression code
x11 <-master11
colnames(x11) <- tolower(colnames(x11))
x11[,2] <- as.numeric(x11[,2])
x11[,3] <- as.numeric(x11[,3])

#set up matrix to correct for spatial autocorrelation in logistic regression for 2011 
x.matrix11 <- earth.dist(x11[,c("latitude", "longitude")])
x.matrix11 <- as.matrix(x.matrix11)
pcnm.x.matrix11 <- pcnm(x.matrix11)
x.env.pcnm11 <- cbind(x11, pcnm.x.matrix11$vectors[,1:round((length(which(pcnm.x.matrix11$values > 0)))/2)])

for(i in 4:ncol(x11)){
  
  x.env.pcnm11 <- cbind(x.env.pcnm11, scale(x.env.pcnm11[,i], scale = sd(x.env.pcnm11[,i], na.rm = TRUE)))
  colnames(x.env.pcnm11)[ncol(x.env.pcnm11)] <- paste("std", colnames(x11)[i], sep = ".")
  
}

colnames(x.env.pcnm11)

#2011 - Look at the column names in the data frame "x.env.pcnm11". 
#Put all the columns named "PCNM..." into the model below, 
#followed by all the environmental variable columns starting with "std."
model.lrm11 <- lrm(cytotype ~ PCNM1 + PCNM2 + PCNM3 + PCNM4 + PCNM5 + PCNM6 + PCNM7 + PCNM8 + std.tmean11 + std.ppt11 + std.vpdmin11 + std.tdmean11, family = binomial(link = "logit"), data=x.env.pcnm11)

#view results, variables with a p-value[Pr(>|z|)] of <0.05 are significant
#the PCNM variables are adjusting for spatial autocorrelation, 
#and allow me to make bolder statements about the variables that are influencing cytotype distribution
summary(model.lrm11)

#2011 logistic regression model results
sink("logistic_regression_results/logistic-regression11.txt")#creates a text file called logistic-regression11.txt in your logistic_regression_results directory
print(model.lrm11)
print(summary(model.lrm11))
sink()

## principle component analysis(PCA) for 1929
cytotypes9 <- bothPts9[ ,1] #holds names of varieites and will be used in PCA plot
bothNum9 <- bothPts9[ ,-1] #remove species names
pca_both9 <- prcomp(bothNum9, center = TRUE, scale. = TRUE) #PCA
print(pca_both9) #print deviations and rotation
summary(pca_both9) #print importance of components
plot(pca_both9, type="l") #plot variances
ncomp <- 2 #specify number of components to load (representing 85% of variation)
#make a PCA plot
ggbiplot(pca_both9, obs.scale =1, var.scale = 1, groups= cytotypes9, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

## principle component analysis(PCA) for 2011
cytotypes11 <- bothPts11[ ,1] #holds names of varieites and will be used in PCA plot
bothNum11 <- bothPts11[ ,-1] #remove species names
pca_both11 <- prcomp(bothNum11, center = TRUE, scale. = TRUE) #PCA
print(pca_both11) #print deviations and rotations
summary(pca_both11) #print importance of components
plot(pca_both11, type="l") #plot variances
ncomp <- 3 #specify number of components to load (representing 99% of variation)
#make a PCA plot
ggbiplot(pca_both11, obs.scale =1, var.scale = 1, groups= cytotypes11, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')


## model-based approaches
# read in default maxent models
rDip9 <- raster("models/diploid1929.grd")
rTetra9 <- raster("models/tetraploid1929.grd")
rBoth9<- raster("models/both1929.grd")

# assessing niche overlap by comparing diploids and tetraploids in 1929
nicheOverlap(rDip9, rTetra9, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip9, rTetra9, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

# read in advanced maxent models
rDipAdv9 <- raster("models/diploidAdv1929.grd")
rTetraAdv9 <- raster("models/tetraploidAdv1929.grd")
rBothAdv9 <- raster("models/bothAdv1929.grd")

# assessing niche overlap by comparing diploids and tetraploids in 1929
nicheOverlap(rDipAdv9, rTetraAdv9, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDipAdv9, rTetraAdv9, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic


## model-based approaches
# read in default maxent models
rDip11 <- raster("models/diploid2011.grd")
rTetra11 <- raster("models/tetraploid2011.grd")
rBoth11<- raster("models/both2011.grd")

# assessing niche overlap by comparing diploids and tetraploids in 2011
nicheOverlap(rDip11, rTetra11, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip11, rTetra11, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

# read in advanced maxent models
rDipAdv11 <- raster("models/diploidAdv2011.grd")
rTetraAdv11 <- raster("models/tetraploidAdv2011.grd")
rBothAdv11 <- raster("models/bothAdv2011.grd")

# assessing niche overlap by comparing diploids and tetraploids in 2011
nicheOverlap(rDipAdv11, rTetraAdv11, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDipAdv11, rTetraAdv11, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in diploid niche from 1929 to 2011
nicheOverlap(rDip9, rDip11, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip9, rDip11, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in tetraploid niche from 1929 to 2011
nicheOverlap(rTetra9, rTetra11, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rTetra9, rTetra11, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in both cytotypes niche from 1929 to 2011
nicheOverlap(rBoth9, rBoth11, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rBoth9, rBoth11, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in diploid niche from 1929 to 2011
nicheOverlap(rDipAdv9, rDipAdv11, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDipAdv9, rDipAdv11, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in tetraploid niche from 1929 to 2011
nicheOverlap(rTetraAdv9, rTetraAdv11, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rTetraAdv9, rTetraAdv11, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

#assessing changes in both cytotypes niche from 1929 to 2011
nicheOverlap(rBothAdv9, rBothAdv11, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rBothAdv9, rBothAdv11, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic