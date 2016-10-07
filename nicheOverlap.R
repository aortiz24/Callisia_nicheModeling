## evaluating niche overlap between taxa

library(dplyr)
library(raster)
library(dismo)
library(ecospat)
library(ENMeval)

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

# extract layer data for each point and add label
dipPts <- raster::extract(predictors, diploid)
dipPts <- cbind.data.frame(species="diploid", dipPts) #add column for diploid
dipPts0 <- raster::extract(predictors0, diploid)
dipPts0 <- cbind.data.frame(species="diploid", dipPts0) #add column for diploid
dipPts1 <- raster::extract(predictors1, diploid)
dipPts1 <- cbind.data.frame(species="diploid", dipPts1) #add column for diploid
tetraPts <- raster::extract(predictors, tetraploid)
tetraPts <- cbind.data.frame(species="tetraploid", tetraPts)
tetraPts<-na.omit(tetraPts)#removing NA values
tetraPts0 <- raster::extract(predictors0, tetraploid)
tetraPts0 <- cbind.data.frame(species="tetraploid", tetraPts0)
tetraPts0<-na.omit(tetraPts0)#removing NA values
tetraPts1 <- raster::extract(predictors1, tetraploid)
tetraPts1 <- cbind.data.frame(species="tetraploid", tetraPts1)
tetraPts1<-na.omit(tetraPts1)#removing NA values

# combine diploid and tetraploid
bothPts <- as.data.frame(rbind(dipPts, tetraPts))
bothPts0 <- as.data.frame(rbind(dipPts0, tetraPts0))
bothPts1 <- as.data.frame(rbind(dipPts1, tetraPts1))

# one-way ANOVA with Tukey's post-hoc (example from altitude)
aov.alt <- aov(alt ~ species, data=bothPts)
summary(aov.alt)
TukeyHSD(aov.alt)

###Using Bioclim weather data
##for loop of one-way ANOVA with Tukey's post-hoc(for all 11 uncorrelated weather variables)
bothPts <- as.data.frame(rbind(dipPts, tetraPts))#save dataset(made previously in script)as object for ANOVA analysis
bothPts #view dataset layout
1:ncol(bothPts) #displays how many columns are in dataset
AVz<- rep(NA,ncol(bothPts)) #creates a table called AVz with the same number of columns as the dataset. When it is created each cell will have an NA, then we will add data from the for loop in this table.
sink("anova_results/ANOVA-Tukey.txt")#creates a text file called ANOVA-Tukey.txt in your anova_results directory
for (i in 2:ncol(bothPts)) {
    column <-names(bothPts[i])
    AVz<-summary(aov(bothPts[,i]~species, data=bothPts))
    tk<-TukeyHSD((aov(bothPts[,i]~species, data=bothPts)))
print(column)
print(AVz)
print(tk)
}
sink()

# principle component analysis(PCA)
bothNum <- bothPts[ ,-1] #remove species names
pca_both <- prcomp(bothNum, center = TRUE, scale. = TRUE) #PCA
print(pca_both) #print deviations and rotations
summary(pca_both) #print importance of components
plot(pca_both, type="l") #plot variances
ncomp <- 8 #specify number of components to load (representing 99% of variation)

## model-based approaches
# read in default maxent models
rDip <- raster("models/diploid.grd")
rTetra <- raster("models/tetraploid.grd")
# assessing niche overlap
nicheOverlap(rDip, rTetra, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip, rTetra, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

###Using PRISM 1930 weather data
##for loop of one-way ANOVA with Tukey's post-hoc(for all uncorrelated weather variables)
bothPts0 <- as.data.frame(rbind(dipPts0, tetraPts0))#save dataset(made previously in script)as object for ANOVA analysis
bothPts0 #view dataset layout
1:ncol(bothPts0) #displays how many columns are in dataset
AVz0<- rep(NA,ncol(bothPts0)) #creates a table called AVz with the same number of columns as the dataset. When it is created each cell will have an NA, then we will add data from the for loop in this table.
sink("anova_results/ANOVA-Tukey0.txt")#creates a text file called ANOVA-Tukey.txt in your anova_results directory
for (i in 2:ncol(bothPts0)) {
  column0 <-names(bothPts0[i])
  AVz0<-summary(aov(bothPts0[,i]~species, data=bothPts0))
  tk0<-TukeyHSD((aov(bothPts0[,i]~species, data=bothPts0)))
  print(column0)
  print(AVz0)
  print(tk0)
}
sink()

# principle component analysis(PCA)
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
# assessing niche overlap
nicheOverlap(rDip0, rTetra0, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip0, rTetra0, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

###Using PRISM 2014 weather data
##for loop of one-way ANOVA with Tukey's post-hoc(for all uncorrelated weather variables)
bothPts1 <- as.data.frame(rbind(dipPts1, tetraPts1))#save dataset(made previously in script)as object for ANOVA analysis
bothPts1 #view dataset layout
1:ncol(bothPts1) #displays how many columns are in dataset
AVz1<- rep(NA,ncol(bothPts1)) #creates a table called AVz with the same number of columns as the dataset. When it is created each cell will have an NA, then we will add data from the for loop in this table.
sink("anova_results/ANOVA-Tukey1.txt")#creates a text file called ANOVA-Tukey.txt in your anova_results directory
for (i in 2:ncol(bothPts1)) {
  column1 <-names(bothPts1[i])
  AVz1<-summary(aov(bothPts1[,i]~species, data=bothPts1))
  tk1<-TukeyHSD((aov(bothPts1[,i]~species, data=bothPts1)))
  print(column1)
  print(AVz1)
  print(tk1)
}
sink()

# principle component analysis(PCA)
bothNum1 <- bothPts1[ ,-1] #remove species names
pca_both1 <- prcomp(bothNum1, center = TRUE, scale. = TRUE) #PCA
print(pca_both1) #print deviations and rotations
summary(pca_both1) #print importance of components
plot(pca_both1, type="l") #plot variances
ncomp <- 2 #specify number of components to load (representing 99% of variation)

## model-based approaches
# read in default maxent models
rDip1 <- raster("models/diploid2014.grd")
rTetra1 <- raster("models/tetraploid2014.grd")
# assessing niche overlap
nicheOverlap(rDip1, rTetra1, stat='D', mask=TRUE, checkNegatives=TRUE) # D statistic
nicheOverlap(rDip1, rTetra1, stat='I', mask=TRUE, checkNegatives=TRUE) # I statistic

# assessing niche equivalency
#nicheEquivalency()
