## evaluating niche overlap between taxa

library(dplyr)
library(raster)
library(dismo)
library(ecospat)
library(ENMeval)

## multi-variate climate space comparisons (standard statistical testing, non-model based)
# import occurrence data and convert to format required by maxent
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

# extract layer data for each point and add label
dipPts <- raster::extract(predictors, diploid)
dipPts <- cbind.data.frame(species="diploid", dipPts) #add column for diploid
tetraPts <- raster::extract(predictors, tetraploid)
tetraPts <- cbind.data.frame(species="tetraploid", tetraPts)
# combine diploid and tetraploid
bothPts <- as.data.frame(rbind(dipPts, tetraPts))

# one-way ANOVA with Tukey's post-hoc (example from altitude)
aov.alt <- aov(alt ~ species, data=bothPts)
summary(aov.alt)
TukeyHSD(aov.alt)

#for loop of one-way ANOVA with Tukey's post-hoc(for all uncorrelated weather variables)
e<-bothPts
for(i in e[,2:12]) {aov(i~ species,data=bothPts)}

# principle component analysis
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

# assessing niche equivalency
#nicheEquivalency()