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

# assessing niche equivalency
#nicheEquivalency()
