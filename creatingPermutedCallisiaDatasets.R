###Creating permuted Callisia Datasets

library(dplyr)
library(raster)
library(dismo)
library(ENMeval)

# create directory for saving the permutation results later
dir.create("permutation_results")

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

##For loop for both cytotypes - 1929 vs 2011
#one dataset will run 100 times with 1929 layers in maxent, and an I statistic will be calculated for each run
#the other dataset will run 100 times with 2011 layers in maxent, and an I statistic will be calculated for each run
#The critical value (the fifth lowest I statistic out of 100) will be used to conclude whether the niches are significantly different for 1929 and 2011
sink("permutation_results/both_cytotypes_permut_vals.csv")#creates a text file called both_cytotypes_permut_vals.csv in your permutation_results directory
for (i in 1:100){
  #making two objects for both cytotypes that are permuted datasets: 
  #x.permuted.object contains half of the both cytotypes occurrence points and will be run in maxent with 1929 layers in for loop
  #x.permuted.object2 contains half of the both cytotypes occurrence points and will be run in maxent with 2011 layers in for loop
  #assign 57 occurrence points from the both cytotypes object to the x.permuted object and do not replace the values
  x.permuted<-replicate(100, {sample(1:nrow(both), size = 57, replace = FALSE)
    #contains the row names of the both cytotypes object in numerical order. The information in these rows will be put into x.permuted. 
    x.permuted <- x.permuted[order(x.permuted)]
    #put the remaining row names of the both cytotypes object into x.permuted2. 
    x.permuted2 <- setdiff(1:nrow(both), x.permuted)
    #contains the row names of the both cytotypes object in numerical order. The information in these rows will be put into x.permuted2. 
    x.permuted2 <- x.permuted2[order(x.permuted2)]
  })
  
  #import specific rows of both cytotypes locality data into x.permuted.object and x.permuted.object2
  #creates paired datasets
  x.permuted.both1 <- both[(x.permuted),]
  x.permuted.both2 <- both[(x.permuted2),]
  
  #1929-runing maxent with jackknife, random seed, and response curves, followed by cross-validation
  permutedMaxBothAdv9 <- maxent(
    x=predictors9,
    p=x.permuted.both1,
    removeDuplicates=TRUE,
    nbg=10000,
    args=c(
      'randomseed=true', #default=false
      'threads=2', #default=1
      'responsecurves=true', #default=false
      'jackknife=true', #default=false
      'replicatetype=crossvalidate',
      'maximumiterations=1000' #default=500
    )
  )
  #1929-creating model
  rPermutedMaxBothAdv9 <- predict(permutedMaxBothAdv9, predictors9)
  
  #2011-maxent with jackknife, random seed, and response curves, followed by cross-validation
  permutedMaxBothAdv11 <- maxent(
    x=predictors11,
    p=x.permuted.both2,
    removeDuplicates=TRUE,
    nbg=10000,
    args=c(
      'randomseed=true', #default=false
      'threads=2', #default=1
      'responsecurves=true', #default=false
      'jackknife=true', #default=false
      'replicatetype=crossvalidate',
      'maximumiterations=1000' #default=500
    )
  )
  #2011-creating model
  rPermutedMaxBothAdv11 <- predict(permutedMaxBothAdv11, predictors11)
  
  #calculate nicheOverlap I statistic
  BothAdvIstat<-nicheOverlap(rPermutedMaxBothAdv9, rPermutedMaxBothAdv11, stat='I', mask=TRUE, checkNegatives=TRUE)
  
  #print the 100 I statistics in the first column in both_cytotypes_permut_vals.csv
  print(BothAdvIstat)
}
sink()

#The critical value (the fifth lowest I statistic out of 100,you only got a value lower than this 5% of the time, P<0.05) will be used to conclude whether the niches are significantly different for 1929 and 2011
PermutIstats <- read.csv("permutation_results/both_cytotypes_permut_vals.csv")
#assign first 100 rows to the R object
PermutIstats <- PermutIstats[1:100,]
#ordering permuted 100 I statistic values from least to greatest
x<-sort(PermutIstats, decreasing = FALSE)
#writes the 100 I statistic values from least to greatest
write.csv(x, file="permutation_results/both_cytotypes_OrderedPermutIstats.csv")
#the critical value is the fifth lowest I statistic out of 100,
#you only get a value lower than this 5% of the time, P<0.05)
#When comparing the both cytotypes niches in 1929 & 2011, the critical value is 0.9204252
