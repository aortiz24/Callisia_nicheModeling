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

#making two objects for diploids that are permuted datasets: 
#x.permutedA contains half of the diploids occurrence points and will be run in maxent with 1929 layers in for loop
#x.permutedA2 contains half of the diploids occurrence points and will be run in maxent with 2011 layers in for loop
#assign 15 occurrence points from the diploid to the x.permutedA and do not replace the values
x.permutedA<-replicate(100, {sample(1:nrow(diploid), size = 15, replace = FALSE)
#contains the row names of the diploid object in numerical order. The information in these rows will be put into x.permutedA 
x.permutedA <- x.permutedA[order(x.permutedA)]
#put the remaining row names of the diploid object into x.permutedA2. 
x.permutedA2 <- setdiff(1:nrow(diploid), x.permutedA)
#contains the row names of the diploid object in numerical order. The information in these rows will be put into x.permutedA2. 
x.permutedA2 <- x.permutedA2[order(x.permutedA2)]
})
  
  ##For loop for diploids - 1929 vs 2011
  #one dataset will run 100 times with 1929 layers in maxent, and an I statistic will be calculated for each run
  #the other dataset will run 100 times with 2011 layers in maxent, and an I statistic will be calculated for each run
  #The critical value (the fifth lowest I statistic out of 100) will be used to conclude whether the niches are significantly different for 1929 and 2011
  sink("permutation_results/diploid_permut_vals.csv")#creates a csv file called diploid_permut_vals.csv in your permutation_results directory
  for (i in 1:100){
  
  #import specific rows of diploid locality data into x.permutedA and x.permutedA2
  #creates paired datasets
  x.permutedA.diploid1 <- diploid[(x.permutedA),]
  x.permutedA.diploid2 <- diploid[(x.permutedA2),]
  
  #1929-runing maxent with jackknife, random seed, and response curves, followed by cross-validation
  permutedMaxDipAdv9 <- maxent(
    x=predictors9,
    p=x.permutedA.diploid1,
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
  rPermutedMaxDipAdv9 <- predict(permutedMaxDipAdv9, predictors9)
  
  #2011-maxent with jackknife, random seed, and response curves, followed by cross-validation
  permutedMaxDipAdv11 <- maxent(
    x=predictors11,
    p=x.permutedA.diploid2,
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
  rPermutedMaxDipAdv11 <- predict(permutedMaxDipAdv11, predictors11)
  
  #calculate nicheOverlap I statistic
  DipAdvIstat<-nicheOverlap(rPermutedMaxDipAdv9, rPermutedMaxDipAdv11, stat='I', mask=TRUE, checkNegatives=TRUE)
  
  #print the 100 I statistics in the first column in diploid_permut_vals.csv
  print(DipAdvIstat)
}
sink()

#The critical value (the fifth lowest I statistic out of 100,you only got a value lower than this 5% of the time, P<0.05) will be used to conclude whether the niches are significantly different for 1929 and 2011
PermutIstats <- read.csv("permutation_results/diploid_permut_vals.csv")
#assign first 100 rows to the R object
PermutIstats <- PermutIstats[1:100,]
#ordering permuted 100 I statistic values from least to greatest
x<-sort(PermutIstats, decreasing = FALSE)
#writes the 100 I statistic values from least to greatest
write.csv(x, file="permutation_results/diploid_OrderedPermutIstats.csv")
#the critical value is the fifth lowest I statistic out of 100,
#you only get a value lower than this 5% of the time, P<0.05)
#When comparing the diploid niches in 1929 & 2011, the critical value is

##For loop for tetraploids - 1929 vs 2011
#one dataset will run 100 times with 1929 layers in maxent, and an I statistic will be calculated for each run
#the other dataset will run 100 times with 2011 layers in maxent, and an I statistic will be calculated for each run
#The critical value (the fifth lowest I statistic out of 100) will be used to conclude whether the niches are significantly different for 1929 and 2011
sink("permutation_results/tetraploid_permut_vals.csv")#creates a csv file called tetraploid_permut_vals.csv in your permutation_results directory
for (i in 1:100){
  #making two objects for tetraploids that are permuted datasets: 
  #x.permutedB contains half of the tetraploids occurrence points and will be run in maxent with 1929 layers in for loop
  #x.permutedB2 contains half of the tetraploids occurrence points and will be run in maxent with 2011 layers in for loop
  #assign 39 occurrence points from the tetraploid object to the x.permutedB and do not replace the values
  x.permutedB<-replicate(100, {sample(1:nrow(tetraploid), size = 39, replace = FALSE)
    #contains the row names of the tetraploid object in numerical order. The information in these rows will be put into x.permutedB. 
    x.permutedB <- x.permutedB[order(x.permutedB)]
    #put the remaining row names of the tetraploid object into x.permutedB2. 
    x.permutedB2 <- setdiff(1:nrow(tetraploid), x.permutedB)
    #contains the row names of the tetraploid object in numerical order. The information in these rows will be put into x.permutedB2. 
    x.permutedB2 <- x.permutedB2[order(x.permutedB2)]
  })
  
  #import specific rows of tetraploid locality data into x.permutedB and x.permutedB2
  #creates paired datasets
  x.permutedB.tetraploid1 <- tetraploid[(x.permutedB),]
  x.permutedB.tetraploid2 <- tetraploid[(x.permutedB2),]
  
  #1929-runing maxent with jackknife, random seed, and response curves, followed by cross-validation
  permutedMaxTetrAdv9 <- maxent(
    x=predictors9,
    p=x.permutedB.tetraploid1,
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
  rPermutedMaxTetrAdv9 <- predict(permutedMaxTetrAdv9, predictors9)
  
  #2011-maxent with jackknife, random seed, and response curves, followed by cross-validation
  permutedMaxTetrAdv11 <- maxent(
    x=predictors11,
    p=x.permutedB.tetraploid2,
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
  rPermutedMaxTetrAdv11 <- predict(permutedMaxTetrAdv11, predictors11)
  
  #calculate nicheOverlap I statistic
  TetrAdvIstat<-nicheOverlap(rPermutedMaxTetrAdv9, rPermutedMaxTetrAdv11, stat='I', mask=TRUE, checkNegatives=TRUE)
  
  #print the 100 I statistics in the first column in tetraploid_permut_vals.csv
  print(TetrAdvIstat)
}
sink()

#The critical value (the fifth lowest I statistic out of 100,you only got a value lower than this 5% of the time, P<0.05) will be used to conclude whether the niches are significantly different for 1929 and 2011
PermutIstats <- read.csv("permutation_results/tetraploid_permut_vals.csv")
#assign first 100 rows to the R object
PermutIstats <- PermutIstats[1:100,]
#ordering permuted 100 I statistic values from least to greatest
x<-sort(PermutIstats, decreasing = FALSE)
#writes the 100 I statistic values from least to greatest
write.csv(x, file="permutation_results/tetraploid_OrderedPermutIstats.csv")
#the critical value is the fifth lowest I statistic out of 100,
#you only get a value lower than this 5% of the time, P<0.05)
#When comparing the tetraploid niches in 1929 & 2011, the critical value is

##For loop for both cytotypes - 1929 vs 2011
#one dataset will run 100 times with 1929 layers in maxent, and an I statistic will be calculated for each run
#the other dataset will run 100 times with 2011 layers in maxent, and an I statistic will be calculated for each run
#The critical value (the fifth lowest I statistic out of 100) will be used to conclude whether the niches are significantly different for 1929 and 2011
sink("permutation_results/both_cytotypes_permut_vals.csv")#creates a csv file called both_cytotypes_permut_vals.csv in your permutation_results directory
for (i in 1:100){
  #making two objects for both cytotypes that are permuted datasets: 
  #x.permutedC contains half of the both cytotypes occurrence points and will be run in maxent with 1929 layers in for loop
  #x.permutedC2 contains half of the both cytotypes occurrence points and will be run in maxent with 2011 layers in for loop
  #assign 57 occurrence points from the both cytotypes object to the x.permutedC and do not replace the values
  x.permutedC<-replicate(100, {sample(1:nrow(both), size = 57, replace = FALSE)
    #contains the row names of the both cytotypes object in numerical order. The information in these rows will be put into x.permutedC. 
    x.permutedC <- x.permutedC[order(x.permutedC)]
    #put the remaining row names of the both cytotypes object into x.permutedC2. 
    x.permutedC2 <- setdiff(1:nrow(both), x.permutedC)
    #contains the row names of the both cytotypes object in numerical order. The information in these rows will be put into x.permutedC2. 
    x.permutedC2 <- x.permutedC2[order(x.permutedC2)]
  })
  
  #import specific rows of both cytotypes locality data into x.permutedC and x.permutedC2
  #creates paired datasets
  x.permutedC.both1 <- both[(x.permutedC),]
  x.permutedC.both2 <- both[(x.permutedC2),]
  
  #1929-runing maxent with jackknife, random seed, and response curves, followed by cross-validation
  permutedMaxBothAdv9 <- maxent(
    x=predictors9,
    p=x.permutedC.both1,
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
    p=x.permutedC.both2,
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
#When comparing the both cytotypes niches in 1929 & 2011, the critical value is
