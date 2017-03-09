## evaluating niche overlap between taxa

library(dplyr)
library(raster)
library(dismo)
library(ecospat)
library(ENMeval)
library(permute)

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

#import binary data for logistic regression
#0 equals diploid
#1 equals tetraploid
binary<- read.csv("CallisiaBinaryData.csv")
binary <- na.omit(binary)
diploidbin <- binary %>%
  filter(Cytotype=="0")
diploidbin <- diploidbin[,c(3,2)]
tetraploidbin <- binary %>%
  filter(Cytotype=="1")
tetraploidbin <- tetraploidbin[,c(3,2)]
#deleting rows whose points are outside of SEstates object
tetraploidbin<- tetraploidbin[-c(10,46,56), ]

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

# extract layer data for each point and add label for logistic regression
dipbinPts0 <- raster::extract(predictors0, diploidbin)
dipbinPts0 <- cbind.data.frame(species="0", dipbinPts0) #add column for diploid
dipbinPts1 <- raster::extract(predictors1, diploidbin)
dipbinPts1 <- cbind.data.frame(species="0", dipbinPts1) #add column for diploid
tetrabinPts0 <- raster::extract(predictors0, tetraploidbin)
tetrabinPts0 <- cbind.data.frame(species="1", tetrabinPts0)
tetrabinPts0<-na.omit(tetrabinPts0)#removing NA values
tetrabinPts1 <- raster::extract(predictors1, tetraploidbin)
tetrabinPts1 <- cbind.data.frame(species="1", tetrabinPts1)
tetrabinPts1<-na.omit(tetrabinPts1)#removing NA values

# combine diploid and tetraploid for logistic regression
bothbinPts0 <- as.data.frame(rbind(dipbinPts0, tetrabinPts0))
bothbinPts1 <- as.data.frame(rbind(dipbinPts1, tetrabinPts1))


###Using PRISM 1930 weather data
##Instead of using an ANOVA, I will use a Logistic Regression
#independent continuous variables:weather layers
#dependent categoric variable:species/ploidy column

#treat species as a categorical variable
bothbinPts1$species <- factor(bothbinPts1$species)

# develop testing and training sets for both cytotypes in 2014
foldbin <- kfold(bothbinPts1, k=5) #split occurence points into 5 sets
bothbin.Test1 <- bothbinPts1[foldbin == 1, ] #take 20% (1/5) for testing
bothbin.Train1 <- bothbinPts1[foldbin != 1, ] #leave 80% for training

#fit logistic regression model to data set
model<-glm(formula= species~ avg_2014_ppt0 + avg_2014_tmax0, data=bothbin.Train1, family="binomial")
#view summary of the results of the logistic regression model
summary(model)

#predict with logistic regression using test data
pred<- predict(model,bothbin.Test1, type="response")
#Misclassification error=when the model gets it wrong
model_pred_species<- rep("0",22)
model_pred_species[pred>0.5]<-"1"
tab<-table(model_pred_species, bothbin.Test1$species)
#predictive values are the numbers in the left column(0 = predicting diploid, 1 = predicting tetraploid)
#the top row of numbers represent the diploids (0) and tetrapliods (1)
# 3 plants are diploids and 15 are tetraploids
#18 correct predictions were made of the 22 possible predictions of the test data
print(tab) 

#1+3=4 predictions are missclasifications
#calculating misclassification rate=the lower, the better
#for this specific case it is about 0.18
1-sum(diag(tab))/sum(tab)

###Using PRISM 2014 weather data
##Instead of using an ANOVA, I will use a Logistic Regression
#independent continuous variables:weather layers
#dependent categoric variable:species/ploidy column

# develop testing and training sets for both cytotypes in 2014
fold <- kfold(bothPts1, k=5) #split occurence points into 5 sets
both.Test1 <- bothPts1[fold == 1, ] #take 20% (1/5) for testing
both.Train1 <- bothPts1[fold != 1, ] #leave 40% for training

#fit logistic regression model to data set
model<-glm(formula= species~., data=both.Train1, family=binomial(link='logit'))
#view results
model 
#view summary of the results of the logistic regression model
summary(model)
##Interpreting results of logistic regression model
#run an ANOVA on the model to analyze the table of deviance
anova(model, test="Chisq")
#McFadden R2 index can be used to assess the model fit
pR2(model)
##Assessing predictive ability of the model
#evaluating the fitting of the model
fitted.results <- predict(model,newdata=subset(both.Test1,select=c(2,3)),type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0) #parameters can be changed
misClasificError <- mean(fitted.results != both.Test1$species)
print(paste('Accuracy',1-misClasificError))
#plot the ROC curve 
p <- predict(model, newdata=both.Test1, type="response")
pr <- prediction(p, both.Test1$species)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
#calculate the AUC 
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

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

