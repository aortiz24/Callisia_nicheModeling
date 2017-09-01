###For checking if the species occurrence points(diploid1) are in the states(SEstates)

#load library
library(sp)
library(dplyr)

# import occurrence data and convert to format required by maxent
Callisia.both <- read.csv(file="CallisiaCompletedData.csv") %>%
  select(Cytotype,Latitude,Longitude)
Callisia.both <- na.omit(Callisia.both)
diploid1 <- Callisia.both %>%
  filter(Cytotype=="2X")
diploid1 <- diploid1[,c(2,3)]
tetraploid1 <- Callisia.both %>%
  filter(Cytotype=="4X")
tetraploid1 <- tetraploid1[,c(2,3)]


##making a Spatial Points DataFrame using the statistical function notation(with a tilde) for diploids
coordinates(diploid1) <-~Longitude+Latitude
crs(diploid1)<-crs(SEstates)
class(diploid1)
#You can use the coordinates to do a spatial query of the polygons in SEstates
ovr<- over(diploid1, SEstates)
#for each point 'ovr' has the matching record from SEstates
head(ovr)
#State names in SEstates dataset in diploid1 object
cntr<-ovr$NAME
#Ask two questions,
#Which points do not match any state?
v<-which(is.na(cntr))
v

##making a Spatial Points DataFrame using the statistical function notation(with a tilde) for tetraploids
coordinates(tetraploid1) <-~Longitude+Latitude
crs(tetraploid1)<-crs(SEstates)
class(tetraploid1)
#You can use the coordinates to do a spatial query of the polygons in SEstates
ovr1<- over(tetraploid1, SEstates)
#for each point 'ovr' has the matching record from SEstates
head(ovr1)
#State names in SEstates dataset in diploid1 object
cntr1<-ovr1$NAME
#Ask two questions,
#Which points do not match any state?
v<-which(is.na(cntr1))
v
#delete rows 10,46,56 from tetraploid1 object and rerun code with tetraploid1 to confirm that there are no points outside of SEstates boundaries
tetraploid1<- tetraploid1[-c(10,46,56), ]
