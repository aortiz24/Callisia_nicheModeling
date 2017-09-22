###Constructing figures for Thesis
#load libraries
library(plyr)
library(gridExtra)

##make test gain table for Callisia
#for diploid in 1929
dip9<- read.csv("models/diploid1929Maxent/maxentResults.csv")
dip9<- dip9[ ,c(9,8,27:29)]

#for tetraploid in 1929
tetra9<- read.csv("models/tetraploid1929Maxent/maxentResults.csv")
tetra9<- tetra9[ ,c(9,8,27:29)]

#for both cytotypes in 1929
both9<- read.csv("models/both1929Maxent/maxentResults.csv")
both9<- both9[ ,c(9,8,27:29)]

#combine rows from all result tables into 1929 test gain table
testGainTable1929 <- rbind(dip9,tetra9,both9)

#for diploid in 2011
dip11<- read.csv("models/diploid2011Maxent/maxentResults.csv")
dip11<- dip11[ ,c(9,8,32:35)]

#for tetraploid in 2011
tetra11<- read.csv("models/tetraploid2011Maxent/maxentResults.csv")
tetra11<- tetra11[ ,c(9,8,32:35)]

#for hybrids in 2011
both11<- read.csv("models/both2011Maxent/maxentResults.csv")
both11<- both11[ ,c(9,8,32:35)]

#combine rows from all result tables into 2011 test gain table
testGainTable2011 <- rbind(dip11,tetra11,both11)

##making publicatoin ready test gain table
#create Taxon column for testGainTable1929 & testGainTable2011
Taxon<-data.frame(Taxon = c("diploids","tetraploids","both cytotypes"))
#add column to testGainTable1929 & testGainTable2011
FinalTestGainTable1929<-cbind(Taxon,testGainTable1929)
FinalTestGainTable2011<-cbind(Taxon,testGainTable2011)

#using plyr,rename columns for testGainTable1929 & testGainTable2011
FinalTestGainTable1929<-rename(FinalTestGainTable1929,c("Test.AUC"="Test AUC","Test.gain"="Full model","Test.gain.with.only.ppt9"="Only precipitation","Test.gain.with.only.tmean9"="Only mean temperature","Test.gain.with.only.vpdmin9"="Only minimum vapor pressure deficit"))
FinalTestGainTable2011<-rename(FinalTestGainTable2011,c("Test.AUC"="Test AUC","Test.gain"="Full model","Test.gain.with.only.ppt11"="Only precipitation","Test.gain.with.only.tmean11"="Only mean temperature","Test.gain.with.only.vpdmax11"="Only maximum vapor pressure deficit","Test.gain.with.only.tdmean11"="Only mean dewpoint temperature"))

#export tables
#1929
pdf(file = "figures/TestGainTable1929.pdf",height = 2, width = 11)
grid.table(FinalTestGainTable1929)
dev.off()

#2011
pdf(file = "figures/TestGainTable2011.pdf",height = 2, width = 14)
grid.table(FinalTestGainTable2011)
dev.off()