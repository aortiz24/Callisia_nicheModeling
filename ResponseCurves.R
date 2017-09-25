###Response Curve Graphs for Callisia graminea
##Diploids 1929 - model must be loaded
# create objects for each variable in model
tDip1 <- response(maxDipAdv9, var = 1) #mean temperature
pDip2 <- response(maxDipAdv9, var = 2) #precipitation
vDip3 <- response(maxDipAdv9, var = 3) #maximum vapor pressure deficit

#open png file
png(filename="figures/diploids_response_curves_1929.png")
# combined figure
par(mfrow=c(1,3))
#plot mean temperature 1929 layer
plot(tDip1, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Mean Temperature (Celsius)", main = "(a)")
#plot precipitation 1929 layer
plot(pDip2, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Precipitation (mm)", main = "(b)")
#plot maximum vapor pressure deficit layer
plot(vDip3, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Maximum Vapor Pressure Deficit", main = "(c)")
dev.off()

##Tetraploids 1929 - model must be loaded
# create objects for each variable in model
tTetra1 <- response(maxTetraAdv9, var = 1) #mean temperature
pTetra2 <- response(maxTetraAdv9, var = 2) #precipitation
vTetra3 <- response(maxTetraAdv9, var = 3) #maximum vapor pressure deficit

#open png file
png(filename="figures/tetraploids_response_curves_1929.png")
# combined figure
par(mfrow=c(1,3))
#plot mean temperature 1929 layer
plot(tTetra1, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Mean Temperature (Celsius)", main = "(a)")
#plot precipitation 1929 layer
plot(pTetra2, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Precipitation (mm)", main = "(b)")
#plot maximum vapor pressure deficit layer
plot(vTetra3, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Maximum Vapor Pressure Deficit", main = "(c)")
dev.off()

##Both Cytotypes 1929 - model must be loaded
# create objects for each variable in model
tBoth1 <- response(maxBothAdv9, var = 1) #mean temperature
pBoth2 <- response(maxBothAdv9, var = 2) #precipitation
vBoth3 <- response(maxBothAdv9, var = 3) #maximum vapor pressure deficit

#open png file
png(filename="figures/BothCytotypes_response_curves_1929.png")
# combined figure
par(mfrow=c(1,3))
#plot mean temperature 1929 layer
plot(tBoth1, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Mean Temperature (Celsius)", main = "(a)")
#plot precipitation 1929 layer
plot(pBoth2, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Precipitation (mm)", main = "(b)")
#plot maximum vapor pressure deficit layer
plot(vBoth3, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Maximum Vapor Pressure Deficit", main = "(c)")
dev.off()

##Diploids 2011 - model must be loaded
# create objects for each variable in model
tDip1 <- response(maxDipAdv11, var = 1) #mean temperature
pDip2 <- response(maxDipAdv11, var = 2) #precipitation
vDip3 <- response(maxDipAdv11, var = 3) #minimum vapor pressure deficit
dDip4 <- response(maxDipAdv11, var = 4) #mean dewpoint temperature

#open png file
png(filename="figures/diploids_response_curves_2011.png")
# combined figure
par(mfrow=c(2,2))
#plot mean temperature 2011 layer
plot(tDip1, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Mean Temperature (Celsius)", main = "(a)")
#plot precipitation 2011 layer
plot(pDip2, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Precipitation (mm)", main = "(b)")
#plot minimum vapor pressure deficit layer
plot(vDip3, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Minimum Vapor Pressure Deficit", main = "(c)")
#plot mean dewpoint temperature layer
plot(dDip4, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Mean Dewpoint Temperature (Celsius)", main = "(d)")
dev.off()

##Tetraploids 2011 - model must be loaded
# create objects for each variable in model
tTetra1 <- response(maxTetraAdv11, var = 1) #mean temperature
pTetra2 <- response(maxTetraAdv11, var = 2) #precipitation
vTetra3 <- response(maxTetraAdv11, var = 3) #minimum vapor pressure deficit
dTetra4 <- response(maxTetraAdv11, var = 4) #mean dewpoint temperature

#open png file
png(filename="figures/tetraploids_response_curves_2011.png")
# combined figure
par(mfrow=c(2,2))
#plot mean temperature 2011 layer
plot(tTetra1, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Mean Temperature (Celsius)", main = "(a)")
#plot precipitation 2011 layer
plot(pTetra2, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Precipitation (mm)", main = "(b)")
#plot minimum vapor pressure deficit layer
plot(vTetra3, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Minimum Vapor Pressure Deficit", main = "(c)")
#plot mean dewpoint temperature layer
plot(dTetra4, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Mean Dewpoint Temperature (Celsius)", main = "(d)")
dev.off()

##Both Cytotypes 2011 - model must be loaded
# create objects for each variable in model
tBoth1 <- response(maxBothAdv11, var = 1) #mean temperature
pBoth2 <- response(maxBothAdv11, var = 2) #precipitation
vBoth3 <- response(maxBothAdv11, var = 3) #minimum vapor pressure deficit
dBoth4 <- response(maxBothAdv11, var = 4) #mean dewpoint temperature

#open png file
png(filename="figures/BothCytotypes_response_curves_2011.png")
# combined figure
par(mfrow=c(2,2))
#plot mean temperature 2011 layer
plot(tBoth1, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Mean Temperature (Celsius)", main = "(a)")
#plot precipitation 2011 layer
plot(pBoth2, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Precipitation (mm)", main = "(b)")
#plot minimum vapor pressure deficit layer
plot(vBoth3, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Minimum Vapor Pressure Deficit", main = "(c)")
#plot mean dewpoint temperature layer
plot(dBoth4, type = "l", col="red", lwd=5, ylim=c(0,1), ylab = "Habitat Suitability", xlab = "Mean Dewpoint Temperature (Celsius)", main = "(d)")
dev.off()
