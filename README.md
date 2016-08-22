# Callisia_nicheModeling

Scripts for creating maps and performing niche modeling in RStudio.

The general workflow is:
•	mapping.R: importing Callisia graminea occurrence data into R, creating simple maps, creating customized shapefiles
•	layerPrep.R: masking/clipping BioClim layers, looking for correlations between layers
•	maxent.R: creating niche models using occurrence data and climate layers
•	nicheOverlap.R: assessing whether niche models for different cytotypes are distinct from each other

There are a number of files and directories that are downloaded or created during this analysis. To perform all commands in layerPrep.R, you should download the BioClim layers and store them in a convenient place (note the paths for using these layers may need to be changed). Clipped layers are included in layers/ for your convenience.
