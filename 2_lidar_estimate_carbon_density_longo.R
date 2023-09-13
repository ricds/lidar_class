# Laboratory of estimating carbon density using LiDAR CHM and Longo 2016 -----------------------------
# Ricardo Dalagnol - 12 Sep 2023
# R version 4.2.1

## Introduction
## Forest carbon density can be estimated using the LiDAR Canopy Height Model (CHM) following the paper of Longo et al. 2016. Aboveground biomass variability across intact and degraded forests in the Brazilian Amazon. Global Biogeochemical Cycles. 10.1002/2016GB005465
## Across the script I use one # for comments on the step being done and ## for additional explanations

## This laboratory covers:
# 1) Creating the CHM as required by the Longo et al. 2016 approach. Here I use the lidR package for the whole processing, but you can make the CHM with other softwares as well and jump to step 2
# 2) Estimating carbon density
## Let's go !

# load libraries
library(lidR) # install.packages("lidR")
library(sp)
library(raster)
library(rgdal)
library(rgl)

# set working directory to the directory of the script
setwd("E:\\lidar_workshop\\")

# set las filename
## This example .las file is a point cloud tile of ~250 x 250 meters
las_fname = "ORIG_NP_T-0262__638250_8876250.las"

# colors for plotting
hgt_col = height.colors(50)





# 1) Creating the CHM as required by the Longo et al. 2016 approach -------
## We first need a traditional CHM (maximum height within the cells) with 1-m resolution, then we re-sample the raster to 50x50 m taking the MEAN value

# read the las file in R inside the 'las' object
las = readLAS(las_fname)
plot(las, legend=T)

# classify ground points
ws = seq(3, 21, 3)
th = seq(0.1, 1.5, length.out = length(ws))
las = classify_ground(las, algorithm = pmf(ws,th))
plot(las, legend=T, color = "Classification")

# normalize the point cloud by subtracting the continuous-DTM (TIN)
nlas = normalize_height(las, tin())

# now we can filter for points with very high height, which are not trees
nlas = filter_poi(nlas, Z > 0 & Z < 100)
plot(nlas, legend=T)

# create a simple 1x1 CHM
## here we don't use na.fill so the borders are not much affected by artificial values
chm = grid_canopy(nlas, res = 1, algorithm = p2r(subcircle = 0.15))
plot(chm, col = hgt_col)

# re-sample the 1-m CHM to 50-m by taking the average value
## this is how the method works, it is based on 50x50 m data
chm_50 = aggregate(chm, fact = 50, fun = mean, na.rm=T)
chm_50 = crop(chm_50, chm) ## This is just so the chm_50 and chm fit exactly the same extent - You can try not running this and check in the overlay what happens!
chm_50
## Check the resolution

# Visualize one on top of the other just to check the overlay
plot(chm_50, main = "50m", col = hgt_col)
plot(chm, main = "1m", add=T, col = hgt_col)

# Visualize one each time
plot(chm, main = "1m", col = hgt_col)
plot(chm_50, main = "50m", col = hgt_col)






# 2) Estimating carbon density --------------------------------------------

# function to calculate mean and SD aboveground carbon density from live trees based on Longo et al 2016 paper Equation S1
# Longo et al. 2016. Aboveground biomass variability across intact and degraded forests in the Brazilian Amazon. Global Biogeochemical Cycles. 10.1002/2016GB005465
acd_calc_func = function(x) {
  mean_acd = 0.054 * (x ^ 1.74)
  sd_acd = 0.17 * (x ^ 1.04)
  y = stack(mean_acd, sd_acd)
  names(y) = c("Mean_ACD", "SD_ACD")
  return(y)
}

# apply the function to calc the ACD
acd = acd_calc_func(chm_50)
plot(acd, col = hgt_col)

# calc summary metrics over the whole area
summary(acd$Mean_ACD[])
summary(acd$SD_ACD[])

## This forest has XXXX +/- XXXXX kg C / m2






# Summary - Read this only after you are done with each step --------------


## Ok, after all this work, let's summarize what we've learned:

# 1) Creating the CHM as required by the Longo et al. 2016 approach
# It requires us to create a traditional 1-m CHM then resample it to 50x50.
# This is the way they modeled the relationship between field and lidar data.
#
# 2) Estimating carbon density: it is really simple to apply the equation.
# This Eq. S7 is for live trees ACD. There are other equations in their paper, including a more complex one
# that uses other LiDAR metrics beyond the CHM ! Perhaps you want to try and do that? :)
#

## Thank you and good luck!

# Ricardo Dalagnol
# ricds@hotmail.com



#
# The End -----------------------------------------------------------------
#

