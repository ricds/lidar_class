# Laboratory of LiDAR CHM processing using lidR -----------------------------
# Ricardo Dalagnol - 12 Sep 2023
# R version 4.2.1

## Introduction
## During the development of this laboratory, I followed the overall guidelines of the lidR book (https://jean-romain.github.io/lidRbook)
## Across the script I use one # for comments on the step being done and ## for additional explanations

## This laboratory covers:
# 1) Point cloud Reading and Plotting
# 2) Ground classification
# 3) Creating a Digital Terrain Model (DTM)
# 4) Height Normalization
# 5) Creating a Canopy Height Model (CHM)
# 6) LiDAR Metrics
## Let's go !

# load or install packages if needed
neededPackages = c("lidR", "sp", "raster", "rgdal", "rgl", "mapview", "gstat")
pkgTest = function(x) { if (x %in% rownames(installed.packages()) == FALSE) { install.packages(x, dependencies= TRUE) }; library(x, character.only = TRUE) }
for (package in neededPackages) { pkgTest(package) }

# set working directory to the directory of the script
setwd("E:\\lidar_workshop\\")

# set las filename
## This example .las file is a point cloud tile of ~250 x 250 meters
## .las and .laz are the main file types to store point cloud data, the 'z' means it is 'zipped'
las_fname = "ORIG_NP_T-0262__638250_8876250.las"

# colors for plotting
hgt_col = height.colors(50)


# 1) Point cloud Reading and Plotting -----------------------------------------

# read the las file in R inside the 'las' object
las = readLAS(las_fname)

# print the basic spatial information of the 'las' object
## Note that it says the CRS (coordinate reference system) of this data, area, total number of points, and density of points
## This is a small area and we have >381.000 points !
print(las)

# print the LAS Header
## This header is standard for las files and can be seen in any software that works with LAS files
## You can check additional info 
print(las@header)

# visualize the data frame inside the las file
## Note that the Classification column is all zeroed, i.e. points are not yet classified (ground, vegetation, etc.)
## You can see additional info on the LAS format in http://www.asprs.org/wp-content/uploads/2010/12/LAS_1_4_r13.pdf
View(las@data)

# plot the point cloud
## This can take a while in slower computers OR if the las file is big. In this example the Las file is not big.
## I suggest maximizing after plotting for better visualization
plot(las, legend=T)

# plot the point cloud with intensity overlay
## What are the differences?
plot(las, color = "Intensity", legend=T)

# load a las catalog object
## The LAScatalog function is meant for optimize processing, but we will use it here just to visualize where is the plot in a map
ctg = readLAScatalog(las_fname)

# visualize where the plot is located in a larger map
## Zoom out after it loads because its a small tile
plot(ctg, mapview = TRUE, map.type = "Esri.WorldImagery")



# 2) Ground classification --------------------------------------------------------------
## LASTools uses a variation of the Axellson 2000 TIN refinement algorithm. Axelsson P (2000) DEM generation from Laser Scanner data using adaptive TIN Models[J]. International Archives of the Photogrammetry, Remote Sensing and Spatial Information Sciences, XXXIII(Pt.B4/1):110-117.
## lidR brings two methods, the PMF (Progressive Morphological Filter) and the CSF (Cloth Simulation Filter), if you access each function you will find the papers

# PMF method
ws = seq(3, 12, 3)
#ws = seq(3, 21, 3) # after you test with the ws 3-12 test swapping this parameter 3-21
th = seq(0.1, 1.5, length.out = length(ws))
las = classify_ground(las, algorithm = pmf(ws,th))
plot(las, color = "Classification")

# filter and visualize only the ground points
gnd = filter_ground(las)
plot(gnd, size = 3, bg = "white", color = "Classification")


# 3) Creating a DTM --------------------------------------------------------------

# calc using three different methods
dtm = grid_terrain(las, res = 1, algorithm = tin())

# visualize
plot(dtm, main="DTM")
plot_dtm3d(dtm)



# 4) Height Normalization --------------------------------------------------------------

# normalize the point cloud by subtracting the continuous-DTM (TIN)
nlas = normalize_height(las, tin())
plot(nlas)
plot(las)
## Compare the nlas with the las, what is the difference?



# 5) Filtering ---------------------------------------------------------------

# for filtering experiment: let's add 20 artificial outliers
set.seed(314)
id = round(runif(20, 0, npoints(nlas)))
set.seed(42)
err = runif(20, -50, 50)
nlas$Z[id] = nlas$Z[id] + err
plot(nlas)

# now we can filter for points with very high height, which are not trees
nlas = filter_poi(nlas, Z > 0 & Z < 100)
plot(nlas)

# noise filtering
nlas_noise = classify_noise(nlas, ivf(1,5))
nlas_noise = filter_poi(nlas_noise, Classification != LASNOISE)
plot(nlas_noise)
nlas = nlas_noise



# 6) Creating a CHM ----------------------------------------------------------
## Remember to use nlas and not las here, because we normalized it ok!

# Pit free method
## Series of sequential height thresholds where Delaunay triangulations are applied to first returns
chm = grid_canopy(nlas, res = 1, pitfree(subcircle = 0.15, thresholds = c(0, 10, 20), max_edge = c(8, 1.5)))
plot(chm, col = hgt_col)
writeRaster(chm, filename = "chm_pitfree.tif", overwrite=T)




# 7) Extracting LiDAR Metrics -----------------------------------------------------------------

# Simple metrics
metrics_cl = cloud_metrics(nlas, func = ~mean(Z))
metrics_cl
metrics_grid = grid_metrics(nlas, func = ~mean(Z), res = 10)
plot(metrics_grid)

# more complex metrics
metrics_cl = cloud_metrics(nlas, func = .stdmetrics)
metrics_cl
metrics_grid = grid_metrics(nlas, func = .stdmetrics, res = 10)
metrics_grid
## we have 56 metrics !
plot(metrics_grid[[1:16]], col = hgt_col)
plot(metrics_grid[[17:32]], col = hgt_col)
plot(metrics_grid[[33:48]], col = hgt_col)
plot(metrics_grid[[39:56]], col = hgt_col)




# Summary - Read this only after you are done with each step --------------


## Ok, after all this work, let's summarize what we've learned:

# 1) Point cloud Reading and Plotting: you can see the data frames and headers inside the LAS files
# with a few commands, the display tools are very easy to use
#
# 2) Ground classification: there are many methods to do this. There is no best method, but one that work well for your data.
#
# 3) Creating a DTM: it is easy to generate a DTM surface, there are more than 3 methods to do this,
# the simplest being the TIN method, TIN method is based on the Delaunay Triangulation, it needs NO parameters (this is very good),
# however it cannot extrapolate outside of the bounding box and its weaker at the edges
#
# 4) Height Normalization: forest height is now normalized and you can see the differences in tree height
#
# 5) Filtering: it is very important to test and use noise filtering to remove spurious points
#
# 6) Creating a CHM: which CHM is better depends on the purpose. P2R with disc is simple but leaves lots of empty cells,
# triangulation smooths a lot, Pitfree is less smoothed than delaunay and fills lots of holes from the P2R with disc
#
# 7) LiDAR Metrics: lots of metrics, their usefulness depend on the application, e.g. they can be very useful for biomass estimates

## This is just the tip of the iceberg ! There is a lot more to learn about LiDAR data processing.
## You can find great step-by-step tutorials in the lidR Book (https://jean-romain.github.io/lidRbook) for lots of lidar processing.

# Ricardo Dalagnol
# ricds@hotmail.com





#
# The End -----------------------------------------------------------------
#































# are you hungry for more lidar processing ? (:



















































# Secret Bonus code on Ground Classification and DTM calculation


## Compare two ground classification algorithms
# PMF method
ws = seq(3, 12, 3)
#ws = seq(3, 21, 3) # after you test with the ws 3-12 test swapping this parameter 3-21
th = seq(0.1, 1.5, length.out = length(ws))
las_ground_pmf = classify_ground(las, algorithm = pmf(ws,th))
plot(las_ground_pmf, color = "Classification")

# CSF method
mycsf = csf(TRUE, 1, 1, time_step = 1)
las_ground_csf = classify_ground(las, mycsf)
plot(las_ground_csf, color = "Classification")




## Compare different DTM interpolation

# calc using three different methods
dtm1 = grid_terrain(las_ground_pmf, res = 1, algorithm = knnidw(k = 6L, p = 2))
dtm2 = grid_terrain(las_ground_pmf, res = 1, algorithm = tin())
dtm3 = grid_terrain(las_ground_pmf, res = 1, algorithm = kriging(k = 10L))

# visualize
plot(dtm1, main="knnidw")
plot(dtm2, main="tin")
plot(dtm3, main="kriging")
plot_dtm3d(dtm1)
plot_dtm3d(dtm2)
plot_dtm3d(dtm3)
## which one looks better?

## now test interpolate the DTM with the two different ground filters
dtm_pmf = grid_terrain(las_ground_pmf, res = 1, algorithm = tin())
dtm_csf = grid_terrain(las_ground_csf, res = 1, algorithm = tin())
plot(dtm_pmf, main="PMF")
plot(dtm_csf, main="CSF")
plot_dtm3d(dtm_pmf)
plot_dtm3d(dtm_csf)
## which one is smoother?



# Secret Bonus codes on CHM creation
## You can play with different ways to calculate the CHM from the Point Cloud

# point2raster, simply assign the highest elevation within a grid (e.g. 1 x 1 m)
chm = grid_canopy(nlas, res = 1, algorithm = p2r())
plot(chm, col = hgt_col)
writeRaster(chm, filename = "chm_p2r_1m.tif", overwrite=T)
mean(chm[], na.rm=T)
## Drawback: Some pixels can be empty if the grid resolution is too fine for the available point density

# Now compare the results with a CHM with 0.5 m
## You can also play around with other grid sizes, for example 30-m Landsat resolution :)
chm = grid_canopy(nlas, res = 0.5, algorithm = p2r())
plot(chm, col = hgt_col)
mean(chm[], na.rm=T)
writeRaster(chm, filename = "chm_p2r.tif", overwrite=T)

# additional processing here is to assume points are instead disks of a know radius (e.g. 15 cm). This is because the laser footprint is not a point, but rather a circular area
chm = grid_canopy(nlas, res = 0.5, algorithm = p2r(subcircle = 0.15))
plot(chm, col = hgt_col)
writeRaster(chm, filename = "chm_p2r_disc.tif", overwrite=T)
mean(chm[], na.rm=T)
## Is it better than the previous one?

# additional argument to this method, na filling with the TIN method, takes a while but no pixel is left empty
chm = grid_canopy(nlas, res = 0.5, algorithm = p2r(subcircle = 0.15, na.fill = tin()))
plot(chm, col = hgt_col)
writeRaster(chm, filename = "chm_p2r_disc_fill.tif", overwrite=T)

# Delanauy Triangulation
## Uses the Delaunay method to interpolate the heights.
## The max edge parameter removes large triangles where you don't have data
chm = grid_canopy(nlas, res = 0.5, algorithm = dsmtin(max_edge = 8))
plot(chm, col = hgt_col)
writeRaster(chm, filename = "chm_triangulation.tif", overwrite=T)
## do a quick test removing the max_edge parameter, just use dsmtin(), what is the difference in result?




## Experience for you: Load all the CHM files in the QGIS and compare them visually. What are the differences of each one? Which one looks 'better'?







#░░░░░░░░▄▄▄▀▀▀▄▄███▄░░░░░░░░░░░░░░
#░░░░░▄▀▀░░░░░░░▐░▀██▌░░░░░░░░░░░░░
#░░░▄▀░░░░▄▄███░▌▀▀░▀█░░░░░░░░░░░░░
#░░▄█░░▄▀▀▒▒▒▒▒▄▐░░░░█▌░░░░░░░░░░░░
#░▐█▀▄▀▄▄▄▄▀▀▀▀▌░░░░░▐█▄░░░░░░░░░░░
#░▌▄▄▀▀░░░░░░░░▌░░░░▄███████▄░░░░░░
#░░░░░░░░░░░░░▐░░░░▐███████████▄░░░
#░░░░░le░░░░░░░▐░░░░▐█████████████▄
#░░░░toucan░░░░░░▀▄░░░▐█████████████▄ 
#░░░░░░has░░░░░░░░▀▄▄███████████████ 
#░░░░░arrived░░░░░░░░░░░░█▀██████░░

