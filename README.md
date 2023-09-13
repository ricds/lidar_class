# lidar_class repository

Repository for class on Airborne LiDAR data processing.

## How to use?

Download the contents into a folder that is easy to access. Or clone repository if you want.

## Repository Contents:

-   **1_lidR_basic_processing.R** : R code with basic explanations on how to load a LAS point cloud file in R, look at data, plot point cloud, filter it, and transform it to a Canopy Height Model or other LiDAR metrics.

-   **2_lidar_estimate_carbon_density_longo.R** : R code to estimate aboveground carbon density based on ALS data.

-   **ORIG_NP_T-0262\_\_638250_8876250.las** : point cloud data from the Brazilian Amazon in LAS format, 250 x 250 meters in the ground. Sample data come from EBA/INPE acquisition.

-   **Decks used in classes in ppt format**

## Pre-requisites

-   [R](https://cran.r-project.org/)

-   [RStudio](https://posit.co/products/open-source/rstudio/)

-   *lidR R package*, but also a few others.\
    Run these lines of code inside R to install required packages:

    ``` r
    neededPackages = c("lidR", "sp", "raster", "rgdal", "rgl", "mapview", "gstat")
    pkgTest = function(x) { if (x %in% rownames(installed.packages()) == FALSE) { install.packages(x, dependencies= TRUE) }; library(x, character.only = TRUE) }
    for (package in neededPackages) { pkgTest(package) }
    ```

## Courses

*This material was used at :*

> LiDAR lecture in Forest course of Remote Sensing graduate program at INPE. {2020, 2021, 2022}. (3 hours)

> LiDAR lecture in GIS graduate program at University of Manchester. {2021, 2022, 2023}. (3 hours)

> LiDAR lecture in Remote Sensing and GIS for Ecology course at UNPESP. 2021. (4 hours)

> LiDAR short course at WorCAP 2023. 2023. (1h30 hours). <https://www.gov.br/inpe/pt-br/eventos/worcap/2023>. \<[ppt link](https://github.com/ricds/lidar_class/raw/main/WORCAP2023_Introduction_Lidar.pptx)\>.

## *Contact*

Ricardo Dal'Agnol da Silva\
University of California Los Angeles (UCLA)\
e-mails: [ricds\@hotmail.com](mailto:ricds@hotmail.com) ; [dalagnol\@ucla.edu](mailto:dalagnol@ucla.edu)\
<https://ricds.wordpress.com/>\
Follow me on Twitter [\@RicardoDalagnol](https://twitter.com/RicardoDalagnol) :)
