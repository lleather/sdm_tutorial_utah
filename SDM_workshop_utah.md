---
title: "Manipulating and Using Spatial Data in R"
author: Lila Leatherman
date: 04/17/2019
output: 
  html_document:
    keep_md: true
  pdf_document: default
---

## Welcome! 

This workshop will show you how to load, manipulate, visualize, and analyze spatial data in R. Our end application will be species distribution models, but you can use similar techniques for other kinds of spatial analysis and modeling.  

Our example today will be species distribution models for aspen (Populus tremuloides) in Utah. See [Elith et al. 2009](https://www.annualreviews.org/doi/abs/10.1146/annurev.ecolsys.110308.120159) for a broader overview on species distribution modeling. 

### A few notes: 

- This workshop was developed in R Markdown. 
As a user of this workshop, this means that if you are working from a PDF version of this report, you are seeing the same content that is in the script. 
As a future user of R Markdown, this means that you can write code, make figures, and annotate your findings all within one document-- that you can then export for multiple uses. Please take the opportunity to explore both the PDF and .Rmd versions of this workshop! 

- This workshop uses the [here](https://github.com/jennybc/here_here) package to create dynamic references to file paths. This means you never have to use setwd() again! By creating an R project in a folder of your choice, here() identifies the source or home directory for the project on your machine, so that you don't have to manually set the location of the files or scripts. 

- If you download or clone the entire github repository, you will have all the data you need to run this model. The repository also already includes prepared data, but this workflow takes you through all the steps you would use to load, prep, and use the data on your own. 

## Learning Objectives 

Over the course of this workshop, you will learn how to: 

- read and write point data
- read and write polygon data
- read and write raster data
- combine different types of spatial data
- crop and mask spatial data
- set and change the projection of a dataset
- manipulate data using dplyr
- create simple visualizations using base R and ggplot2

Additionally, we will learn some foundational concepts behind species distribution modeling.

## Set up

First, we need to install all of the packages we need for this project. 


```r
# LOAD LIBRARIES
  #install.packages("here)
    library(here)
```

```
## here() starts at /Users/lilaleatherman/Documents/github/sdm_tutorial_utah
```

```r
  #install.packages("tidyverse")
    library(tidyverse)
```

```
## ── Attaching packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──
```

```
## ✔ ggplot2 3.1.0       ✔ purrr   0.3.2  
## ✔ tibble  2.1.1       ✔ dplyr   0.8.0.1
## ✔ tidyr   0.8.3       ✔ stringr 1.4.0  
## ✔ readr   1.1.1       ✔ forcats 0.3.0
```

```
## Warning: package 'tibble' was built under R version 3.5.2
```

```
## Warning: package 'tidyr' was built under R version 3.5.2
```

```
## Warning: package 'purrr' was built under R version 3.5.2
```

```
## Warning: package 'dplyr' was built under R version 3.5.2
```

```
## Warning: package 'stringr' was built under R version 3.5.2
```

```
## ── Conflicts ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
  #install.packages("rgeos")  
    library(rgeos)
```

```
## rgeos version: 0.4-1, (SVN revision 579)
##  GEOS runtime version: 3.6.1-CAPI-1.10.1 
##  Linking to sp version: 1.3-1 
##  Polygon checking: TRUE
```

```r
  #install.packages("rgbif")
    library(rgbif)
  #install.packages("maps")
    library(maps)
```

```
## 
## Attaching package: 'maps'
```

```
## The following object is masked from 'package:purrr':
## 
##     map
```

```r
  #install.packages("maptools")
    library(maptools)  
```

```
## Loading required package: sp
```

```
## Checking rgeos availability: TRUE
```

```r
  #install.packages("raster")
    library(raster)
```

```
## 
## Attaching package: 'raster'
```

```
## The following object is masked from 'package:dplyr':
## 
##     select
```

```
## The following object is masked from 'package:tidyr':
## 
##     extract
```

```r
  #install.packages("rgdal") 
    library(rgdal) 
```

```
## rgdal: version: 1.3-6, (SVN revision 773)
##  Geospatial Data Abstraction Library extensions to R successfully loaded
##  Loaded GDAL runtime: GDAL 2.1.3, released 2017/20/01
##  Path to GDAL shared files: /Library/Frameworks/R.framework/Versions/3.5/Resources/library/rgdal/gdal
##  GDAL binary built with GEOS: FALSE 
##  Loaded PROJ.4 runtime: Rel. 4.9.3, 15 August 2016, [PJ_VERSION: 493]
##  Path to PROJ.4 shared files: /Library/Frameworks/R.framework/Versions/3.5/Resources/library/rgdal/proj
##  Linking to sp version: 1.3-1
```

```r
  #install.packages("biomod2")
    library(biomod2)
```

```
## biomod2 3.3-18 loaded.
## 
## Type browseVignettes(package='biomod2') to access directly biomod2 vignettes.
```

```r
  #install.packages("ggplot2")
    library(ggplot2)
  
#set chunk options for writing
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

#reassign functions that are masked
extract <- raster::extract
select <- dplyr::select
```

## Load and organize spatial data 

First, we're going to load some background spatial data for our project. 

First, we'll get a polygon representing the state of Utah. There are several different packages that can help you acquire administrative boundaries, but here we're using getData() from the raster package.


```r
#a way to get state data
states_list <- c('Utah')
states_all <- getData("GADM",country="USA",level=1, path = here("data/UT"))
UT.shp <- states_all[states_all$NAME_1 %in% states_list,]

#inspect
#UT.shp
plot(UT.shp)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

This loads and plots a shapefile for Utah. If you downloaded the whole directory, this file has already been saved and provided for you, but you can save it again for practice. 


```r
#create directory for output
dir.create(here("data/UT"))

#export
writeOGR(UT.shp, here("data/UT"), layer = "UT", driver = "ESRI Shapefile", overwrite = TRUE)

#load back in 
UT.shp <- shapefile(here("data/UT/UT.shp"))
```

Next, we're going to load in our climate data. These data are commonly used for species distribution modeling and represent environmental variables that combine different climatic variables into variables that are more meaningful for species. These data were downloaded from http://worldclim.org/version2 at 30s resolution and cropped to Utah. These data are in Raster stack format - a collection of raster layers. (What is this analagous to in Arc?)

Here, I'm using the  stack(), commented out below, can also be used to read in a raster stack, a .grd file. Alternately, the command readRDS can also be used to load the data that I have saved in .Rdata format, which is specific to R. 


```r
envStack <- stack(here("data/climate/envStack_init.grd"))
#envStack <- readRDS(here("data/climate/envStack_init.RData"))

#inspect
plot(envStack)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

We have environmental data, but we need to crop them only to the spatial extent that we're interested in: in this case, Utah. 


```r
#crop to Utah extent
envStack_UT <- crop(envStack, UT.shp)

#inspect - this just crops to the spatial extent of the object, not the outline of the polygon
plot(envStack_UT)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
# crop to Utah boundaries
envStack_UT <- mask(envStack, UT.shp)
plot(envStack_UT)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

```r
#export
writeRaster(envStack_UT, file = here("data/climate/envStack_UT.grd"), options = "INTERLEAVE=BAND", overwrite=TRUE)
saveRDS(envStack_UT, file = here("data/climate/envStack_UT.RData"))
```

Looks good! We can also export our final raster stack. In this case, we can write the file as a raster, or save the data as a .Rdata file. Be careful which file format you save a raster stack in-- even though some file formats can be used for both single and multiple layer raster data (e.g., .tif), these formats do not preserve the names of the layers in a raster stack. 

### Load occurrence data - from GBIF

Next, we'll load occurrence data for our species of interest from a few different sources. First, we download data from the Global Biodiversity Information Facility (GBIF) for our species of interest. 


```r
#function also doesn't work if the gbif website is down
gbif.POTR <- occ_search(scientificName = "Populus tremuloides", 
                           return = "data", 
                           hasCoordinate = TRUE, 
                           hasGeospatialIssue = FALSE, 
                           limit = 200000, 
                           country = "US", stateProvince = c("Utah"), 
                           fields = c("name", "decimalLongitude", "decimalLatitude"))

colnames(gbif.POTR) <- c("name", "lon", "lat")

#export
dir.create(path = here("/data/occurrence/GBIF"))
write.csv(gbif.POTR, "./data/occurrence/GBIF/gbif.POTR.csv", row.names = FALSE)

#read back in 
gbif.POTR <- read.csv("./data/occurrence/GBIF/gbif.POTR.csv", stringsAsFactors = FALSE)

###convert to shapefile 

#create object to turn into a shapefile
gbif.POTR.shp <- gbif.POTR 

#set the fields that represent the coordinates
coordinates(gbif.POTR.shp) <- ~ lon + lat

#convert to spatial points daa frame
gbif.POTR.shp <- SpatialPointsDataFrame(coords = coordinates(gbif.POTR.shp), 
                                       data = data.frame(gbif.POTR.shp))

# Same as "define projection" in Arc
proj4string(gbif.POTR.shp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") 

#check that CRS is the same as enviro data
identicalCRS(gbif.POTR.shp, envStack)
```

```
## [1] TRUE
```

```r
#inspect - we can just use plot() to plot the spatial data
plot(gbif.POTR.shp)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-6-1.png)<!-- -->
And we have points! 

Unfortunately, looks like we have some extraneous points that we didn't expect here, so we'll need to crop them out.


```r
#make sure we only have points from UT
#basically, subsetting the points to the ones that fall  have to do this spatial operation on a spatial object
gbif.POTR.shp <- gbif.POTR.shp[UT.shp, ]
plot(gbif.POTR.shp)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

```r
#replace file with the subsetted points in our area of interest
#we can access just the data frame portion of the shapefile
gbif.POTR <- gbif.POTR.shp@data %>%
  dplyr::select(-optional)

#export again so we have the most up-to-date version saved! both for our collaborators, and FUTURE US
write.csv(gbif.POTR, "./data/occurrence/GBIF/gbif.POTR.csv", row.names = FALSE)
```

Let's plot these on a map!


```r
#to plot points, run this whole chunk 

maps::map(database = "state", regions = "utah")
    points(gbif.POTR.shp)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


### Load occurrence data - from the Forest Service FIA

I downloaded and processed these data from the Forest Service website, which you can access here: https://apps.fs.usda.gov/fia/datamart/ . 


```r
# Normally, you will just load one version of the data to manipulate. But here, I'm showing you a coupel different ways to load in these data.

#load shapefile
fia.POTR.shp <- shapefile(here("data/occurrence/FIA/fia.POTR.shp"))

#load .csv
fia.POTR <- read.csv(here("data/occurrence/FIA/FIA_POTR_UT_presAbs.csv"))
```

Let's look at these data: 


```r
# plot with base R, again needs to be a shapefile to plot like this
# color by presence / absence recorded
plot(UT.shp)
    points(fia.POTR.shp, pch = 21, bg = "white", cex = 0.5)
    points(fia.POTR.shp[fia.POTR$presAbs == 1,], pch = 21, bg = "dodgerblue")
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-10-1.png)<!-- -->
The FIA data include points both where POTR was observed, and where it was not observed.

### Make the two data sources play together!


```r
# join the data

# Make the fields align!
#The GBIF data only represent places where POTR was observed, but we need to add a field that indicates this.
#Also, the FIA data do not have a "name" column, so we remove it here
gbif.POTR <- gbif.POTR %>%
  mutate(presAbs = 1) %>%
  select(-name)

# bind the data frames together 
POTR.dat <- bind_rows(fia.POTR, gbif.POTR) %>%
  select(lon, lat, presAbs) # make sure we get only the fields present in both datasets 

# #### alternately, you can join two shapefiles like so: 
# # still have to create presAbs field for the gbif data
# gbif.POTR.shp@data <- data.frame(gbif.POTR.shp@data[c("lon", "lat")],
#                                  presAbs = rep(1, nrow(gbif.POTR.shp)))
# 
# #check that the projections are the same for these data sets
# identicalCRS(gbif.POTR.shp, fia.POTR.shp)
# 
# #bind
# POTR.dat.shp <- rbind(gbif.POTR.shp, fia.POTR.shp)

# plot all together
# ggplot can be used to plot 

# another way to get spatial data - ready for plotting in ggplot. 
# this is just a list of points that can be rendered as a polygon by ggplot. 
# I like ggplot, which is part of the "tidyverse", because I think it's a little easier to understand than base R!
UT <- map_data("state", region = "utah")

ggplot() +
      geom_polygon(data = UT, aes(x=long, y = lat, group = group), fill = NA, color = "black") + 
      geom_point(data = POTR.dat %>% arrange(presAbs), 
                 aes(x = lon, y = lat, color = factor(presAbs))) +
      coord_fixed(1.3) + 
      labs(title = "FIA and GBIF PRESENCE/ABSENCE DATA - Populus tremuloides")
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
#export your data!
write.csv(POTR.dat, here("data/occurrence/full_POTR.csv"), row.names = FALSE )
```

### Load forest mask 

Another task you might want to do is to only look at one raster, within the extent of another raster. We won't be using this today, but I have provided this as an example. We have a layer that represents areas of forest in Utah, which was downloaded and prepped from : https://swregap.org/data/landcover/


```r
#load mask layer
forest_mask <- raster(here("data/ut_landcover/ut_forestmask.tif"))

#inspect
plot(forest_mask)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
#make sure mask and layer to be masked have same extent and projection (and resolution?) - otherwise, it won't work!
compareRaster(forest_mask, envStack_UT)
```

```
## [1] TRUE
```

```r
#perform the operation
envStack_mask <- mask(x = envStack_UT, mask = forest_mask, maskvalue = 0)

#inspect
plot(envStack_mask)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

So now, if we wanted it-- we only have environmental data for areas where there is forest in Utah. 

### Extract environmental data to points

For species distribution modeling, we need to create a data frame that has the values for environmental variables at each of our points. We can do an extract operation to get this information. We do this using the extract() function in the raster package. You can extract either by specifying the coordinates, or by using the shapefile 


```r
## extract enviro data to points 

#by specifying coordinates
env.dat <- raster::extract(x = envStack_UT, y = POTR.dat[,c("lon", "lat")])
head(env.dat)
```

```
##      Annual.Mean.Temperature Mean.Diurnal.Range Isothermality
## [1,]                9.487500           14.79167      37.44726
## [2,]               10.350000           11.41667      31.62512
## [3,]                9.054167           14.00833      34.33415
## [4,]                8.079166           13.24167      36.17942
## [5,]                6.712500           16.45833      43.08464
## [6,]                9.108334           15.66667      34.89235
##      Temperature.Seasonality Max.Temperature.of.Warmest.Month
## [1,]                887.3636                             28.7
## [2,]                904.1747                             28.3
## [3,]                952.7078                             28.7
## [4,]                849.0461                             26.4
## [5,]                795.2562                             25.6
## [6,]               1031.5254                             30.2
##      Min.Temperature.of.Coldest.Month Temperature.Annual.Range
## [1,]                            -10.8                     39.5
## [2,]                             -7.8                     36.1
## [3,]                            -12.1                     40.8
## [4,]                            -10.2                     36.6
## [5,]                            -12.6                     38.2
## [6,]                            -14.7                     44.9
##      Mean.Temperature.of.Wettest.Quarter
## [1,]                            8.233333
## [2,]                           21.183334
## [3,]                           20.133333
## [4,]                           18.383333
## [5,]                           16.283333
## [6,]                           16.616667
##      Mean.Temperature.of.Driest.Quarter
## [1,]                          20.750000
## [2,]                          14.083333
## [3,]                          -1.633333
## [4,]                          11.383333
## [5,]                           9.916667
## [6,]                          -3.950000
##      Mean.Temperature.of.Warmest.Quarter
## [1,]                            20.75000
## [2,]                            21.90000
## [3,]                            20.96667
## [4,]                            18.88333
## [5,]                            16.86667
## [6,]                            21.50000
##      Mean.Temperature.of.Coldest.Quarter Annual.Precipitation
## [1,]                          -1.1666666                  309
## [2,]                          -0.2500001                  287
## [3,]                          -2.6333334                  212
## [4,]                          -1.8500000                  397
## [5,]                          -2.6666667                  316
## [6,]                          -3.9499998                  273
##      Precipitation.of.Wettest.Month Precipitation.of.Driest.Month
## [1,]                             34                            16
## [2,]                             35                            11
## [3,]                             27                            12
## [4,]                             45                            17
## [5,]                             49                            16
## [6,]                             30                            16
##      Precipitation.Seasonality Precipitation.of.Wettest.Quarter
## [1,]                  19.20286                               98
## [2,]                  28.65912                               97
## [3,]                  32.08870                               79
## [4,]                  25.09788                              129
## [5,]                  34.61431                              115
## [6,]                  21.77701                               85
##      Precipitation.of.Driest.Quarter Precipitation.of.Warmest.Quarter
## [1,]                              63                               63
## [2,]                              44                               78
## [3,]                              39                               67
## [4,]                              72                              105
## [5,]                              53                               96
## [6,]                              49                               71
##      Precipitation.of.Coldest.Quarter
## [1,]                               68
## [2,]                               69
## [3,]                               39
## [4,]                               86
## [5,]                               65
## [6,]                               49
```

```r
str(env.dat)
```

```
##  num [1:3861, 1:19] 9.49 10.35 9.05 8.08 6.71 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : NULL
##   ..$ : chr [1:19] "Annual.Mean.Temperature" "Mean.Diurnal.Range" "Isothermality" "Temperature.Seasonality" ...
```

```r
# # by using the shapefile which is already spatial 
# env.dat <- raster::extract(x = envStack_UT, y = POTR.dat.shp)
# head(env.dat)


plot(env.dat[,2] ~ factor(POTR.dat$presAbs))
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

## Prep the model

### Prep the data for the model

Unlike the previous steps we've completed, this step is more exclusive to species distribution modeling. You can run the following chunks of code, which are specific to the species distribution modeling process


```r
POTR.mod.dat <- BIOMOD_FormatingData(resp.var = as.numeric(POTR.dat$presAbs),
                                     resp.xy = POTR.dat[, c("lon", "lat")],
                                     #resp.var = POTR.dat.shp, # for input: can use shapefile with the presence-absence response in the @data slot
                                     expl.var = stack(envStack_UT),
                                     #eval.resp.var = ,
                                     #PA.strategy = "random", 
                                     #PA.nb.rep = 0, 
                                     #PA.nb.absences = 0,
                                     resp.name = "Populus.tremuloides")
```

```
## 
## -=-=-=-=-=-=-=-=-=-= Populus.tremuloides Data Formating -=-=-=-=-=-=-=-=-=-=
## 
## > No pseudo absences selection !
##       ! No data has been set aside for modeling evaluation
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Done -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
```

```r
POTR.mod.dat
```

```
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-= 'BIOMOD.formated.data' -=-=-=-=-=-=-=-=-=-=-=-=-=
## 
## sp.name =  Populus.tremuloides
## 
## 	 530 presences,  3331 true absences and  0 undifined points in dataset
## 
## 
## 	 19 explanatory variables
## 
##  Annual.Mean.Temperature Mean.Diurnal.Range Isothermality  
##  Min.   :-1.133          Min.   : 8.258     Min.   :28.58  
##  1st Qu.: 4.292          1st Qu.:12.950     1st Qu.:36.24  
##  Median : 6.562          Median :14.283     Median :38.01  
##  Mean   : 6.516          Mean   :14.183     Mean   :38.06  
##  3rd Qu.: 8.783          3rd Qu.:15.417     3rd Qu.:39.91  
##  Max.   :15.442          Max.   :19.967     Max.   :45.62  
##  Temperature.Seasonality Max.Temperature.of.Warmest.Month
##  Min.   : 674.5          Min.   :14.50                   
##  1st Qu.: 793.8          1st Qu.:22.00                   
##  Median : 835.2          Median :25.30                   
##  Mean   : 843.5          Mean   :25.13                   
##  3rd Qu.: 889.7          3rd Qu.:28.30                   
##  Max.   :1156.1          Max.   :34.90                   
##  Min.Temperature.of.Coldest.Month Temperature.Annual.Range
##  Min.   :-19.00                   Min.   :27.30           
##  1st Qu.:-13.40                   1st Qu.:34.80           
##  Median :-12.20                   Median :37.60           
##  Mean   :-12.09                   Mean   :37.22           
##  3rd Qu.:-10.80                   3rd Qu.:39.60           
##  Max.   : -3.10                   Max.   :50.70           
##  Mean.Temperature.of.Wettest.Quarter Mean.Temperature.of.Driest.Quarter
##  Min.   :-6.283                      Min.   :-7.550                    
##  1st Qu.: 6.050                      1st Qu.:-1.467                    
##  Median :12.867                      Median : 7.350                    
##  Mean   :11.103                      Mean   : 6.798                    
##  3rd Qu.:16.650                      3rd Qu.:13.917                    
##  Max.   :24.050                      Max.   :23.267                    
##  Mean.Temperature.of.Warmest.Quarter Mean.Temperature.of.Coldest.Quarter
##  Min.   : 8.267                      Min.   :-8.917                     
##  1st Qu.:14.533                      1st Qu.:-5.150                     
##  Median :17.217                      Median :-3.600                     
##  Mean   :17.224                      Mean   :-3.390                     
##  3rd Qu.:19.950                      3rd Qu.:-1.700                     
##  Max.   :26.717                      Max.   : 5.350                     
##  Annual.Precipitation Precipitation.of.Wettest.Month
##  Min.   :155.0        Min.   :21.00                 
##  1st Qu.:312.0        1st Qu.:36.00                 
##  Median :374.0        Median :42.00                 
##  Mean   :372.2        Mean   :42.47                 
##  3rd Qu.:431.0        3rd Qu.:48.00                 
##  Max.   :639.0        Max.   :73.00                 
##  Precipitation.of.Driest.Month Precipitation.Seasonality
##  Min.   : 5.00                 Min.   : 6.004           
##  1st Qu.:14.00                 1st Qu.:14.976           
##  Median :20.00                 Median :20.832           
##  Mean   :20.26                 Mean   :21.254           
##  3rd Qu.:25.00                 3rd Qu.:26.762           
##  Max.   :40.00                 Max.   :40.002           
##  Precipitation.of.Wettest.Quarter Precipitation.of.Driest.Quarter
##  Min.   : 58.0                    Min.   : 23.00                 
##  1st Qu.: 99.0                    1st Qu.: 57.00                 
##  Median :114.0                    Median : 72.00                 
##  Mean   :115.4                    Mean   : 71.33                 
##  3rd Qu.:131.0                    3rd Qu.: 86.00                 
##  Max.   :202.0                    Max.   :124.00                 
##  Precipitation.of.Warmest.Quarter Precipitation.of.Coldest.Quarter
##  Min.   : 41.00                   Min.   : 23.00                  
##  1st Qu.: 79.00                   1st Qu.: 66.00                  
##  Median : 91.00                   Median : 87.00                  
##  Mean   : 94.84                   Mean   : 88.26                  
##  3rd Qu.:107.00                   3rd Qu.:110.00                  
##  Max.   :202.00                   Max.   :177.00                  
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
```

```r
BIOMOD_ModelingOptions() # need to install java in order to run Maxent.Phillips; we won't be doing this today because it can be pretty finicky!
```

```
## 
## -=-=-=-=-=-=-=-=-=-=-=-=  'BIOMOD.Model.Options'  -=-=-=-=-=-=-=-=-=-=-=-=
## 
## 
## GLM = list( type = 'quadratic',
##             interaction.level = 0,
##             myFormula = NULL,
##             test = 'AIC',
##             family = binomial(link = 'logit'),
##             mustart = 0.5,
##             control = glm.control(epsilon = 1e-08, maxit = 50
## , trace = FALSE) ),
## 
## 
## GBM = list( distribution = 'bernoulli',
##             n.trees = 2500,
##             interaction.depth = 7,
##             n.minobsinnode = 5,
##             shrinkage = 0.001,
##             bag.fraction = 0.5,
##             train.fraction = 1,
##             cv.folds = 3,
##             keep.data = FALSE,
##             verbose = FALSE,
##             perf.method = 'cv',
##             n.cores = 1),
## 
## GAM = list( algo = 'GAM_mgcv',
##             type = 's_smoother',
##             k = -1,
##             interaction.level = 0,
##             myFormula = NULL,
##             family = binomial(link = 'logit'),
##             method = 'GCV.Cp',
##             optimizer = c('outer','newton'),
##             select = FALSE,
##             knots = NULL,
##             paraPen = NULL,
##             control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
## , maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
## , rank.tol = 1.49011611938477e-08
## , nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
## , optim = list(factr=1e+07)
## , newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
## , outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE
## , efs.lspmax = 15, efs.tol = 0.1, keepData = FALSE, scale.est = fletcher
## , edge.correct = FALSE) ),
## 
## 
## CTA = list( method = 'class',
##             parms = 'default',
##             cost = NULL,
##             control = list(xval = 5, minbucket = 5, minsplit = 5
## , cp = 0.001, maxdepth = 25) ),
## 
## 
## ANN = list( NbCV = 5,
##             size = NULL,
##             decay = NULL,
##             rang = 0.1,
##             maxit = 200),
## 
## SRE = list( quant = 0.025),
## 
## FDA = list( method = 'mars',
##             add_args = NULL),
## 
## MARS = list( type = 'simple',
##              interaction.level = 0,
##              myFormula = NULL,
##              nk = NULL,
##              penalty = 2,
##              thresh = 0.001,
##              nprune = NULL,
##              pmethod = 'backward'),
## 
## RF = list( do.classif = TRUE,
##            ntree = 500,
##            mtry = 'default',
##            nodesize = 5,
##            maxnodes = NULL),
## 
## MAXENT.Phillips = list( path_to_maxent.jar = '/Users/lilaleatherman/Documents/github/sdm_tutorial_utah',
##                memory_allocated = 512,
##                background_data_dir = 'default',
##                maximumbackground = 'default',
##                maximumiterations = 200,
##                visible = FALSE,
##                linear = TRUE,
##                quadratic = TRUE,
##                product = TRUE,
##                threshold = TRUE,
##                hinge = TRUE,
##                lq2lqptthreshold = 80,
##                l2lqthreshold = 10,
##                hingethreshold = 15,
##                beta_threshold = -1,
##                beta_categorical = -1,
##                beta_lqp = -1,
##                beta_hinge = -1,
##                betamultiplier = 1,
##                defaultprevalence = 0.5),
## 
## MAXENT.Tsuruoka = list( l1_regularizer = 0,
##                         l2_regularizer = 0,
##                         use_sgd = FALSE,
##                         set_heldout = 0,
##                         verbose = FALSE)
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
```

```r
#myBiomodOptions <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = "maxent/maxent.jar"))
```

### Run the model


```r
POTR.mod <- BIOMOD_Modeling(data = POTR.mod.dat, 
                            #models = c('GLM','GAM','ANN','RF','MAXENT.Tsuruoka'),  
                            models = c('GBM','ANN','RF'), 
                            #SaveObj = TRUE,
                            #models.options = myBiomodOptions,
                            # , DataSplit = 80
                            VarImport = 1)
```

```
## 
## 
## Loading required library...
## 
## Checking Models arguments...
## 
## Creating suitable Workdir...
## 
## 	> No weights : all observations will have the same weight
## 
## 
## -=-=-=-=-=-=-=-=-= Populus.tremuloides Modeling Summary -=-=-=-=-=-=-=-=-=
## 
##  19  environmental variables ( Annual.Mean.Temperature Mean.Diurnal.Range Isothermality Temperature.Seasonality Max.Temperature.of.Warmest.Month Min.Temperature.of.Coldest.Month Temperature.Annual.Range Mean.Temperature.of.Wettest.Quarter Mean.Temperature.of.Driest.Quarter Mean.Temperature.of.Warmest.Quarter Mean.Temperature.of.Coldest.Quarter Annual.Precipitation Precipitation.of.Wettest.Month Precipitation.of.Driest.Month Precipitation.Seasonality Precipitation.of.Wettest.Quarter Precipitation.of.Driest.Quarter Precipitation.of.Warmest.Quarter Precipitation.of.Coldest.Quarter )
## Number of evaluation repetitions : 1
## Models selected : GBM ANN RF 
## 
## Total number of model runs : 3 
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## 
## 
## -=-=-=- Run :  Populus.tremuloides_AllData 
## 
## 
## -=-=-=--=-=-=- Populus.tremuloides_AllData_Full 
## 
## Model=Generalised Boosting Regression 
## 	 2500 maximum different trees and  3  Fold Cross-Validation
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

```
## 
## 	Evaluating Model stuff...
## 	Evaluating Predictor Contributions... 
## 
## Model=Artificial Neural Network 
## 	 5 Fold Cross Validation + 3 Repetitions
## 	Model scaling...
## 	Evaluating Model stuff...
## 	Evaluating Predictor Contributions... 
## 
## Model=Breiman and Cutler's random forests for classification and regression
## 	Evaluating Model stuff...
## 	Evaluating Predictor Contributions... 
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Done -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
```

```r
POTR.mod
```

```
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-=-= BIOMOD.models.out -=-=-=-=-=-=-=-=-=-=-=-=-=-=
## 
## Modeling id : 1555528742
## 
## Species modeled : Populus.tremuloides
## 
## Considered variables : Annual.Mean.Temperature Mean.Diurnal.Range 
## Isothermality Temperature.Seasonality Max.Temperature.of.Warmest.Month 
## Min.Temperature.of.Coldest.Month Temperature.Annual.Range 
## Mean.Temperature.of.Wettest.Quarter Mean.Temperature.of.Driest.Quarter 
## Mean.Temperature.of.Warmest.Quarter Mean.Temperature.of.Coldest.Quarter 
## Annual.Precipitation Precipitation.of.Wettest.Month 
## Precipitation.of.Driest.Month Precipitation.Seasonality 
## Precipitation.of.Wettest.Quarter Precipitation.of.Driest.Quarter 
## Precipitation.of.Warmest.Quarter Precipitation.of.Coldest.Quarter
## 
## 
## Computed Models :  Populus.tremuloides_AllData_Full_GBM 
## Populus.tremuloides_AllData_Full_ANN Populus.tremuloides_AllData_Full_RF
## 
## 
## Failed Models :  none
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
```

So now we've built a suite of species distribution models! This graph is one method of evaluating the boosted regression tree model, abbreviated in the model formulation as "GBM," or "generalized boosted model." We can plot and visualize these in a couple different ways. 

Below I use tidyverse techniques to manipulate and plot these evaluations in a more attractive way, using the ggplot package. 

## View and evaluate the model

### Get model evaluations


```r
POTR.mod.eval <- get_evaluations(POTR.mod)
    dimnames(POTR.mod.eval)
```

```
## [[1]]
## [1] "KAPPA" "TSS"   "ROC"  
## 
## [[2]]
## [1] "Testing.data" "Cutoff"       "Sensitivity"  "Specificity" 
## 
## [[3]]
## [1] "GBM" "ANN" "RF" 
## 
## [[4]]
## [1] "Full"
## 
## [[5]]
## Populus.tremuloides_AllData 
##                   "AllData"
```

```r
    POTR.mod.eval["TSS","Testing.data",,,]
```

```
##   GBM   ANN    RF 
## 0.637 0.588 0.869
```

```r
    POTR.mod.eval["KAPPA","Testing.data",,,]
```

```
##   GBM   ANN    RF 
## 0.432 0.338 0.649
```

```r
    POTR.mod.eval["ROC","Testing.data",,,]
```

```
##   GBM   ANN    RF 
## 0.880 0.811 0.947
```

```r
# prep and plot these data a little more attractively
    
POTR.mod.eval <- 
     data.frame(POTR.mod.eval[1:3, 1:4, ,,]) %>%
       rownames_to_column("metric") %>%
       gather(-metric, key = "mod", value = "val") %>%
       mutate(mod = gsub(pattern = "Testing.data", "Testingdata", mod)) %>%
       separate(mod, c("val_type", "model"))

#inspect evaluation statistics
print(
POTR.mod.eval %>%
  filter(val_type == "Testingdata") %>%
  ggplot(aes(x = metric, y = val, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "model evaluation - on testing data",
       x = "model type")
) 
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

So, this has reorganized the data so that we are now plotting model evaluation statistics for each model type, side by side. 

We can do the same thing to look at the importance of different variables in the model. Here, this means how much each environmental variable contributes to the model.

### Get variable importance 


```r
get_variables_importance(POTR.mod)
```

```
## , , Full, AllData
## 
##                                       GBM   ANN    RF
## Annual.Mean.Temperature             0.467 0.341 0.127
## Mean.Diurnal.Range                  0.001 0.007 0.011
## Isothermality                       0.002 0.016 0.010
## Temperature.Seasonality             0.004 0.004 0.012
## Max.Temperature.of.Warmest.Month    0.030 0.040 0.067
## Min.Temperature.of.Coldest.Month    0.001 0.025 0.009
## Temperature.Annual.Range            0.001 0.000 0.021
## Mean.Temperature.of.Wettest.Quarter 0.004 0.000 0.015
## Mean.Temperature.of.Driest.Quarter  0.006 0.000 0.017
## Mean.Temperature.of.Warmest.Quarter 0.008 0.070 0.079
## Mean.Temperature.of.Coldest.Quarter 0.015 0.247 0.053
## Annual.Precipitation                0.012 0.100 0.026
## Precipitation.of.Wettest.Month      0.000 0.031 0.006
## Precipitation.of.Driest.Month       0.002 0.145 0.014
## Precipitation.Seasonality           0.005 0.006 0.011
## Precipitation.of.Wettest.Quarter    0.002 0.021 0.011
## Precipitation.of.Driest.Quarter     0.002 0.000 0.022
## Precipitation.of.Warmest.Quarter    0.008 0.002 0.017
## Precipitation.of.Coldest.Quarter    0.008 0.035 0.014
```

We can also prep and plot this a little more attractively. 
 

```r
# investigate variable importance
# still some parsing to do here to visualize
var.imp <- (get_variables_importance(POTR.mod))
var.imp <- data.frame(var.imp)
# var.imp$var <- names(envStack_aoi)
# var.imp

var.imp <-
  var.imp %>%
  mutate(var = names(envStack_UT))%>%
  gather(-var, key = "mod", value = "val") %>%
  separate(mod, c("model_type", "model_level", "data_amt")) 

#plot var importance
print(
var.imp %>%
  #mutate(val = ifelse(model_type == "RF", val*10, val)) %>% #multiply RF values by 10 to compare better? idk why these are so low
  ggplot(aes(y = val, x = reorder(var, val), fill = val, group = interaction(data_amt, model_type))) +
  geom_bar(stat = "identity", position = "dodge") + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title = "Variable importance, by model",
       x = "variable") + 
  facet_grid(~model_type)
)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

### Create an ensemble model-- averaging all of the models together

Lastly in our species distribution modeling process, we create an ensemble model that averages all the models together. Here, this means we make a united map or raster of our results. First we'll build the model, then we'll project it. 


```r
myBiomodEM <- BIOMOD_EnsembleModeling(
                  modeling.output = POTR.mod,
                  chosen.models = 'all',
                  em.by= 'all',
                  eval.metric = c('TSS'),
                  eval.metric.quality.threshold = c(0.6),
                  prob.mean = T,
                  prob.cv = T,
                  prob.ci = T,
                  prob.ci.alpha = 0.05,
                  prob.median = T,
                  committee.averaging = T,
                  prob.mean.weight = T,
                  prob.mean.weight.decay = 'proportional')
```

```
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-= Build Ensemble Models -=-=-=-=-=-=-=-=-=-=-=-=-=
## 
##    ! all models available will be included in ensemble.modeling
##    > Evaluation & Weighting methods summary :
##       TSS over 0.6
## 
## 
##   > mergedAlgo_mergedRun_mergedData ensemble modeling
##    ! Models projections for whole zonation required...
## 	> Projecting Populus.tremuloides_AllData_Full_GBM ...
## 	> Projecting Populus.tremuloides_AllData_Full_RF ...
## 
##    > Mean of probabilities...
## 			Evaluating Model stuff...
##    > Coef of variation of probabilities...
## 			Evaluating Model stuff...
##    > Confidence Interval...
## 			Evaluating Model stuff...
## 			Evaluating Model stuff...
##    > Median of probabilities...
## 			Evaluating Model stuff...
##    >  Committee averaging...
## 			Evaluating Model stuff...
##    > Probabilities weighting mean...
## 		  original models scores =  0.637 0.869
## 		  final models weights =  0.423 0.577
## 			Evaluating Model stuff...
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Done -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
```

```r
myBiomodEM
```

```
## 
## -=-=-=-=-=-=-=-=-=-=-= 'BIOMOD.EnsembleModeling.out' -=-=-=-=-=-=-=-=-=-=-=
## 
## sp.name : Populus.tremuloides
## 
## expl.var.names : Annual.Mean.Temperature Mean.Diurnal.Range Isothermality 
## Temperature.Seasonality Max.Temperature.of.Warmest.Month 
## Min.Temperature.of.Coldest.Month Temperature.Annual.Range 
## Mean.Temperature.of.Wettest.Quarter Mean.Temperature.of.Driest.Quarter 
## Mean.Temperature.of.Warmest.Quarter Mean.Temperature.of.Coldest.Quarter 
## Annual.Precipitation Precipitation.of.Wettest.Month 
## Precipitation.of.Driest.Month Precipitation.Seasonality 
## Precipitation.of.Wettest.Quarter Precipitation.of.Driest.Quarter 
## Precipitation.of.Warmest.Quarter Precipitation.of.Coldest.Quarter
## 
## 
## models computed: 
## Populus.tremuloides_EMmeanByTSS_mergedAlgo_mergedRun_mergedData, Populus.tremuloides_EMcvByTSS_mergedAlgo_mergedRun_mergedData, Populus.tremuloides_EMciInfByTSS_mergedAlgo_mergedRun_mergedData, Populus.tremuloides_EMciSupByTSS_mergedAlgo_mergedRun_mergedData, Populus.tremuloides_EMmedianByTSS_mergedAlgo_mergedRun_mergedData, Populus.tremuloides_EMcaByTSS_mergedAlgo_mergedRun_mergedData, Populus.tremuloides_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
```

```r
get_evaluations(myBiomodEM)
```

```
## $Populus.tremuloides_EMmeanByTSS_mergedAlgo_mergedRun_mergedData
##       Testing.data Cutoff Sensitivity Specificity
## KAPPA        0.603    285      88.491      87.932
## TSS          0.789    246      93.208      85.590
## ROC          0.937    183      98.868      80.036
## 
## $Populus.tremuloides_EMcvByTSS_mergedAlgo_mergedRun_mergedData
##       Testing.data Cutoff Sensitivity Specificity
## KAPPA           NA     NA          NA          NA
## TSS             NA     NA          NA          NA
## ROC             NA     NA          NA          NA
## 
## $Populus.tremuloides_EMciInfByTSS_mergedAlgo_mergedRun_mergedData
##       Testing.data Cutoff Sensitivity Specificity
## KAPPA        0.493  105.0      65.660      89.853
## TSS          0.581   22.0      74.340      83.458
## ROC          0.810   25.5      73.962      84.239
## 
## $Populus.tremuloides_EMciSupByTSS_mergedAlgo_mergedRun_mergedData
##       Testing.data Cutoff Sensitivity Specificity
## KAPPA        0.533    501      86.981      84.839
## TSS          0.795    423      99.623      79.736
## ROC          0.927    424      99.623      79.856
## 
## $Populus.tremuloides_EMmedianByTSS_mergedAlgo_mergedRun_mergedData
##       Testing.data Cutoff Sensitivity Specificity
## KAPPA        0.603    285      88.491      87.932
## TSS          0.789    246      93.208      85.590
## ROC          0.937    183      98.868      80.036
## 
## $Populus.tremuloides_EMcaByTSS_mergedAlgo_mergedRun_mergedData
##       Testing.data Cutoff Sensitivity Specificity
## KAPPA        0.616    747      92.075      87.631
## TSS          0.797    747      92.075      87.631
## ROC          0.927    750      92.075      87.631
## 
## $Populus.tremuloides_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData
##       Testing.data Cutoff Sensitivity Specificity
## KAPPA        0.616  284.0      90.943      87.932
## TSS          0.815  212.0      97.736      83.729
## ROC          0.940  227.5      96.604      85.080
```

### Project the model spatially


```r
myBiomodProj <- BIOMOD_Projection(modeling.output = POTR.mod,
                    new.env = stack(envStack_UT),
                    proj.name = 'current' ,
                    selected.models = 'all' , # will return separate projections for each model 
                    binary.meth = 'TSS' ,
                    compress = 'xz' ,
                    clamping.mask = F,
                    output.format = '.grd' )
```

```
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-= Do Models Projections -=-=-=-=-=-=-=-=-=-=-=-=-=
## 
## 	> Building clamping mask
## 
## 	> Projecting Populus.tremuloides_AllData_Full_GBM ...
## 	> Projecting Populus.tremuloides_AllData_Full_ANN ...
## 	> Projecting Populus.tremuloides_AllData_Full_RF ...
## 
## 	> Building TSS binaries
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Done -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
```

```r
myBiomodProj
```

```
## 
## -=-=-=-=-=-=-=-=-=-=-=-= 'BIOMOD.projection.out' -=-=-=-=-=-=-=-=-=-=-=-=
## 
## Projection directory : Populus.tremuloides/current
## 
## 
## sp.name : Populus.tremuloides
## 
## expl.var.names : Annual.Mean.Temperature Mean.Diurnal.Range Isothermality 
## Temperature.Seasonality Max.Temperature.of.Warmest.Month 
## Min.Temperature.of.Coldest.Month Temperature.Annual.Range 
## Mean.Temperature.of.Wettest.Quarter Mean.Temperature.of.Driest.Quarter 
## Mean.Temperature.of.Warmest.Quarter Mean.Temperature.of.Coldest.Quarter 
## Annual.Precipitation Precipitation.of.Wettest.Month 
## Precipitation.of.Driest.Month Precipitation.Seasonality 
## Precipitation.of.Wettest.Quarter Precipitation.of.Driest.Quarter 
## Precipitation.of.Warmest.Quarter Precipitation.of.Coldest.Quarter
## 
## 
## modeling id : 1555528742 ( 
## Populus.tremuloides/Populus.tremuloides.1555528742.models.out )
## 
## models projected : 
## Populus.tremuloides_AllData_Full_GBM, Populus.tremuloides_AllData_Full_ANN, Populus.tremuloides_AllData_Full_RF
## 
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
```

```r
plot(myBiomodProj)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

We can also plot these one at a time.


```r
plot(myBiomodProj, str.grep = 'RF' )
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

If we want to plot the map on its own, we can extract the results to a new object.This is the output that you can save and manipulate. 


```r
myCurrentProj <- get_predictions(myBiomodProj)

#combine projections from all models
presentResult <- calc(myCurrentProj,fun = median); #Choose whatever descriptive statistic you'd like
plot(presentResult)
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
#myCurrentProj
plot(myCurrentProj[[1]], main = "Aspen (Populus Tremuloides) in Utah")
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-22-2.png)<!-- -->

### Save output


```r
dir.create(here("model_output"))
dir.create(here("model_output/rasters"))

# save raster stack
writeRaster(myCurrentProj, filename = here("model_output/rasters/biomod_out.grd"), options = "INTERLEAVE=BAND", overwrite = TRUE)

#load and inspect
present_In <- stack(here("model_output/rasters/biomod_out.grd"))
present_In
```

```
## class       : RasterStack 
## dimensions  : 601, 601, 361201, 3  (nrow, ncol, ncell, nlayers)
## resolution  : 0.008333333, 0.008333333  (x, y)
## extent      : -114.05, -109.0417, 37, 42.00833  (xmin, xmax, ymin, ymax)
## coord. ref. : +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 
## names       : Populus.tremuloides_AllData_Full_GBM, Populus.tremuloides_AllData_Full_ANN, Populus.tremuloides_AllData_Full_RF 
## min values  :                                   30,                                   19,                                   0 
## max values  :                                  617,                                  353,                                 938
```

## Make a map using ggplot

I have included a supplement that will go into more detail, but you can also make maps using ggplot2 that I think are more straightforward to customize. 


```r
library(rasterVis)
library(viridis)

#choose one layer
RF_model <- present_In$Populus.tremuloides_AllData_Full_RF

p = gplot(RF_model) +
  geom_tile(aes(fill = value)) + 
  scale_fill_viridis(na.value = "white") + 
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  coord_fixed(ratio = 1.3) # sets the xy resolution to a constant value

p
```

![](SDM_workshop_utah_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

