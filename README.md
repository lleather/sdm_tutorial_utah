---
title: "Manipulating and Using Spatial Data in R"
author: Lila Leatherman
date: 04/17/2019
output: html_notebook
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

```{r}
# LOAD LIBRARIES
  #install.packages("here)
    library(here)
  #install.packages("tidyverse")
    library(tidyverse)
  #install.packages("rgeos")  
    library(rgeos)
  #install.packages("rgbif")
    library(rgbif)
  #install.packages("maps")
    library(maps)
  #install.packages("maptools")
    library(maptools)  
  #install.packages("raster")
    library(raster)
  #install.packages("rgdal") 
    library(rgdal) 
  #install.packages("biomod2")
    library(biomod2)
  #install.packages("ggplot2")
    library(ggplot2)
  
#set chunk options for writing
knitr::opts_chunk$set(message = FALSE, warning = FALSE)


```

## Load and organize spatial data 

First, we're going to load some background spatial data for our project. 

First, we'll get a polygon representing the state of Utah. There are several different packages that can help you acquire administrative boundaries, but here we're using getData() from the raster package.

```{r}

#a way to get state data
states_list <- c('Utah')
states_all <- getData("GADM",country="USA",level=1, path = here("data/UT"))
UT.shp <- states_all[states_all$NAME_1 %in% states_list,]

#inspect
#UT.shp
plot(UT.shp)


```

This loads and plots a shapefile for Utah. If you downloaded the whole directory, this file has already been saved and provided for you, but you can save it again for practice. 

```{r , results = FALSE }
#create directory for output
dir.create(here("data/UT"))

#export
writeOGR(UT.shp, here("data/UT"), layer = "UT", driver = "ESRI Shapefile", overwrite = TRUE)

#load back in 
UT.shp <- shapefile(here("data/UT/UT.shp"))

```

Next, we're going to load in our climate data. These data are commonly used for species distribution modeling and represent environmental variables that combine different climatic variables into variables that are more meaningful for species. These data were downloaded from http://worldclim.org/version2 at 30s resolution and cropped to Utah. These data are in Raster stack format - a collection of raster layers. (What is this analagous to in Arc?)

Here, I'm using the  stack(), commented out below, can also be used to read in a raster stack, a .grd file. Alternately, the command readRDS can also be used to load the data that I have saved in .Rdata format, which is specific to R. 

```{r}
envStack <- stack(here("data/climate/envStack_init.grd"))
#envStack <- readRDS(here("data/climate/envStack_init.RData"))

#inspect
plot(envStack)

```

We have environmental data, but we need to crop them only to the spatial extent that we're interested in: in this case, Utah. 

```{r}
#crop to Utah extent
envStack_UT <- crop(envStack, UT.shp)

#inspect - this just crops to the spatial extent of the object, not the outline of the polygon
plot(envStack_UT)

# crop to Utah boundaries
envStack_UT <- mask(envStack, UT.shp)
plot(envStack_UT)

#export
writeRaster(envStack_UT, file = here("data/climate/envStack_UT.grd"), options = "INTERLEAVE=BAND", overwrite=TRUE)
saveRDS(envStack_UT, file = here("data/climate/envStack_UT.RData"))

```

Looks good! We can also export our final raster stack. In this case, we can write the file as a raster, or save the data as a .Rdata file. Be careful which file format you save a raster stack in-- even though some file formats can be used for both single and multiple layer raster data (e.g., .tif), these formats do not preserve the names of the layers in a raster stack. 

### Load occurrence data - from GBIF

Next, we'll load occurrence data for our species of interest from a few different sources. First, we download data from the Global Biodiversity Information Facility (GBIF) for our species of interest. 

```{r}
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

#inspect - we can just use plot() to plot the spatial data
plot(gbif.POTR.shp)
```
And we have points! 

Unfortunately, looks like we have some extraneous points that we didn't expect here, so we'll need to crop them out.

```{r}
#make sure we only have points from UT
#basically, subsetting the points to the ones that fall  have to do this spatial operation on a spatial object
gbif.POTR.shp <- gbif.POTR.shp[UT.shp, ]
plot(gbif.POTR.shp)

#replace file with the subsetted points in our area of interest
#we can access just the data frame portion of the shapefile
gbif.POTR <- gbif.POTR.shp@data %>%
  select(-optional)

#export again so we have the most up-to-date version saved! both for our collaborators, and FUTURE US
write.csv(gbif.POTR, "./data/occurrence/GBIF/gbif.POTR.csv", row.names = FALSE)

```

Let's plot these on a map!

```{r}
#to plot points, run this whole chunk 

maps::map(database = "state", regions = "utah")
    points(gbif.POTR.shp)
```


### Load occurrence data - from the Forest Service FIA

I downloaded and processed these data from the Forest Service website, which you can access here: https://apps.fs.usda.gov/fia/datamart/ . 

```{r}

# Normally, you will just load one version of the data to manipulate. But here, I'm showing you a coupel different ways to load in these data.

#load shapefile
fia.POTR.shp <- shapefile(here("data/occurrence/FIA/fia.POTR.shp"))

#load .csv
fia.POTR <- read.csv(here("data/occurrence/FIA/FIA_POTR_UT_presAbs.csv"))

```

Let's look at these data: 

```{r}
# plot with base R, again needs to be a shapefile to plot like this
# color by presence / absence recorded
plot(UT.shp)
    points(fia.POTR.shp, pch = 21, bg = "white", cex = 0.5)
    points(fia.POTR.shp[fia.POTR$presAbs == 1,], pch = 21, bg = "dodgerblue")
```
The FIA data include points both where POTR was observed, and where it was not observed.

### Make the two data sources play together!

```{r}
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

#export your data!
write.csv(POTR.dat, here("data/occurrence/full_POTR.csv"), row.names = FALSE )

```

### Load forest mask 

Another task you might want to do is to only look at one raster, within the extent of another raster. We won't be using this today, but I have provided this as an example. We have a layer that represents areas of forest in Utah, which was downloaded and prepped from : https://swregap.org/data/landcover/

```{r}
#load mask layer
forest_mask <- raster(here("data/ut_landcover/ut_forestmask.tif"))

#inspect
plot(forest_mask)

#make sure mask and layer to be masked have same extent and projection (and resolution?) - otherwise, it won't work!
compareRaster(forest_mask, envStack_UT)

#perform the operation
envStack_mask <- mask(x = envStack_UT, mask = forest_mask, maskvalue = 0)

#inspect
plot(envStack_mask)

```

So now, if we wanted it-- we only have environmental data for areas where there is forest in Utah. 

### Extract environmental data to points

For species distribution modeling, we need to create a data frame that has the values for environmental variables at each of our points. We can do an extract operation to get this information. We do this using the extract() function in the raster package. You can extract either by specifying the coordinates, or by using the shapefile 

```{r}
## extract enviro data to points 

#by specifying coordinates
env.dat <- raster::extract(x = envStack_UT, y = POTR.dat[,c("lon", "lat")])
head(env.dat)
str(env.dat)

# # by using the shapefile which is already spatial 
# env.dat <- raster::extract(x = envStack_UT, y = POTR.dat.shp)
# head(env.dat)


plot(env.dat[,2] ~ factor(POTR.dat$presAbs))


```

## Prep the model

### Prep the data for the model

Unlike the previous steps we've completed, this step is more exclusive to species distribution modeling. You can run the following chunks of code, which are specific to the species distribution modeling process

```{r}
POTR.mod.dat <- BIOMOD_FormatingData(resp.var = as.numeric(POTR.dat$presAbs),
                                     resp.xy = POTR.dat[, c("lon", "lat")],
                                     #resp.var = POTR.dat.shp, # for input: can use shapefile with the presence-absence response in the @data slot
                                     expl.var = stack(envStack_UT),
                                     #eval.resp.var = ,
                                     #PA.strategy = "random", 
                                     #PA.nb.rep = 0, 
                                     #PA.nb.absences = 0,
                                     resp.name = "Populus.tremuloides")
POTR.mod.dat

BIOMOD_ModelingOptions() # need to install java in order to run Maxent.Phillips; we won't be doing this today because it can be pretty finicky!
#myBiomodOptions <- BIOMOD_ModelingOptions(MAXENT.Phillips = list(path_to_maxent.jar = "maxent/maxent.jar"))
```

### Run the model

```{r}
POTR.mod <- BIOMOD_Modeling(data = POTR.mod.dat, 
                            #models = c('GLM','GAM','ANN','RF','MAXENT.Tsuruoka'),  
                            models = c('GBM','ANN','RF'), 
                            #SaveObj = TRUE,
                            #models.options = myBiomodOptions,
                            # , DataSplit = 80
                            VarImport = 1)

POTR.mod
```

So now we've built a suite of species distribution models! This graph is one method of evaluating the boosted regression tree model, abbreviated in the model formulation as "GBM," or "generalized boosted model." We can plot and visualize these in a couple different ways. 

Below I use tidyverse techniques to manipulate and plot these evaluations in a more attractive way, using the ggplot package. 

## View and evaluate the model

### Get model evaluations

```{r}

POTR.mod.eval <- get_evaluations(POTR.mod)
    dimnames(POTR.mod.eval)
    POTR.mod.eval["TSS","Testing.data",,,]
    POTR.mod.eval["KAPPA","Testing.data",,,]
    POTR.mod.eval["ROC","Testing.data",,,]
    
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

So, this has reorganized the data so that we are now plotting model evaluation statistics for each model type, side by side. 

We can do the same thing to look at the importance of different variables in the model. Here, this means how much each environmental variable contributes to the model.

### Get variable importance 

```{r}
get_variables_importance(POTR.mod)

```

We can also prep and plot this a little more attractively. 
 
```{r}
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

### Create an ensemble model-- averaging all of the models together

Lastly in our species distribution modeling process, we create an ensemble model that averages all the models together. Here, this means we make a united map or raster of our results. First we'll build the model, then we'll project it. 

```{r}

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

myBiomodEM
get_evaluations(myBiomodEM)

```

### Project the model spatially

```{r}

myBiomodProj <- BIOMOD_Projection(modeling.output = POTR.mod,
                    new.env = stack(envStack_UT),
                    proj.name = 'current' ,
                    selected.models = 'all' , # will return separate projections for each model 
                    binary.meth = 'TSS' ,
                    compress = 'xz' ,
                    clamping.mask = F,
                    output.format = '.grd' )

myBiomodProj

plot(myBiomodProj)

```

We can also plot these one at a time.

```{r}
plot(myBiomodProj, str.grep = 'RF' )

```

If we want to plot the map on its own, we can extract the results to a new object.This is the output that you can save and manipulate. 

```{r}

myCurrentProj <- get_predictions(myBiomodProj)

#combine projections from all models
presentResult <- calc(myCurrentProj,fun = median); #Choose whatever descriptive statistic you'd like
plot(presentResult)

myCurrentProj
plot(myCurrentProj)

```

### Save output

```{r}
dir.create(here("model_output"))
dir.create(here("model_output/rasters"))

# save raster stack
writeRaster(myCurrentProj, filename = here("model_output/rasters/biomod_out.grd"), options = "INTERLEAVE=BAND", overwrite = TRUE)

#load and inspect
present_In <- stack(here("model_output/rasters/biomod_out.grd"))
present_In

```

## Make a map using ggplot

I have included a supplement that will go into more detail, but you can also make maps using ggplot2 that I think are more straightforward to customize. 

```{r}
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

