### LOAD PACKAGES -----------------------------------------------------------------------------------------
library(sf)
library(ggmap)
#library(rgdal)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(raster)
library(terra)
library(lubridate)
library(akgfmaps)
library(biomod2)
library(gam)
library(purrr)
library(effsize)
library(dismo)
library(gbm) 
library(pROC)
library(ggrepel)
library(geosphere)
library(usdm)
library(statip)
library(ggpattern)


### SET SPATIAL DETAILS ---------------------------------------------------------
crs.latlon <- "epsg:4326" #lat lon crs

map.crs <- "EPSG:3338"

in.crs = "+proj=longlat +datum=NAD83"

# LOAD SPATIAL LAYERS -----------------------------------------------------------
region_layers <- akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto")

survey_gdb <- "./Data/SAP_layers.gdb"

# Bristol Bay strata multipolygon
BB_strata <- sf::st_read(survey_gdb, layer = "BristolBaySurveyStrata") %>%
              # can also use sf::st_read() to read layers in as sf objects
              vect() %>%
            project(., map.crs)


st_read("./Data/Closure areas/RKCSA_sub.shp") %>%
  vect() -> RKCSA_sub

st_read("./Data/Closure areas/RKCSA.shp") %>%
  vect() -> RKCSA

st_read("./Data/Closure areas/area512.shp") %>%
  vect() -> area512
st_read("./Data/Closure areas/fullBLZ.shp") %>%
  vect() -> fullBLZ

st_read("./Data/Closure areas/BLZwestof162.shp") %>%
  vect() -> westBLZ

st_read("./Data/Closure areas/NBBTCA.shp") %>%
  vect() -> NBBTCA


# OTHER PARAMETERS --------------------------------------------------------------
lm_iter <- 8
thres <- 0.6626053

readRDS("./Models/nb_model.rda") -> nb_model
theta <- nb_model$theta # dispersion parameter from negative binomial model
