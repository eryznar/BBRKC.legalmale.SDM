### LOAD PACKAGES -----------------------------------------------------------------------------------------
library(sf)
library(ggmap)
library(rgdal)
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

library(dismo)
library(gbm) 
library(pROC)
library(ggrepel)
library(geosphere)
library(usdm)

### SET SPATIAL DETAILS ---------------------------------------------------------
crs.latlon <- "epsg:4326" #lat lon crs

map.crs <- "EPSG:3338"

in.crs = "+proj=longlat +datum=NAD83"

# LOAD SPATIAL LAYERS -----------------------------------------------------------
region_layers <- akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto")

survey_gdb <- "./Data/SAP_layers.gdb"

readOGR(dsn=survey_gdb,layer="BristolBaySurveyStrata") %>%
  vect(crs = crs.latlon) %>%
  project(map.crs) -> BB_strata


st_read("./Data/Closure areas/RKCSA_sub.shp") %>%
  vect() -> RKCSA_sub

st_read("./Data/Closure areas/RKCSA.shp") %>%
  vect() -> RKCSA

st_read("./Data/Closure areas/area512.shp") %>%
  vect() -> area512


# OTHER PARAMETERS --------------------------------------------------------------
lm_iter <- 8
thres <- 0.6626053