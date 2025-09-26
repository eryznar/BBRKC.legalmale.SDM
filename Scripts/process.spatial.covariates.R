#PURPOSE:
# To compile spatial data to date for BBRKC legal male SDMs and process as needed

#AUTHOR:
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)

### LOAD PROCESSING PARAMETERS --------------------------------------------------
source("./Scripts/load.libs.params.R")

### PROCESS SPATIAL COVARIATES ---------------------------------------------------------------------------------------
slope <- rast("./Data/Slope.grd")
currentE <- rast("./Data/BristolBay.uEast.tif")
currentN <- rast("./Data/BristolBay.vNorth.tif")
bt <- rast("./Data/BristolBay.bottomtemp.tif")
tidemax <- rast("./Data/Tmax.grd")
ice<- rast("./Data/ice.tif") #Jan-Feb and Mar-Apr
lm_sap<- rast("./Data/lm_sap_CPUE.EBS.tif") %>%
  mask(BB_strata)
sed <- rast("./Data/EBS_phi_1km.grd") #already has crs defined, same across all years
F_bycatch <- rast("./Data/gfbycatch_fall97.23.tif")
depth <- rast("./Data/efh_bathy_1km.tif")  %>%
  mask(BB_strata)

# Read in averaged min/max temp rasters and BBRKC abundance
bt.mm <- rast("./Data/max.min.bottomtempNEW.tif") # bottom temperature
cE.mm <- rast("./Data/max.min.currentENEW.tif") # east current
cN.mm <- rast("./Data/max.min.currentNNEW.tif") # east current
ice.mm <- rast("./Data/max.min.iceNEW.tif") # Jan/Feb and Mar/Apr ice
lm.mm <- rast("./Data/max.min.SAPNEW.tif")

#Set input crs for rasters that don't already have them
crs(ice) <- in.crs
crs(bt.mm) <- map.crs
crs(cE.mm) <- map.crs
crs(cN.mm) <- map.crs
crs(ice.mm) <- map.crs
crs(lm.mm) <- map.crs

#Project rasters to the same crs
slope2 <- terra::project(slope, map.crs)
tidemax2 <- terra::project(tidemax, map.crs)

ice2 <- terra::project(ice, map.crs)
sed2 <- terra::project(sed, map.crs)
lm_sap2 <- terra::project(lm_sap, map.crs)
F_bycatch2 <- terra::project(F_bycatch, map.crs)

bt.mm2 <- terra::project(bt.mm, map.crs)
cE.mm2 <- terra::project(cE.mm, map.crs)
cN.mm2 <- terra::project(cN.mm, map.crs)
ice.mm2 <- terra::project(ice.mm, map.crs)
lm.mm2 <- terra::project(lm.mm, map.crs)

#Set extents and crop rasters that don't match the coarsest raster
rast_ext <- c(-1500000, -170000, 539823, 1600000)

slope3 <- crop(slope2, rast_ext)
tidemax3 <- crop(tidemax2, rast_ext)
ice3 <- crop(ice2, rast_ext)
sed3 <- crop(sed2, rast_ext)
lm_sap3 <- crop(lm_sap2, rast_ext)
depth2 <- crop(depth, rast_ext)
F_bycatch3 <- crop(F_bycatch2, rast_ext)
bt2 <- crop(bt, rast_ext)
currentN2 <- crop(currentN, rast_ext)
currentE2 <- crop(currentE, rast_ext)

bt.mm3 <- crop(bt.mm2, rast_ext)
cE.mm3 <- crop(cE.mm2, rast_ext)
cN.mm3 <- crop(cN.mm2, rast_ext)
ice.mm3 <- crop(ice.mm2, rast_ext)
lm.mm3 <- crop(lm.mm2, rast_ext)


# Resample rasters to match resolution
resample(ice3, sed3) -> ice4
resample(bt2, sed3) -> bt3
resample(currentE2, sed3) -> currentE3
resample(currentN2, sed3) -> currentN3

resample(lm_sap3, sed3) -> lm_sap4
resample(F_bycatch3, sed3) -> F_bycatch4
resample(depth2, sed3) -> depth3

resample(bt.mm3, sed3) -> bt.mm4
resample(cE.mm3, sed3) -> cE.mm4
resample(cN.mm3, sed3) -> cN.mm4
resample(ice.mm3, sed3) -> ice.mm4
resample(lm.mm3, sed3) -> lm.mm4


#Create raster stack of Environmental variables with the same crs, resolution, and extent, crop to Bristol Bay
BB_extent <- c(-900000, -190000, 539823, 1050000)

years <- c(1997:2019, 2021:2023)

# Repeat stationary vars by number of years
slope4 <- rep(slope3, length(years))
tidemax4 <- rep(tidemax3, length(years))
depth3 <- rep(depth3, length(years))

sed3 <- rep(sed3, length(years))

# bt.mm4 <- rep(bt.mm4, length(years))
# cE.mm4 <- rep(cE.mm4, length(years))
# cN.mm4 <- rep(cN.mm4, length(years))
# ice.mm4 <- rep(ice.mm4, length(years))
# lm.mm4 <- rep(lm.mm4, length(years))
# 

# Set names for stationary vars
names(slope4) <- paste("Slope", years)
names(tidemax4) <- paste("Tidemax", years)
names(depth3) <- paste("Depth", years)
names(sed3) <- paste("Sed", years)

names(bt.mm4) <- paste(names(bt.mm4), 2023)
names(cE.mm4) <- paste(names(cE.mm4), 2023)
names(cN.mm4) <- paste(names(cN.mm4), 2023)
names(ice.mm4) <- paste(names(ice.mm4), 2023)
names(lm.mm4) <- paste(names(lm.mm4), 2023)

c(bt.mm4, cE.mm4, cN.mm4, ice.mm4, lm.mm4) -> mm_rasts

c(bt3, ice4, sed3, lm_sap4, depth3, F_bycatch4, slope4, currentE3, currentN3,
  tidemax4, mm_rasts) %>%
  terra::mask(vect(region_layers$survey.area)) %>%
  terra::mask(BB_strata) %>%
  crop(BB_extent) -> lm_preds.mm 
subset(lm_preds.mm, grep(" min| max", names(lm_preds.mm), invert = TRUE)) -> lm_preds

writeRaster(lm_preds, "./Data/Fall_lm.preds.tif", overwrite = TRUE)
writeRaster(lm_preds.mm, "./Data/Fall_lm.preds.mm.tif", overwrite = TRUE)

