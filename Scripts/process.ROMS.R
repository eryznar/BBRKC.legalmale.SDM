#PURPOSE:
# To process ROMs bottom temperature data

#AUTHOR:
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)

### LOAD PROCESSING PARAMETERS --------------------------------------------------
source("./Scripts/load.libs.params.R")

## prepare bathymetry data of Bering Sea region
library(rgdal)
library(dplyr)
library(ggplot2)
library(raster)
library(sf)
library(stringr)
library(crabber)
library(magrittr)
library(terra)
library(ncmeta)
library(RNetCDF)
library(tsibble)
library(lubridate)
library(stars)
library(cubelyr)
library(crabber)

# Mapping crs and layers
map.crs <- "EPSG:3338"
in.crs = "+proj=longlat +datum=NAD83"


# Shapefiles for management boundaries and closure areas
survey_gdb<- "./Data/SAP_layers.gdb"

readOGR(dsn=survey_gdb,layer="BristolBaySurveyStrata") %>%
  vect(crs = in.crs) %>%
  project(map.crs) -> BB_strata

BB_strata <- st_as_sf(BB_strata)

map_layers <- akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto")



# BOTTOM TEMPERATURE ------------------------------------------------------------------------------
  # Specify files (TRY download.files using URLs so you don't have to download manually)
  files <- c(c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_1995-1999_average_temp_bottom5m.nc"),
             c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2000-2004_average_temp_bottom5m.nc"),
             c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2005-2009_average_temp_bottom5m.nc"),
             c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2010-2014_average_temp_bottom5m.nc"),
             c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2015-2019_average_temp_bottom5m.nc"),
             c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2020-2024_average_temp_bottom5m.nc"))
  
  
  bt1 <- stars::read_ncdf(files[1], var="temp",
                        curvilinear = c("lon_rho",
                                        "lat_rho")) %>%
    c(., stars::read_ncdf(files[3], var="temp",
                          curvilinear = c("lon_rho",
                                          "lat_rho")),
      along = 3) %>%
    c(., stars::read_ncdf(files[5], var="temp",
                          curvilinear = c("lon_rho",
                                          "lat_rho")),
      along = 3) %>%
    st_set_crs(., "EPSG:4326") %>% #assign CRS
    st_transform(.,map.crs) %>%
    .[BB_strata] 
  
   bt.00_04 <- stars::read_ncdf(files[2], var="temp",
                                curvilinear = c("lon_rho",
                                                "lat_rho")) %>%
               st_set_crs(., "EPSG:4326") %>% #assign CRS
               st_transform(.,map.crs) %>%
               .[BB_strata]
   attr(bt.00_04,"dimensions")[[3]]$values <- st_get_dimension_values(bt.00_04, "ocean_time")
   
   
   bt.10_14 <- stars::read_ncdf(files[4], var="temp",
                                curvilinear = c("lon_rho",
                                                "lat_rho")) %>%
     st_set_crs(., "EPSG:4326") %>% #assign CRS
     st_transform(.,map.crs) %>%
     .[BB_strata] 
   attr(bt.10_14,"dimensions")[[3]]$values <- st_get_dimension_values(bt.10_14, "ocean_time")
   
  
   bt.20_24 <-  stars::read_ncdf(files[6], var="temp",
                             curvilinear = c("lon_rho",
                                             "lat_rho")) %>%
      st_set_crs(., "EPSG:4326") %>% #assign CRS
      st_transform(.,map.crs) %>%
      .[BB_strata] 
    attr(bt.20_24,"dimensions")[[3]]$values <- st_get_dimension_values(bt.20_24, "ocean_time")
  
    t <-  c(bt1, bt.00_04, bt.10_14, bt.20_24, along = 3)
    
    # Pull out dates
      roms.dates <- tibble(roms_dates = st_get_dimension_values(t,
                                                                "ocean_time"),
                           yw = yearweek(roms_dates),
                           year = year(roms_dates),
                           month = month(roms_dates),
                           period = case_when((month %in% 1:2) ~ "Jan/Feb",
                                              (month %in% 3:4) ~ "Mar/Apr",
                                              (month %in% 5:6) ~ "May/Jun",
                                              (month %in% 7:8) ~ "Jul/Aug",
                                              (month %in% 9:10) ~ "Sep/Oct",
                                              (month %in% 11:12) ~ "Nov/Dec"))
      
      roms.dates %>%
        filter(year %in% 1997:2024) -> roms.dates2
    
    index<- which(st_get_dimension_values(t, "ocean_time") %in% roms.dates2$roms_dates)
    
    bt.dat <- data.frame()
    
    for(ii in 1:length(index)){
        
        print(paste0("week ", (1:length(index))[ii], "/", length(index)))
        roms.init <- t[,,,index[ii]]
        
        time.stamp <- st_get_dimension_values(roms.init, "ocean_time")
        
        yr <- year(time.stamp)
        month <- month(time.stamp)
        
        if(month %in% 1:2){
          period = "Jan/Feb"
        }else if(month %in% 3:4){
          period = "Mar/Apr"
        }else if(month %in% 5:6){
          period = "May/Jun"
        }else if(month %in% 7:8){
          period = "Jul/Aug"
        }else if(month %in% 9:10){
          period = "Sep/Oct"
        }else{
          period = "Nov/Dec"
        }
        
        
        dat <- data.frame(as.data.frame(roms.init),
                          year = yr,
                          month = month,
                          period = period) %>%
          dplyr::select(!ocean_time) %>%
          rename(lon = xi_rho, lat = eta_rho) %>%
          group_by(year, period, lon, lat) %>%
          reframe(mean.temp = mean(as.numeric(temp))) %>%
          na.omit()
        
        rbind(bt.dat, dat) -> bt.dat
        # 
        # %>%
        #   sf::st_as_sf(coords = c(x = "lon", y = "lat"),
        #                crs = map.crs)
        # 
        # lm_df %>%
        #   filter(year == yr) %>%
        #   st_as_sf(.,coords = c("longitude","latitude"),crs = 4326) %>%
        #   st_transform(.,map.crs) -> lm.init
        # 
        # mapped_init <-
        #   st_extract(roms.init, lm.init) %>%
        #   st_set_geometry(NULL) %>%
        #   bind_cols(.,survey_init_df)%>%
        #   dplyr::rename(roms_temp = temp) %>%
        #   mutate(roms_temp = as.numeric(roms_temp))
      }
  
    saveRDS(bt.dat, "./Data/EFH Legal Male/ROMs.data/BristolBay.bottomtemp.rda")
    
# CURRENT EAST ----------------------------------------------------------------------
  # Specify files
  files<- c(c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_1995-1999_average_uEast_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2000-2004_average_uEast_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2005-2009_average_uEast_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2010-2014_average_uEast_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2015-2019_average_uEast_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2020-2024_average_uEast_bottom5m.nc"))
    
 # Read in files
  uEast1 <- stars::read_ncdf(files[1], var="uEast") %>%
        c(., stars::read_ncdf(files[3], var="uEast"),
        along = 3) %>%
      c(., stars::read_ncdf(files[5], var="uEast"),
        along = 3)
  
  uEast.00_04 <- stars::read_ncdf(files[2], var="uEast") 
  attr(uEast.00_04,"dimensions")[[3]]$values <- st_get_dimension_values(uEast.00_04, "ocean_time")

  uEast.10_14 <- stars::read_ncdf(files[4], var="uEast") 
  attr(uEast.10_14,"dimensions")[[3]]$values <- st_get_dimension_values(uEast.10_14, "ocean_time")
  
  uEast.20_24 <- stars::read_ncdf(files[6], var="uEast") 
  attr(uEast.20_24,"dimensions")[[3]]$values <- st_get_dimension_values(uEast.20_24, "ocean_time")
  
  uEast2 <- c(uEast1, uEast.00_04, uEast.10_14, uEast.20_24, along =3)
  
 # Current doesn't have lat/lon dimensions, specify dimensions from bt file
  st_dimensions(uEast2) <- st_dimensions(t)
  uEast <- uEast2 %>%
            .[BB_strata]
  
  # Pull out time information
  roms.dates <- tibble(roms_dates = st_get_dimension_values(uEast,"ocean_time"),
                       yw = yearweek(roms_dates),
                       year = year(roms_dates),
                       month = month(roms_dates),
                       period = case_when((month %in% 1:2) ~ "Jan/Feb",
                                          (month %in% 3:4) ~ "Mar/Apr",
                                          (month %in% 5:6) ~ "May/Jun",
                                          (month %in% 7:8) ~ "Jul/Aug",
                                          (month %in% 9:10) ~ "Sep/Oct",
                                          (month %in% 11:12) ~ "Nov/Dec"))
  
  # Filter dates by years of interest, create time index to subset ROMs data
  roms.dates %>%
    filter(year %in% 1997:2023) -> roms.dates2
  
  index<- which(st_get_dimension_values(uEast, "ocean_time") %in% roms.dates2$roms_dates)
  
  # Loop through ROMs stars objects, pull out weekly data, average by time period
  uEast.dat <- data.frame()
  
  for(ii in 1:length(index)){
    
    print(paste0("week ", (1:length(index))[ii], "/", length(index)))
    roms.init <- uEast[,,,index[ii]]
    
    time.stamp <- st_get_dimension_values(roms.init, "ocean_time")
    
    yr <- year(time.stamp)
    month <- month(time.stamp)
    
    if(month %in% 1:2){
      period = "Jan/Feb"
    }else if(month %in% 3:4){
      period = "Mar/Apr"
    }else if(month %in% 5:6){
      period = "May/Jun"
    }else if(month %in% 7:8){
      period = "Jul/Aug"
    }else if(month %in% 9:10){
      period = "Sep/Oct"
    }else{
      period = "Nov/Dec"
    }
    
    
    dat <- data.frame(as.data.frame(roms.init),
                      year = yr,
                      month = month,
                      period = period) %>%
      dplyr::select(!ocean_time) %>%
      rename(lon = xi_rho, lat = eta_rho) %>%
      group_by(year, period, lon, lat) %>%
      reframe(mean.uEast = mean(as.numeric(uEast))) %>%
      na.omit()
    
    rbind(uEast.dat, dat) -> uEast.dat
   
  }
  
  saveRDS(uEast.dat, "./Data/EFH Legal Male/ROMs.data/BristolBay.uEast.rda")
  
  
# CURRENT North -----------------------------------------------------------------------------
  # Specify files
  files<- c(c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_1995-1999_average_vNorth_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2000-2004_average_vNorth_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2005-2009_average_vNorth_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2010-2014_average_vNorth_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2015-2019_average_vNorth_bottom5m.nc"),
            c("./Data/EFH Legal Male/ROMS.data/B10K-K20_CORECFS_2020-2024_average_vNorth_bottom5m.nc"))
  
  # Read in files
  vNorth1 <- stars::read_ncdf(files[1], var="vNorth") %>%
    c(., stars::read_ncdf(files[3], var="vNorth"),
      along = 3) %>%
    c(., stars::read_ncdf(files[5], var="vNorth"),
      along = 3)
  
  vNorth.00_04 <- stars::read_ncdf(files[2], var="vNorth") 
  attr(vNorth.00_04,"dimensions")[[3]]$values <- st_get_dimension_values(vNorth.00_04, "ocean_time")
  
  vNorth.10_14 <- stars::read_ncdf(files[4], var="vNorth") 
  attr(vNorth.10_14,"dimensions")[[3]]$values <- st_get_dimension_values(vNorth.10_14, "ocean_time")
  
  vNorth.20_24 <- stars::read_ncdf(files[6], var="vNorth") 
  attr(vNorth.20_24,"dimensions")[[3]]$values <- st_get_dimension_values(vNorth.20_24, "ocean_time")
  
  vNorth2 <- c(vNorth1, vNorth.00_04, vNorth.10_14, vNorth.20_24, along =3)
  
  # Current doesn't have lat/lon dimensions, specify dimensions from bt file
  st_dimensions(vNorth2) <- st_dimensions(t)
  vNorth <- vNorth2 %>%
    .[BB_strata]
  
  # Pull out time information
  roms.dates <- tibble(roms_dates = st_get_dimension_values(vNorth,"ocean_time"),
                       yw = yearweek(roms_dates),
                       year = year(roms_dates),
                       month = month(roms_dates),
                       period = case_when((month %in% 1:2) ~ "Jan/Feb",
                                          (month %in% 3:4) ~ "Mar/Apr",
                                          (month %in% 5:6) ~ "May/Jun",
                                          (month %in% 7:8) ~ "Jul/Aug",
                                          (month %in% 9:10) ~ "Sep/Oct",
                                          (month %in% 11:12) ~ "Nov/Dec"))
  
  # Filter dates by years of interest, create time index to subset ROMs data
  roms.dates %>%
    filter(year %in% 1997:2023) -> roms.dates2
  
  index<- which(st_get_dimension_values(vNorth, "ocean_time") %in% roms.dates2$roms_dates)
  
  # Loop through ROMs stars objects, pull out weekly data, average by time period
  vNorth.dat <- data.frame()
  
  for(ii in 1:length(index)){
    
    print(paste0("week ", (1:length(index))[ii], "/", length(index)))
    roms.init <- vNorth[,,,index[ii]]
    
    time.stamp <- st_get_dimension_values(roms.init, "ocean_time")
    
    yr <- year(time.stamp)
    month <- month(time.stamp)
    
    if(month %in% 1:2){
      period = "Jan/Feb"
    }else if(month %in% 3:4){
      period = "Mar/Apr"
    }else if(month %in% 5:6){
      period = "May/Jun"
    }else if(month %in% 7:8){
      period = "Jul/Aug"
    }else if(month %in% 9:10){
      period = "Sep/Oct"
    }else{
      period = "Nov/Dec"
    }
    
    
    dat <- data.frame(as.data.frame(roms.init),
                      year = yr,
                      month = month,
                      period = period) %>%
      dplyr::select(!ocean_time) %>%
      rename(lon = xi_rho, lat = eta_rho) %>%
      group_by(year, period, lon, lat) %>%
      reframe(mean.vNorth = mean(as.numeric(vNorth))) %>%
      na.omit()
    
    rbind(vNorth.dat, dat) -> vNorth.dat
    
  }
  
  saveRDS(vNorth.dat, "./Data/EFH Legal Male/ROMs.data/BristolBay.vNorth.rda")
  
  
# INTERPOLATE ROMs ONTO REGULAR GRIDS ----------------------------------------------

readRDS("./Data/EFH Legal Male/ROMs.data/BristolBay.bottomtemp.rda") -> bt.dat
  
  
  yr = c(1997:2023)
  
 interp.ROMs <- function(dat, yrs, prd){
   for(ii in 1:length(yrs)){
     dat %>%
       filter(year== yrs[1], period == prd) %>%
       st_as_sf(., coords = c("lon", "lat"), crs= map.crs) %>%
       rename(val = names(.)[3])-> dat2

     # Inverse distance weighting
     idw_fit <- gstat::gstat(formula = val~1, locations = dat2, nmax = 4)

     # Predict station points
     stn.predict <- predict(idw_fit, dat2)

     # Generate extrapolation grid
     extrap.box = c(xmn = -179.5, xmx = -157, ymn = 50, ymx = 68)
     extrap.box = c(ymn = 54.25, ymx = 59.25, xmn = -167.5, xmx = -158)
     grid.cell = c(0.001, 0.001)
     
     sp_extrap.raster <- raster::raster(xmn = extrap.box['xmn'],
                                        xmx=extrap.box['xmx'],
                                        ymn=extrap.box['ymn'],
                                        ymx=extrap.box['ymx'],
                                        ncol=(extrap.box['xmx']-extrap.box['xmn'])/grid.cell,
                                        nrow=(extrap.box['ymx']-extrap.box['ymn'])/grid.cell,
                                        crs = raster::crs("+proj=longlat")) %>% 
       raster::projectRaster(crs = raster::crs(dat2))
     
     res(sp_extrap.raster) <- 100

     # Predict, rasterize, mask
     extrap.grid <- predict(idw_fit, as(sp_extrap.raster, "SpatialPoints")) %>%
       sf::st_as_sf() %>%
       sf::st_transform(crs = raster::crs(dat2)) %>%
       stars::st_rasterize() %>%
       sf::st_join(BB_strata, join = st_intersects)
     
     r<- rast(extrap.grid)
     BB_extent <- c(-900000, -190000, 539823, 1050000)
     crop(r, BB_extent) -> rr
     plot(rr$var1.pred_lyr.1)
     # 
     # Set up a new IDW for ordinary kriging ----
     idw_vgm_fit <- gstat::gstat(formula = val ~ 1,
                                 locations = dat2,
                                 nmax = Inf)

     # Ordinary Kriging: Stein's Matern VGM----
     ste.vgfit <- gstat::fit.variogram(variogram(idw_vgm_fit),
                                       vgm(c("Gau")))

     ste_fit <- gstat::gstat(formula = val ~ 1,
                             locations = dat2,
                             model = ste.vgfit,
                             nmax = Inf)
     extrap.grid <- predict(ste_fit, as(sp_extrap.raster, "SpatialGrid")) %>%
                              sf::st_as_sf() %>%
                              sf::st_transform(crs = raster::crs(dat2)) %>%
                              stars::st_rasterize() %>%
                              sf::st_join(BB_strata, join = st_intersects)

     r <- rast(extrap.grid)$var1.pred_lyr.1
     
     names(r) <- paste0(prd, " ", names(dat)[5], " ", yrs[ii])
     
    suppressWarnings(c(roms.rast, r) -> roms.rast)
    
    print(names(r))
   }
   
   return(roms.rast)
 }
 
 yrs <- 1997:2023
 
interp.ROMs(bt.dat, 2020, "Jan/Feb") -> hh
 
 roms.rast <- rast()
 c("Jan/Feb", "Mar/Apr", "May/Jun", "Jul/Aug") %>%
    purrr::map(~interp.ROMs(bt.dat, yrs, .x)) -> bt.out
 
 roms.rast <- rast()
 c("Sep/Oct", "Nov/Dec") %>%
   purrr::map(~interp.ROMs(bt.dat, 1997:2022, .x)) -> bt.out.FALL
 
 
 roms.rast <- rast()
 c("Jan/Feb", "Mar/Apr", "May/Jun", "Jul/Aug") %>%
   purrr::map(~interp.ROMs(uEast.dat, yrs, .x)) -> uEast.out
 
 roms.rast <- rast()
 c("Sep/Oct", "Nov/Dec") %>%
   purrr::map(~interp.ROMs(uEast.dat, 1997:2022, .x)) -> uEast.out.FALL
 
 
 roms.rast <- rast()
 c("Jan/Feb", "Mar/Apr", "May/Jun", "Jul/Aug") %>%
   purrr::map(~interp.ROMs(vNorth.dat, yrs, .x)) -> vNorth.out
 
 roms.rast <- rast()
 c("Sep/Oct", "Nov/Dec") %>%
   purrr::map(~interp.ROMs(vNorth.dat, 1997:2022, .x)) -> vNorth.out.FALL
 
 
 c(rast(bt.out), rast(bt.out.FALL)) -> bt.rast
 c(rast(uEast.out), rast(uEast.out.FALL)) -> uEast.rast
 c(rast(vNorth.out), rast(vNorth.out.FALL)) -> vNorth.rast
 
 
 writeRaster(bt.rast, "./Data/EFH Legal Male/ROMs.data/BristolBay.bottomtemp.tif", overwrite = TRUE)
 writeRaster(uEast.rast, "./Data/EFH Legal Male/ROMs.data/BristolBay.uEast.tif", overwrite = TRUE)
 writeRaster(vNorth.rast, "./Data/EFH Legal Male/ROMs.data/BristolBay.vNorth.tif", overwrite = TRUE)
 
 
 
 
 
  
  bt.dat %>%
    filter(year ==1997, period == "Jan/Feb")-> ll
  
 r <- rasterize(ll[, 3:4], sp_extrap.raster, ll[,5], fun=mean)
  