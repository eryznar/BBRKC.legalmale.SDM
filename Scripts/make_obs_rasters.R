source("./Scripts/load.libs.params.R")

# Read in response data
lm_df <- read.csv("./Data/legalmale_direct.fish.csv") %>%
  filter(season == "F")

# Read in spatial covariates
Fall_lm.preds <- rast("./Data/Fall_lm.preds.tif")
samp <- Fall_lm.preds$`Jan/Feb mean.temp 1997`

x <- raster(ncol=36, nrow=18, xmn=-1000, xmx=1000, ymn=-100, ymx=900)

x <- rast(ext(samp), resolution=res(samp))
crs(x) <- crs(samp)


as.data.frame(crds(samp)) -> xy

xy$lon <- xy$x
xy$lat <- xy$y

rasterFromXYZ(xy[,-4], crs = crs(samp)) -> xx
xx <- rep(rast(xx), 27)
names(xx) <- paste0("x", 1997:2023)
rasterFromXYZ(xy[,-3], crs = crs(samp)) -> yy
yy <- rep(rast(yy), 27)
names(yy) <- paste0("y", 1997:2023)

c(xx, yy) -> xy_rasts
ext(xy_rasts) <- ext(samp)

resample(xy_rasts, Fall_lm.preds) -> xy_rasts
c(xy_rasts, Fall_lm.preds) -> Fall_lm.preds

writeRaster(xy_rasts, "./Data/xy_rasts.tif", overwrite = TRUE)
writeRaster(Fall_lm.preds, "./Data/Fall_lm.predsxy.tif", overwrite = TRUE)

predict_yr <- 1997:2023
obs_rast <- c()
for(ii in 1:length(predict_yr)){
  lm_df %>%
    filter(year == predict_yr[ii]) %>%
    dplyr::select(longitude, latitude, catch_pp) %>%
    dplyr::rename(x = longitude, y= latitude) %>%
    st_as_sf(., coords = c("x", "y"), crs = in.crs) %>%
    st_transform(., crs(samp)) -> pp
  
  as.data.frame(cbind(st_coordinates(pp), pp))-> ll
  
  ll <- as.matrix(ll[,1:3])
  
  
  # set up a raster geometry, here deriving an extent from your data
  e <- ext(samp)
  # set up the raster, for example
  r <- rast(e, ncol=ncol(samp), nrow=nrow(samp), crs = map.crs)
  terra::project(r, crs(samp)) -> r
  
  
  # you need to provide a function 'fun' for when there are multiple points per cell
  oo <- rasterize(ll[,1:2], r, ll[,3], fun=mean)
  
  names(oo) <- paste("obs", predict_yr[ii])
  
  obs_rast <- c(obs_rast, oo)
}

obs_rast <- rast(obs_rast)
writeRaster(obs_rast, "./Data/obs_rasters.tif", overwrite = TRUE)
