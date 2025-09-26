#PURPOSE: To generate rasters of influential temperature covariates for fall legal male SDM of average
#         of minimum and maximum 5 values per grid cell across years

# AUTHOR: 
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)


### LOAD PACKAGES -----------------------------------------------------------------------------------------
source("./Scripts/load.libs.params.R")

### LOAD DATA -------------------------------------------------------------------
Fall_lm.preds <- rast("./Data/EFH Legal Male/Fall_lm.preds.tif")

btemp <- subset(Fall_lm.preds, grep("mean.temp ", names(Fall_lm.preds)))
#btemp2 <- subset(btemp, grep("May/Jun |Mar/Apr ", names(btemp), invert = TRUE)) # drop collinear
ice <- subset(Fall_lm.preds, grep("Ice ", names(Fall_lm.preds)))
currentE <- subset(Fall_lm.preds, grep("uEast ", names(Fall_lm.preds)))
#currentE2 <- subset(currentE, grep("Jul/Aug ", names(currentE), invert = TRUE)) # drop collinear
currentN <- subset(Fall_lm.preds, grep("vNorth ", names(Fall_lm.preds)))
SAP <- subset(Fall_lm.preds, grep("SAP Count ", names(Fall_lm.preds)))



SOtemp <- subset(Fall_lm.preds, grep("Sep_Oct ", names(Fall_lm.preds)))
JAtemp <- subset(Fall_lm.preds, grep("Jul_Aug ", names(Fall_lm.preds)))
JFice <- subset(Fall_lm.preds, grep("Jan_Feb Ice ", names(Fall_lm.preds)))
MAice <- subset(Fall_lm.preds, grep("Mar_Apr Ice ", names(Fall_lm.preds)))
SAP <- subset(Fall_lm.preds, grep("SAP Count ", names(Fall_lm.preds)))


yrs <- c(1997:2019, 2021:2023)
yrs2 <- c(1998:2019, 2021:2023)

# BOTTOM TEMPERATURE -------------------------------------------------------------
    btemp_df <- data.frame()
    btemp.max_df <- data.frame()
    btemp.min_df <- data.frame()
    btemp_rasts <- c()
 
 prd <- c("Jan/Feb",  "Mar/Apr", "May/Jun", "Jul/Aug", "Sep/Oct", "Nov/Dec")
 
 for(jj in 1:length(prd)){
   btemp2 <- subset(btemp, grep(prd[jj], names(btemp)))
   for(ii in 1:nlyr(btemp2)){
  
     df <- cbind(crds(btemp2[[ii]]), as.data.frame(btemp2[[ii]]), yrs[ii])
     names(df)[3:4] <- c("val", "year")
     
     
     rbind(df, btemp_df) -> btemp_df
     
   }
   
   btemp_df %>%
     filter(year %in% yrs2) %>%
     group_by(x, y) %>%
     top_n(5, val) %>%
     mutate(period = prd[jj]) -> btemp.max
   
   btemp.max %>%
     group_by(x, y) %>%
     reframe(avg.max = mean(val)) %>%
     rast() -> btemp.max.rast
   
   btemp_df %>%
     na.omit() %>%
     filter(year %in% yrs2) %>%
     group_by(x, y) %>%
     top_n(-5, val) %>%
     mutate(period = prd[jj]) -> btemp.min
   
   btemp.min %>%
     na.omit() %>%
     group_by(x, y) %>%
     reframe(avg.min = mean(val)) %>%
     rast()  -> btemp.min.rast
   
   btemp.max_df <- rbind(btemp.max_df, btemp.max)
   btemp.min_df <- rbind(btemp.min_df, btemp.min)
  
   rasts <- c(btemp.max.rast, btemp.min.rast)
   names(rasts) <- paste0(prd[jj], c(" BT max", " BT min"))
   
   btemp_rasts <- c(btemp_rasts, rasts)
   
   print(prd[jj])
 }   
   
  btemp.rr <- rast(btemp_rasts)  
  
# CURRENT E ----------------------------------------------------------------------
 cE_df <- data.frame()
 cE_rasts <- c()
 
 prd <- c("Jan/Feb",  "Mar/Apr", "May/Jun", "Jul/Aug", "Sep/Oct", "Nov/Dec")
 
 for(jj in 1:length(prd)){
   cE3 <- subset(currentE, grep(prd[jj], names(currentE)))
   for(ii in 1:nlyr(cE3)){
     
     df <- cbind(crds(cE3[[ii]]), as.data.frame(cE3[[ii]]), yrs[ii], prd[jj])
     names(df)[3:5] <- c("val", "year", "period")
     
     
     rbind(df, cE_df) -> cE_df
     
   }
   
   cE_df %>%
     na.omit() %>%
     rename(cEast = val) %>%
     right_join(., btemp.max_df %>% filter(period == prd[jj]), by = c("x", "y", "year", "period")) %>%
     dplyr::select(!val) -> cE.max
   
   cE.max %>%
     group_by(x, y) %>%
     reframe(avg.max = mean(cEast)) %>%
     rast() -> cE.max.rast
   
   cE_df %>%
     na.omit() %>%
     rename(cEast = val) %>%
     right_join(., btemp.min_df %>% filter(period == prd[jj]), by = c("x", "y", "year", "period")) %>%
     dplyr::select(!val) -> cE.min
   
   cE.min %>%
     group_by(x, y) %>%
     reframe(avg.min = mean(cEast)) %>%
     rast()  -> cE.min.rast
   
   rasts <- c(cE.max.rast, cE.min.rast)
   names(rasts) <- paste0(prd[jj], c(" currentE max", " currentE min"))
   
   cE_rasts <- c(cE_rasts, rasts)
   
   print(prd[jj])
 }   
 
 cE.rr <- rast(cE_rasts)
 
# CURRENT N ----------------------------------------------------------------------
 cN_df <- data.frame()
 cN_rasts <- c()
 
 prd <- c("Jan/Feb",  "Mar/Apr", "May/Jun", "Jul/Aug", "Sep/Oct", "Nov/Dec")
 
 for(jj in 1:length(prd)){
   cN2 <- subset(currentN, grep(prd[jj], names(currentN)))
   for(ii in 1:nlyr(cN2)){
     
     df <- cbind(crds(cN2[[ii]]), as.data.frame(cN2[[ii]]), yrs[ii], prd[jj])
     names(df)[3:5] <- c("val", "year", "period")
     
     
     rbind(df, cN_df) -> cN_df
     
   }
   
   cN_df %>%
     na.omit() %>%
     rename(cNorth = val) %>%
     right_join(., btemp.max_df %>% filter(period == prd[jj]), by = c("x", "y", "year", "period")) %>%
     dplyr::select(!val) -> cN.max
   
   cN.max %>%
     group_by(x, y) %>%
     reframe(avg.max = mean(cNorth)) %>%
     rast() -> cN.max.rast
   
   cN_df %>%
     na.omit() %>%
     rename(cNorth = val) %>%
     right_join(., btemp.min_df %>% filter(period == prd[jj]), by = c("x", "y", "year", "period")) %>%
     dplyr::select(!val) -> cN.min
   
   cN.min %>%
     group_by(x, y) %>%
     reframe(avg.min = mean(cNorth)) %>%
     rast()  -> cN.min.rast
   
   rasts <- c(cN.max.rast, cN.min.rast)
   names(rasts) <- paste0(prd[jj], c(" currentN max", " currentN min"))
   
   cN_rasts <- c(cN_rasts, rasts)
   
   print(prd[jj])
 }   
 
 cN.rr <- rast(cN_rasts)
 
# ICE ----------------------------------------------------------------------------
 ice_df <- data.frame()
 ice_rasts <- c()
 
 prd <- c("Jan_Feb",  "Mar_Apr")
 
 for(jj in 1:length(prd)){
   ice2 <- subset(ice, grep(prd[jj], names(ice)))
   for(ii in 1:nlyr(ice2)){
     
     df <- cbind(crds(ice2[[ii]]), as.data.frame(ice2[[ii]]), yrs[ii], prd[jj])
     names(df)[3:5] <- c("val", "year", "period")
     
     
     rbind(df, ice_df) -> ice_df
     
   }
   
   ice_df %>%
     na.omit() %>%
     filter(year %in% yrs2) %>%
     group_by(x, y) %>%
     top_n(5, val) %>%
     mutate(period = prd[jj]) -> ice.max
   
   ice.max %>%
     na.omit() %>%
     group_by(x, y) %>%
     reframe(avg.min = mean(val)) %>%
     rast()  -> ice.max.rast
   
   ice_df %>%
     na.omit() %>%
     filter(year %in% yrs2) %>%
     group_by(x, y) %>%
     top_n(-5, val) %>%
     mutate(period = prd[jj]) -> ice.min
   
   ice.min %>%
     na.omit() %>%
     group_by(x, y) %>%
     reframe(avg.min = mean(val)) %>%
     rast()  -> ice.min.rast
   
   rasts <- c(ice.max.rast, ice.min.rast)
   names(rasts) <- paste0(prd[jj], c(" ice max", " ice min"))
   
   ice_rasts <- c(ice_rasts, rasts)
   
   print(prd[jj])
 }   
 
 ice.rr <- rast(ice_rasts)
 names(ice.rr) <- c("Jan/Feb ice max", "Jan/Feb ice min", "Mar/Apr ice max", "Mar/Apr ice min")
 
# BBRKC summer survey CPUE --------------------------------------------------------------------
  SAP_df <- data.frame()
 yrs <- c(1996:2019, 2021:2023)
 
  for(ii in 1:nlyr(SAP)){
    
    df <- cbind(crds(SAP[[ii]]), as.data.frame(SAP[[ii]]), yrs[ii])
    names(df)[3:4] <- c("val", "year")
    
    
    rbind(df, SAP_df) -> SAP_df
    
  }
 
  SAP_df %>%
    replace_na(list(year = 1996)) %>%
    filter(year %in% yrs2) -> SAP_df
 
  rbind(SAP_df %>%
          rename(CPUE = val) %>%
          right_join(., btemp.max_df %>% filter(period == "May/Jun"), by = c("x", "y", "year")) %>%
          dplyr::select(!val),
        SAP_df %>%
          rename(CPUE = val) %>%
          right_join(., btemp.max_df %>% filter(period == "Jul/Aug"), by = c("x", "y", "year")) %>%
          dplyr::select(!val)) %>%
    group_by(x, y) %>%
    reframe(avg.max = mean(CPUE)) %>%
    rast() -> SAP.max.rast
  
  rbind(SAP_df %>%
          rename(CPUE = val) %>%
          right_join(., btemp.min_df %>% filter(period == "May/Jun"), by = c("x", "y", "year")) %>%
          dplyr::select(!val),
        SAP_df %>%
          rename(CPUE = val) %>%
          right_join(., btemp.min_df %>% filter(period == "Jul/Aug"), by = c("x", "y", "year")) %>%
          dplyr::select(!val)) %>%
    group_by(x, y) %>%
    reframe(avg.min = mean(CPUE)) %>%
    rast() -> SAP.min.rast
  
  SAP.rast <- c(SAP.max.rast, SAP.min.rast)
  names(SAP.rast) <- c("SAP Count max", "SAP Count min")
  
# COMBINE ALL MAX/MIN RASTS ------------------------------------------------------
  
  writeRaster(btemp.rr, "./Data/EFH Legal Male/max.min.bottomtempNEW.tif", overwrite = TRUE)
  writeRaster(cE.rr, "./Data/EFH Legal Male/max.min.currentENEW.tif", overwrite = TRUE)
  writeRaster(cN.rr, "./Data/EFH Legal Male/max.min.currentNNEW.tif", overwrite = TRUE)
  writeRaster(ice.rr, "./Data/EFH Legal Male/max.min.iceNEW.tif", overwrite = TRUE)
  writeRaster(SAP.rast, "./Data/EFH Legal Male/max.min.SAPNEW.tif", overwrite = TRUE)
  