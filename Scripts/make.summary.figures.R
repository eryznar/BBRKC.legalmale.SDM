# PURPOSE:
# To generate several summary plots for Bristol Bay red king crab legal male SDMs

# AUTHOR: 
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)


### LOAD PACKAGES -----------------------------------------------------------------------------------------
source("./Scripts/load.libs.params.R")

### LOAD DATA ----------------------------------------------------------------------------------------------
  # Read in spatial covariates
  Fall_lm.preds <- rast("./Data/Fall_lm.preds.tif")
  
  # Read in response data
  lm_df <- read.csv("./Data/legalmale_direct.fish.csv")
  
  # Read in training/testing data
  lm_F_train <- read.csv("./Output/lm_F_train.csv") %>%
    filter(iter == lm_iter)
  lm_F_test <- read.csv("./Output/lm_F_test.csv") %>%
    filter(iter == lm_iter)
  
  # Best models 
  lm.F.modelb <- readRDS("./Models/lm.modelb.F.8.rda")
  lm.F.modelp <- readRDS("./Models/lm.modelp.F.8.rda")

### LOAD FUNCTIONS ---------------------------------------------------------------------------------------
  # Prediction raster using raster stack of predictors
  MakePredictionRaster<- function(preds, model_b, model_p, resp_data, seas, predict_yr){
  
  # Plot 1: prediction rasters -----------------------------------------------------------------      
  # Loop through years, generate prediction rasters, make them into data frames
  
  spatpred_df <- data.frame()
  
  for(ii in 1:length(predict_yr)){
    
    preds <- subset(preds, grep("GFbycatch.F", names(preds), invert = TRUE)) # dropping gf bycatch bc it came out least important and the raster makes the spatial predictions look wonky
    
    # Specify past year based on present year and prediction period for subsetting
    pst <- ifelse(predict_yr[ii] == 2021, c(predict_yr[ii]-2, predict_yr[-length(predict_yr[ii])]), 
                  c(predict_yr[ii]-1, predict_yr[-length(predict_yr[ii])]))
    
    
    preds_past <- subset(preds, grep("Sep/Oct |Nov/Dec ", names(preds)))
    preds_pres <- subset(preds, grep("Sep/Oct |Nov/Dec ", names(preds), invert = TRUE))
    
    
    preds2 <- c(terra::subset(preds_past, grep(pst, names(preds_past), value = T)),
                terra::subset(preds_pres, grep(predict_yr[ii], names(preds_pres))))
    
    
    labs <- c("Sep_Oct_BT", "Nov_Dec_BT", "Sep_Oct_currentE", "Nov_Dec_currentE", "Sep_Oct_currentN",
              "Nov_Dec_currentN", 
              "Jan_Feb_BT", "Mar_Apr_BT", "May_Jun_BT", "Jul_Aug_BT", 
              "Jan_Feb_Ice", "Mar_Apr_Ice", 
              "Sed", "SAP_Count", "Depth", "Slope", 
              "Jan_Feb_currentE", "Mar_Apr_currentE", "May_Jun_currentE", "Jul_Aug_currentE",
              "Jan_Feb_currentN", "Mar_Apr_currentN", "May_Jun_currentN", "Jul_Aug_currentN",
              "Tidemax") #raster labs
    
    
    
    # Set names of preds
    names(preds2) <- labs
    
    # Create data frame for wind (since no spatial layer) to feed models by year
    data2 <- resp_data %>%
      filter(year == predict_yr[ii]) %>%
      mutate(fishery = as.factor(fishery))
    
    # Generate spatial model predictions for Bristol Bay management area extent using raster predictors
    #PA
    BB_strata_buff <- buffer(BB_strata, width = -1000)
    
    spatpred_b <- suppressWarnings(predict(preds2, # raster stack 
                                           model_b, # fitted model
                                           n.trees=model_b$gbm.call$best.trees, # see help
                                           #factors = f,
                                           type="response", 
                                           ext = panel_extent) %>%
                                     #mask(BB_strata, touches = FALSE) %>%
                                     mask(preds2$Jan_Feb_Ice) %>%
                                     mask(BB_strata_buff, touches = FALSE) %>%
                                     na.omit())
    
    
    #Abund
    spatpred_p <- suppressWarnings(predict(preds2, # raster stack 
                                           model_p, # fitted model
                                           n.trees=model_p$gbm.call$best.trees, # see help
                                           #factors = ff,
                                           type="response", 
                                           ext = panel_extent) %>%
                                     #mask(BB_strata, touches = FALSE)  %>%
                                     mask(preds2$Jan_Feb_Ice) %>%
                                     mask(BB_strata_buff, touches = FALSE) %>%
                                     na.omit())
    
    # Multiply binomial raster by abundance raster (part of delta model framework)
    spatpred <- spatpred_b * spatpred_p
    
    thres <- mean(na.omit(values(spatpred_b$lyr1)))
    
    # # Create data frame of spatial predictions (for pretty mapping)
    out <- cbind(crds(spatpred$lyr1), as.data.frame(spatpred$lyr1)) %>%
      dplyr::rename("catch_pp" = lyr1) %>%
      mutate(year = predict_yr[ii])
    
    # out <- cbind(crds(spatpred_b$lyr1), as.data.frame(spatpred_b$lyr1) %>% rename(PA = lyr1),
    #              as.data.frame(spatpred_p$lyr1) %>% rename(CPUE = lyr1)) %>%
    #   mutate(PA = ifelse(PA < thres, 0, 1),
    #          catch_pp = CPUE * PA,
    #          year = predict_yr[ii]) %>%
    #   dplyr::select(x, y, catch_pp, year)
    
    rbind(spatpred_df, out) -> spatpred_df
    
    pred_raster <- spatpred
  }
  
  # Plot prediction rasters by year
  # Specify plot boundaries
  plot.boundary.untrans <- data.frame(y = c(54.25, 59.25), 
                                      x = c(-167.5, -158)) # plot boundary unprojected
  
  plot.boundary <- plot.boundary.untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y) # plot boundary projected
  
  # Plot
  ggplot() +
    geom_tile(data = spatpred_df, aes(x = x, y = y, fill = catch_pp))+
    ggplot2::geom_sf(data = st_as_sf(area512), 
                     fill = NA, 
                     color = "white",
                     linewidth = 1)+
    ggplot2::geom_sf(data = st_as_sf(BB_strata), 
                     fill = NA, 
                     color = "black",
                     linewidth = 2)+
    
    ggplot2::geom_sf(data = region_layers$akland, 
                     fill = "grey70", 
                     color = "black")+
    
    ggplot2::geom_sf(data = st_as_sf(RKCSA),
                     fill = NA,
                     color = "white",
                     linewidth = 1)+
    ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                     fill = NA,
                     color = "white",
                     linewidth = 1)+
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    ggplot2::scale_x_continuous(name = "", 
                                breaks = c(-166, -162, -158), labels = paste0(c(166, 164, 162), "°W")) +
    scale_fill_viridis_c(option = "magma", name = "Count")+
    facet_wrap(~year)+
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.position = "bottom",
          panel.grid.major = element_blank(),
          legend.direction = "horizontal",
          plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 12)) -> pred_plot

  ggsave(plot = pred_plot, "./Figures/pred_fall.lm.dist.png", height = 10,
         width = 10, units = "in")

  
  # Plot 2: AVG EFH percentile maps -----------------------------------------------------------------------
  # # Avg habitat across all years
  # spatpred_df2 <- spatpred_df %>%
  #   group_by(x, y) %>%
  #   reframe(catch_pp = mean(catch_pp))
  # 
  # # Find plotting breaks
  # quantiles = c(.05, .25, .5, .75)
  # quants<-sort(unique(c(0,quantiles,1)))
  # 
  # threshold <- 0.0513
  # 
  # sample <- stats::na.omit(spatpred_df2$catch_pp)
  # sample[sample <= threshold] <- NA
  # perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
  # perc.breaks[1]<-0
  # perc.breaks[length(perc.breaks)]<-Inf
  # 
  # 
  # # Make raster again
  #  spatpred_df2 %>%
  #     rast() %>%
  #     raster() -> spatpred
  #   
  # # Set crs
  #   crs(spatpred) <- map.crs
  #   
  # # Cut the prediction map by the breaks
  #   perc.map <- raster::cut(spatpred, breaks = perc.breaks)
  #   
  # # set up the factor maps
  #   perc.vals <- raster::getValues(perc.map)
  #   perc.vals[perc.vals == 1] <- NA
  #   
  # # convert the raster to polygons
  #   percpoly0 <- stars::st_as_stars(perc.map)
  #   percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
  #   percpoly2 <- percpoly[percpoly$layer != 1, ]
  #   
  # # we'll need a new outline
  #   perc.dummy.raster <- raster::raster(perc.map)
  #   perc.vals2 <- is.na(perc.vals) == F
  #   perc.dummy.raster <- raster::setValues(perc.dummy.raster, values = perc.vals2)
  #   
  #   percdummy0 <- stars::st_as_stars(perc.dummy.raster)
  #   percdummy <- sf::st_cast(sf::st_as_sf(percdummy0, merge = TRUE))
  #   percdummy2 <- sf::st_transform(percdummy, sf::st_crs(map.crs))
  #   
  # # Dropping the smallest areas
  #   percdummy.poly <- sf::st_cast(percdummy2, "POLYGON")
  #   areas <- sf::st_area(percdummy.poly)
  #   
  #   outside <- order(areas, decreasing = T)[1]
  #   toosmall <- which(as.numeric(areas) < 10^8)
  #   
  #   perc.x <- percdummy2$layer[-c(outside, toosmall)]
  #   perc.y <- percdummy2$geometry[-c(outside, toosmall)]
  #   percdummy3 <- sf::st_sf(perc.x, perc.y) %>%
  #             vect() %>%
  #             crop(BB_strata)
  #   
  # # Set up plot boundary
  #   plot.boundary.untrans <- data.frame(y = c(54.25, 59.25), 
  #                                       x = c(-167.5, -158)) # plot boundary unprojected
  #   
  #   plot.boundary <- plot.boundary.untrans %>%
  #     sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  #     sf::st_transform(crs = map.crs) %>%
  #     sf::st_coordinates() %>%
  #     as.data.frame() %>%
  #     dplyr::rename(x = X, y = Y) # plot boundary projected
  #   
  # # Set up title
  # if(seas == "F"){
  #   title = "Fall Red King Crab Legal Male Encounter Probability"
  # } else{
  #   title = "Winter Red King Crab Legal Male Encounter Probability"
  # }
  #   
  # # Map
  #   ggplot2::ggplot() +
  #     #ggplot2::geom_sf(data = survey.sf, fill = "grey95")+
  #     ggplot2::geom_sf(data = percpoly2, ggplot2::aes(fill = as.factor(layer)), col = NA) +
  #     ggplot2::geom_sf(data = st_as_sf(percdummy3),fill=NA, size = .3) +
  #     ggplot2::geom_sf(data = st_as_sf(area512), 
  #                      fill = NA, 
  #                      color = "purple",
  #                      linewidth = 1.5)+
  #     ggplot2::geom_sf(data = st_as_sf(BB_strata), 
  #                      fill = NA, 
  #                      color = "black",
  #                      linewidth = 2)+
  #     ggplot2::geom_sf(data = region_layers$akland, 
  #                      fill = "grey70", 
  #                      color = "black")+
  #     ggplot2::geom_sf(data = st_as_sf(RKCSA),
  #                      fill = NA,
  #                      color = "red",
  #                      linewidth = 1.5)+
  #     ggtitle(title)+
  #     ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
  #                      fill = NA,
  #                      color = "red",
  #                      linewidth = 1.5)+
  #     coord_sf(xlim = plot.boundary$x,
  #              ylim = plot.boundary$y)+
  #     viridis::scale_fill_viridis(discrete = T, name = "Percentiles", labels = c("95%", "75%", "50%", "25%")) +
  #     ggplot2::theme_bw() +
  #     ggplot2:: theme(
  #       panel.border = ggplot2::element_rect(color = "black", fill = NA),
  #       panel.background = ggplot2::element_rect(fill = NA, color = "black"),
  #       legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
  #       #legend.position = legend.pos,
  #       panel.grid.major = element_blank(),
  #       axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 12),
  #       legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12),
  #       legend.position = "bottom", plot.title = element_text(size = 18),
  #       plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> avg.perc_rast
  #   
  # ggsave(plot = avg.perc_rast, "./Figures/EFH Legal Male/fall.lm.avg.perc_rast.png", 
  #        height = 10, width = 10, units = "in")
  #   
  
  # Plot 3: EFH percentile avg habitat across for last five years ------------------------------------------------
  plot_yrs <- predict_yr[(length(predict_yr)-4):length(predict_yr)]
  
  # Set up list to store maps
  spatpred_list <- list()
  perc.dat <- data.frame()
  
  for(ii in 1:length(plot_yrs)){
    spatpred_df_yr <-  spatpred_df %>%
      filter(year %in% plot_yrs[ii]) %>%
      group_by(x, y) %>%
      reframe(catch_pp = mean(catch_pp))
    
    # Find plotting breaks
    quantiles = c(.05, .25, .5, .75)
    quants<-sort(unique(c(0,quantiles,1)))
    
    threshold <- 0.0513
    
    sample <- stats::na.omit(spatpred_df_yr$catch_pp)
    sample[sample <= threshold] <- NA
    perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
    perc.breaks[1]<-0
    perc.breaks[length(perc.breaks)]<-Inf
    
    # Make raster again
    spatpred_df_yr %>%
      rast() %>%
      raster() -> spatpred_yr
    
    # Set crs
    crs(spatpred_yr) <- map.crs
    
    # Cut the prediction map by the breaks
    perc.map <- raster::cut(spatpred_yr, breaks = perc.breaks)
    
    # set up the factor maps
    perc.vals <- raster::getValues(perc.map)
    perc.vals[perc.vals == 1] <- NA
    
    # convert the raster to polygons
    percpoly0 <- stars::st_as_stars(perc.map)
    percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
    percpoly2 <- percpoly[percpoly$layer != 1, ]
    
    # we'll need a new outline
    perc.dummy.raster <- raster::raster(perc.map)
    perc.vals2 <- is.na(perc.vals) == F
    perc.dummy.raster <- raster::setValues(perc.dummy.raster, values = perc.vals2)
    
    percdummy0 <- stars::st_as_stars(perc.dummy.raster)
    percdummy <- sf::st_cast(sf::st_as_sf(percdummy0, merge = TRUE))
    percdummy2 <- sf::st_transform(percdummy, sf::st_crs(map.crs))
    
    # Dropping the smallest areas
    percdummy.poly <- sf::st_cast(percdummy2, "POLYGON")
    areas <- sf::st_area(percdummy.poly)
    
    outside <- order(areas, decreasing = T)[1]
    toosmall <- which(as.numeric(areas) < 10^8)
    
    perc.x <- percdummy2$layer[-c(outside, toosmall)]
    perc.y <- percdummy2$geometry[-c(outside, toosmall)]
    percdummy3 <- sf::st_sf(perc.x, perc.y) %>%
      vect() %>%
      crop(BB_strata)
    
    # Calculate convex hulls for response data sampling distribution by year
    resp_data %>%
      filter(year == plot_yrs[ii], catch_pp > 0) %>%
      sf::st_as_sf(coords = c("longitude", "latitude"), crs = in.crs) %>%
      sf::st_transform(sf::st_crs(map.crs)) %>%
      vect(.) %>%
      mask(., BB_strata)%>%
      sf::st_as_sf() -> pres_hull_df
    
    # resp_data %>%
    #   filter(year == plot_yrs[ii], catch_pp == 0) %>%
    #   sf::st_as_sf(coords = c("longitude", "latitude"), crs = in.crs) %>%
    #   sf::st_transform(sf::st_crs(map.crs)) %>%
    #   vect() %>%
    #   extract(BB_strata)%>%
    #   sf::st_as_sf() -> abs_hull_df
    
    pres_hull<- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(pres_hull_df))), dist = 15000), dTolerance = 5000)
    #abs_hull<- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(abs_hull_df))), dist = 15000), dTolerance = 5000)
    
    
    # Set up plot boundary
    plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                        x = c(-167.5, -158)) # plot boundary unprojected
    
    plot.boundary <- plot.boundary.untrans %>%
      sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
      sf::st_transform(crs = map.crs) %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      dplyr::rename(x = X, y = Y) # plot boundary projected
    
    # Set up title
    if(seas == "F"){
      title = "Fall Red King Crab Legal Male Encounter Probability"
    } else{
      title = "Winter Red King Crab Legal Male Encounter Probability"
    }
    
    # Set up year label size
    size = ifelse(plot_yrs[ii] < max(plot_yrs), 5, 10)
    lw = ifelse(plot_yrs[ii] < max(plot_yrs), 1, 1.5)
    #col1 <- ifelse(plot_yrs[ii] < max(plot_yrs), NA, "red")
    col1 <- "red"
    #col2 <- ifelse(plot_yrs[ii] < max(plot_yrs), NA, "purple")
    col2 <- "purple"
    
    # Set up year label location
    year_untrans <- data.frame(lab = plot_yrs[ii], x = -158.3, y = 55.3)
    
    year_lab <- year_untrans %>%
      sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
      sf::st_transform(crs = map.crs) %>%
      cbind(st_coordinates(.)) %>%
      as.data.frame()
    
    
    # Map
    ggplot2::ggplot() +
      #ggplot2::geom_sf(data = survey.sf, fill = "grey95")+
      ggplot2::geom_sf(data = percpoly2, ggplot2::aes(fill = as.factor(layer)), col = NA) +
      ggplot2::geom_sf(data = st_as_sf(percdummy3),fill=NA, size = .3) +
      ggplot2::geom_sf(data = st_as_sf(area512),
                       fill = NA,
                       color = col2,
                       linewidth = lw)+
      ggplot2::geom_sf(data = pres_hull, 
                       fill = NA, 
                       color = alpha("white", 0.85),
                       linewidth = 1)+
      # ggplot2::geom_sf(data = abs_hull, 
      #                  fill = NA, 
      #                  color = alpha("white", 0.8),
      #                  linewidth = 1,
      #                  linetype = "dashed")+
      ggplot2::geom_sf(data = st_as_sf(BB_strata),
                       fill = NA,
                       color = "black",
                       linewidth = 2)+
      ggplot2::geom_sf(data = region_layers$akland,
                       fill = "grey70",
                       color = "black")+
      ggplot2::geom_sf(data = st_as_sf(RKCSA),
                       fill = NA,
                       color = col1,
                       linewidth = lw)+
      geom_text(data = year_lab, aes(x=X, y=Y, label= lab), fontface = "bold", size=size) +
      #ggtitle(title)+
      ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                       fill = NA,
                       color = col1,
                       linewidth = lw)+
      coord_sf(xlim = plot.boundary$x,
               ylim = plot.boundary$y)+
      viridis::scale_fill_viridis(discrete = T, name = "Percentiles", labels = c("95%", "75%", "50%", "25%")) +
      ggplot2::theme_bw() +
      ggplot2:: theme(
        panel.border = ggplot2::element_rect(color = "black", fill = NA),
        panel.background = ggplot2::element_rect(fill = NA, color = "black"),
        legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
        #legend.position = legend.pos,
        panel.grid.major = element_blank(),
        axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10),
        legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11),
        legend.position = "bottom", plot.title = element_text(size = 18),
        plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> spatpred_list[[ii]]
    
    #ggsave(spatpred_list[ii], paste0("./Figures/EFH Legal Male/GIF/lm.perc.", plot_yrs[ii], ".png",  height = 10, width = 10, units = "in")
    
    
    # Calculate percent in area 512 and RKCSA by year
    spatpred_df_yr %>%
      st_as_sf(., coords = c("x", "y"), crs = map.crs)  %>%
      vect() -> sp.dat
    
    
    sum(sp.dat$catch_pp) -> total
    
    sp.dat %>%
      mask(., area512) -> a512tot
    
    sp.dat %>%
      mask(., RKCSA_sub) -> RKCSAtot
    
    dat <- data.frame(total = total, a512perc = sum(a512tot$catch_pp)/total, RKCSAperc=sum(RKCSAtot$catch_pp)/total, year = plot_yrs[ii])
    
    perc.dat <- rbind(perc.dat, dat)
    
  }
  
  
  
  # Arrange plots
  ggarrange(ggarrange(spatpred_list[[1]] + rremove("axis.text") + rremove("ticks"),
                      spatpred_list[[2]] + rremove("axis.text") + rremove("ticks"),
                      spatpred_list[[3]] + rremove("axis.text") + rremove("ticks"),
                      spatpred_list[[4]] + rremove("axis.text") + rremove("ticks"),
                      nrow=2, ncol=2, legend = "none"),
            spatpred_list[[5]], heights = c(1,1.15), nrow= 2) %>%
    annotate_figure(top=text_grob(title, size=20, face="bold")) -> yr.perc_rast
  
  # Save
  ggsave(plot = yr.perc_rast, "./Figures/fall.lm.yr.perc_rast.png", height=10, width=7.1, units="in")
  
  
  return(list(pred_raster = pred_raster, spatpred_df = spatpred_df, pred_plot = pred_plot, percpoly2 = percpoly2, percdummy3=percdummy3,
               perc.dat = perc.dat))
  
}

  # Summary dot plot of PA and top 10 along with csv
  MakeDotPlot <- function(resp_data, seas, predict_yr){
  
  presence_df <- data.frame()
  absence_df <- data.frame()
  highdensity_df <- data.frame()
  
  for (ii in 1:length(predict_yr)){
    # Filter response data by season and plotting years
    resp_data %>%
      filter(season == seas, year == predict_yr[ii]) -> resp_data2
    
    # Specify high density quantile
    hd <- stats::quantile(resp_data2 %>% filter(catch_pp >0) %>% dplyr::pull(catch_pp), .9)
    
    # Specify high density, presence, and absence dfs
    presence = resp_data2[resp_data2[, "catch_pp"] > 0, ]
    absence = resp_data2[resp_data2[, "catch_pp"] == 0, ]
    highdensity = resp_data2[resp_data2[, "catch_pp"] >= hd, ]
    
    
    rbind(presence_df, presence) -> presence_df
    rbind(absence_df, absence) -> absence_df
    rbind(highdensity_df, highdensity) -> highdensity_df
    
  }
  
  # Transform data frame crs', make into sf objects
  presence_df %>%
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = in.crs) %>%
    sf::st_transform(sf::st_crs(map.crs)) %>%
    vect(.) %>%
    mask(., BB_strata)%>%
    sf::st_as_sf() -> pres_df
  
  absence_df %>%
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = in.crs) %>%
    sf::st_transform(sf::st_crs(map.crs)) %>%
    vect(.) %>%
    mask(., BB_strata) %>%
    sf::st_as_sf()-> abs_df
  
  highdensity_df %>%
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = in.crs) %>%
    sf::st_transform(sf::st_crs(map.crs)) %>%
    vect(.) %>%
    mask(., BB_strata) %>%
    sf::st_as_sf()-> hd_df
  
  # Set up plot boundary
  plot.boundary.untrans <- data.frame(y = c(54.25, 59.25), 
                                      x = c(-167.5, -158)) # plot boundary unprojected
  
  plot.boundary <- plot.boundary.untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y) # plot boundary projected
  
  # Specify dot features
  abs.name = "absent"
  pres.name = "present"
  hd.name = "top 10%"
  abs.col = "#4B0055" 
  pres.col = "#009B95" 
  hd.col = "#FDE333"
  abs.shape = 16
  pres.shape = 1
  hd.shape = 16
  abs.size = 2
  pres.size = 2
  hd.size = 2
  pres.fac <- 2
  abs.fac <- 1
  hd.fac <- 3
  
  # Now go through and set up the dot locations and add them to legend
  if (is.data.frame(abs_df)) {
    leg.name <- abs.name
    leg.col <- abs.col
    leg.shape <- abs.shape
    leg.size <- abs.size
    abs.fac <- 1
  } else {
    abs.fac <- 0
  }
  
  if (is.data.frame(pres_df)) {
    leg.name <- c(leg.name, pres.name)
    leg.col <- c(leg.col, pres.col)
    leg.shape <- c(leg.shape, pres.shape)
    leg.size <- c(leg.size, pres.size)
    pres.fac <- abs.fac + 1
  }
  
  if (is.data.frame(hd_df)) {
    
    leg.name <- c(leg.name, hd.name)
    leg.col <- c(leg.col, hd.col)
    leg.shape <- c(leg.shape, hd.shape)
    leg.size <- c(leg.size, hd.size)
    hd.fac <- pres.fac + 1
  }
  
  rbind(abs_df %>% mutate(type = "absent"),
        pres_df %>% mutate(type = "present"),
        hd_df %>% mutate(type = "hd")) -> PA_df
  
  # Map
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = PA_df %>% filter(type == "present"), 
                     size = pres.size, ggplot2::aes(color = factor(pres.fac)), shape = pres.shape, stroke = .8)+
    ggplot2::geom_sf(data = PA_df %>% filter(type == "hd"), size = hd.size, 
                     shape = hd.shape, ggplot2::aes(color = factor(hd.fac)))+
    ggplot2::geom_sf(data = PA_df %>% filter(type == "absent"), 
                     alpha = .25, size = abs.size, shape = abs.shape, ggplot2::aes(color = factor(abs.fac)))+
    ggplot2::geom_sf(data = st_as_sf(area512),
                     fill = NA,
                     color = "purple",
                     linewidth = 1.5)+
    ggplot2::geom_sf(data = st_as_sf(BB_strata), 
                     fill = NA, 
                     color = "black",
                     linewidth = 2)+
    ggplot2::geom_sf(data = region_layers$akland, 
                     fill = "grey70", 
                     color = "black")+
    ggplot2::geom_sf(data = st_as_sf(RKCSA),
                     fill = NA,
                     color = "red",
                     linewidth = 1.5)+
    #labs(title = "Fall Directed Fishery Sampling Distribution")+
    ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                     fill = NA,
                     color = "red",
                     linewidth = 1.5)+
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    ggplot2::scale_color_manual(name = NULL, values = c("#009B95", "#FDE333", "#4B0055"), labels = c("present", "top 10%", "absent")) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = c(1, 16, 16), size = leg.size)))+
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      plot.title = element_text(size = 24),
      plot.subtitle = element_text(size = 14),
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 20), legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> sum_dotplot
  
  ggsave(plot = sum_dotplot, "./Figures/fall.lm.dotplot.png", 
         height = 10, width = 10, units = "in")
  
  
  # Summary table
  right_join(
    PA_df %>%
      group_by(year, type) %>%
      reframe(n = n()),
    PA_df %>%
      group_by(year) %>%
      reframe(Total = n()), by = "year") -> sum_table
  
  return(list(sum_dotplot = sum_dotplot, sum_table = sum_table))
}

### RUN FUNCTIONS ----------------------------------------------------------------------------------------
  # Set function parameters
  preds <- Fall_lm.preds
  model_b <- lm.F.modelb
  model_p <- lm.F.modelp
  resp_data <- lm_df
  seas <- "F"
  predict_yr <- c(1998:2019, 2021:2023)

  # Run functions
  MakeDotPlot(lm_df, "F", c(1997:2019, 2021:2023)) -> dot_out

  MakePredictionRaster(preds, model_b, model_p, resp_data, seas, predict_yr) -> out


### EVALUATE MODELS ----------------------------------------------------------------------------------
  lm_F_test %>%
    mutate(catch_pp = catch_pp) -> test2
  
  # Calculate AUC -------------------------------------------------------------
  pred.b <- predict.gbm(lm.F.modelb, # fitted model to predict
                        test2, # data to predict to
                        n.trees=model_b$gbm.call$best.trees, # see help
                        type="response") # predict probabilities
  
  obs <- test2$PA 
  
  d <- cbind(obs, pred.b)
  pres <- d[d[,1]==1, 2]
  abs <- d[d[,1]==0, 2]
  e <- dismo::evaluate(p=pres, a=abs)
  AUC <- e@auc
  

  # Calculate RMSE ------------------------------------------------------------
  pred.p <- predict.gbm(lm.F.modelp, # fitted model to predict
                        test2 %>% filter(catch_pp>0), # data to predict to
                        n.trees=model_p$gbm.call$best.trees, # see help
                        type="response") # predict probabilities
  
  obs <- test2$catch_pp[which(test2$catch_pp >0)]
  pred <- pred.p
  
  RMSE <- sqrt(mean((obs-pred.p)^2))
  
  # Calculate rho -------------------------------------------------------------
  d <- as.data.frame(cbind(obs, pred))
  
  cor.test(d$obs, d$pred, method = "spearman") -> cor_out
  
  vals.ranked <- data.frame(cbind(rank(d$obs, ties.method = 'average'),
                                  rank(d$pred, ties.method = 'average')))
  
  colnames(vals.ranked) <- c('obs', 'preds')
  rho <- cov(vals.ranked) / (sd(vals.ranked$obs) * sd(vals.ranked$preds))
  
  # Calculate PDE --------------------------------------------------------------
  PDE <- (model_p$self.statistics$null-model_p$self.statistics$resid)/model_p$self.statistics$null
  
  # make model evaluation summary data frame
  eval_df <- data.frame(RMSE = RMSE, AUC = AUC, rho = cor_out$estimate, PDE = PDE)
  
  # Var inf plot -------------------------------------------------------------
  summary(model_b) -> sum_b
  summary(model_p) -> sum_p
  
  rbind(sum_b, sum_p) %>%
    group_by(var) %>%
    dplyr::summarise(sum_inf = mean(rel.inf),
                     se = (sd(rel.inf))/sqrt(2)) -> inf_sum
  
  top <- inf_sum[order(inf_sum$sum_inf, decreasing = TRUE), ] %>%
    slice(1:6) %>%
    mutate(var = case_when((var == "Depth")~ "Depth",
                           (var == "Jan_Feb_BT") ~ "Jan/Feb BT",
                           (var == "Tidemax") ~ "Tidal maximum",
                           (var == "Jul_Aug_BT") ~ "Jul/Aug BT",
                           (var == "SAP_Count") ~ "BBRKC survey CPUE",
                           (var == "Sep_Oct_BT") ~ "Sep/Oct BT"))
  
  xlabs = c("Depth", "Maximum tidal\ncurrent", "July/August\nbottom temperature", "January/February\nbottom temperature",  
            "September/October\nbottom temperature", "BBRKC survey\nCPUE") 
  
  
  ggplot(top, aes(reorder(var, -sum_inf), sum_inf)) +
    geom_bar(stat = "identity", color = "black", fill = c("#0E3F5C","#236474","#43898C","#6CAFA3","#9CD6BA","#D1FBD4")) +
    theme_bw() +
    geom_errorbar(aes(ymin=sum_inf - se, ymax=sum_inf+se),size=0.8, width=0.3)+
    xlab(element_blank()) +
    ylab("Relative influence (%)") +
    scale_x_discrete(labels = xlabs)+
    #annotate(geom = "text", label = paste("Spearman's \u03C1 = ", round(cor_out$estimate, 2)), 
    annotate(geom = "text", label = paste("Spearman's \u03C1 = ", round(cor_out$estimate, 2)), 
             x = 3.5, y = 30, fontface = "bold", color = "#4B0055", size = 7)+
    annotate(geom = "text", label = paste("AUC = ", round(AUC, 2)), 
             x = 3.5, y = 26, fontface = "bold", color = "#2A6D7A", size = 7)+
    annotate(geom = "text", label = paste("PDE = ", round(PDE, 2)), 
             x = 3.5, y = 28, fontface = "bold", color = "#4B0055", size = 7)+
    theme(axis.title = element_text(size = 20), axis.text.y = element_text(size = 17),
          #axis.text.x=element_text(size = 17)) -> var_inf_plot
          axis.text.x=element_text(angle=40, hjust=1, size = 17)) -> var_inf_plot
  
  ggsave(plot = var_inf_plot, paste0("./Figures/Fall.LM.varinfplot.png"), 
         height=8, width=10, units="in")
  
  
  # Variable importance information for top 6 model_p vars -----------------------------
  # Model p
  summary(model_p) %>%
    slice(1:6) -> top_p
  
  top_p %>%
    dplyr::select(var) %>%
    pull() %>%
    map_df(~plot(model_p, .x, return.grid = TRUE)) %>%
    reshape::melt(id.vars = "y", variable_name = "var") %>%
    filter(is.na(value) == "FALSE") %>% 
    right_join(., top_p) %>%
    dplyr::rename(Fitted = y, Value = value) -> top_p_df
  
  unique(top_p_df$var) -> top_p_names
  
  top_p_names2 <- c("July/August bottom temperature", "November/December north current", "Sediment", "Slope", 
                    "Depth", "BBRKC survey CPUE")
  
  data.frame(names = top_p_names2, rel.inf = top_p$rel.inf) %>%
    mutate(lab = paste0(names, " (", round(rel.inf, 1), "%)")) -> p_labs
  
  ggplot(top_p_df, aes(x = Value, y = Fitted)) +
    geom_line(linewidth = 1, color = "#4B0055") +
    facet_wrap(~factor(var, levels = top_p_names, labels = p_labs$lab),  labeller = label_wrap_gen(width=25), 
               scales = "free_x", ncol = 3, nrow = 3)+
    ylab("Fitted response")+
    xlab("Observed value")+
    theme_bw() -> p_response
  
  ggsave(plot = p_response, paste0("./Figures/Covariate.Fall.p.response.png"), 
         height=5, width=7, units="in")
  
  # Model b
  summary(model_b) %>%
    slice(1:6) -> top_b
  
  top_b %>%
    dplyr::select(var) %>%
    pull() %>%
    map_df(~plot(model_b, .x, return.grid = TRUE)) %>%
    reshape::melt(id.vars = "y", variable_name = "var") %>%
    filter(is.na(value) == "FALSE") %>% 
    right_join(., top_b) %>%
    dplyr::rename(Fitted = y, Value = value) -> top_b_df
  
  unique(top_b_df$var) -> top_b_names
  
  top_b_names2 <- c("Depth", "Maximum tidal current", "September/October bottom temperature",
                    "January/February bottom temperature", 
                    "May/June north current", "BBRKC survey CPUE")
  
  data.frame(names = top_b_names2, rel.inf = top_b$rel.inf) %>%
    mutate(lab = paste0(names, " (", round(rel.inf, 1), "%)")) -> b_labs
  
  ggplot(top_b_df, aes(x = Value, y = Fitted)) +
    geom_line(linewidth = 1, color = "#009B95") +
    ylab("Fitted response")+
    xlab("Observed value")+
    facet_wrap(~factor(var, levels = top_b_names, labels = b_labs$lab),labeller = label_wrap_gen(width=25),
               scales = "free_x", ncol = 3, nrow = 3)+
    theme_bw() -> b_response
  
  ggsave(plot = b_response, paste0("./Figures/Covariate.Fall.b.response.png"), 
         height=5, width=7, units="in")

### CATCH BY YEAR ---------------------------------------------------------------------------------------
yrs <- expand.grid(predict_year = c(1997:2019, 2021:2023),
                   source = unique(lm_df$source))

ts_dat <- lm_df %>%
  filter(year >= 1997, season == "F") %>%
  rename(predict_year = year)%>%
  group_by(predict_year, source) %>%
  reframe(N = n(),
          var = ((var(catch_pp))/N),
          sd = sqrt(var),
          ci = 1.96 * sd,
          total = (sum(catch_pp)),
          avg = mean(catch_pp)) %>%
  right_join(., yrs)

ggplot(ts_dat %>% filter(predict_year < 2023), mapping = aes(x = predict_year, y = avg, color = source)) +
  geom_ribbon(ts_dat, mapping =aes(ymin = avg - ci, ymax = avg+ci, fill= source, color = NULL),
              alpha = 0.4)+
  geom_line(aes(color = source), linewidth = 1, na.rm = FALSE) +
  geom_point(ts_dat %>% filter(predict_year == 2023),
             mapping = aes(x=predict_year, y = avg, color = source),
             size = 2.5)+
  geom_errorbar(ts_dat %>% filter(predict_year == 2023),
                mapping =aes(ymin = avg - ci, ymax = avg+ci, color = source), linewidth = 1)+
  theme_bw() +
  xlab("Year") +
  ylab("Catch per pot") +
  scale_x_continuous(breaks = seq(min(ts_dat$predict_year), max(ts_dat$predict_year), by = 2),
                     labels= seq(min(ts_dat$predict_year), max(ts_dat$predict_year), by = 2))+
  #ggtitle(paste(mat_sex, "bycatch"))+
  #scale_y_continuous(breaks = breaks)+
  scale_color_manual(values = c("#4B0055", "#FDE333", "#009B95"),
                     name = "",
                     labels = c("Logbook", "Observer - bycatch", "Observer - BBRKC"))+
  scale_fill_manual(values = c("#4B0055", "#FDE333", "#009B95"),
                    name = "",
                    labels = c("Logbook", "Observer - bycatch", "Observer - BBRKC"))+
  ggplot2::theme(
    panel.border = ggplot2::element_rect(color = "black", fill = NA),
    panel.background = ggplot2::element_rect(fill = NA, color = "black"),
    legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
    legend.position = "bottom",
    axis.text = ggplot2::element_text(size = 18),
    axis.title = ggplot2::element_text(size = 20),
    legend.text = ggplot2::element_text(size = 20), legend.title = ggplot2::element_text(size = 12),
    plot.background = ggplot2::element_rect(fill = "white", color = "white"))  -> catch_ppts

ggsave(plot =catch_ppts, paste0("./Figures/Fall.LM.catchppts.png"), 
       height=8, width=10, units="in")

### PREDICTOR MAPS -----------------------------------------------------------------------------
  #Subset present and past predictors
  preds <- terra::subset(Fall_lm.preds, grep(2023, names(Fall_lm.preds), value = T))
  predspst <- terra::subset(Fall_lm.preds, grep(2022, names(Fall_lm.preds), value = T))

  # Dropping May/Jun BT, Mar/Apr BT, and Jul/Aug E current for collinearity
  preds2 <- c(preds$`Jan/Feb mean.temp 2023`, preds$`Jul/Aug mean.temp 2023`, predspst$`Sep/Oct mean.temp 2022`, predspst$`Nov/Dec mean.temp 2022`,
              preds$`Jan_Feb Ice 2023`, preds$`Mar_Apr Ice 2023`,
              preds$`Jan/Feb mean.uEast 2023`, preds$`Mar/Apr mean.uEast 2023`, preds$`May/Jun mean.uEast 2023`, predspst$`Sep/Oct mean.uEast 2022`, predspst$`Nov/Dec mean.uEast 2022`,
              preds$`Jan/Feb mean.vNorth 2023`, preds$`Mar/Apr mean.vNorth 2023`, preds$`May/Jun mean.vNorth 2023`, preds$`Jul/Aug mean.vNorth 2023`, predspst$`Sep/Oct mean.vNorth 2022`, predspst$`Nov/Dec mean.vNorth 2022`,
              preds$`Tidemax 2023`, preds$`Sed 2023`, preds$`Depth 2023`, preds$`Slope 2023`,
              preds$GFbycatch.F2023, preds$`SAP Count 2023`)    
  
  names(preds2) <- c("Jan/Feb bottom temperature", "Jul/Aug bottom temperature", "Sep/Oct bottom temperature", "Nov/Dec bottom temperature",
                     "Jan/Feb Ice", "Mar/Apr Ice", 
                     "Jan/Feb east current", "Mar/Apr east current", "May/Jun east current", "Sep/Oct east current", "Nov/Dec east current",
                     "Jan/Feb north current", "Mar/Apr north current", "May/Jun north current", "Jul/Aug north current", "Sep/Oct north current", "Nov/Dec north current",
                     "Maximum tidal current", "Sediment", "Depth", "Slope",
                     "Groundfish fishery BBRKC bycatch CPUE", "BBRKC survey CPUE")
  
  #Set labels for plotting
  leg.labs <- c("°C", "°C",  "°C","°C", 
                "% cover", "% cover", 
                "m/sec", "m/sec", "m/sec", "m/sec", "m/sec",
                "m/sec", "m/sec", "m/sec", "m/sec", "m/sec", "m/sec",
                "cm/sec", "phi", "Meters", "°",
                "", "") 
  
  
  main.labs <- names(preds2)


  # Specify plot boundaries
  plot.boundary.untrans <- data.frame(y = c(54.25, 59.5), 
                                      x = c(-167.5, -158)) # plot boundary unprojected
  
  plot.boundary <- plot.boundary.untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y) # plot boundary projected
  
  #Set up empty list to hold plots
  pred_list <- list()
  
  #Run for loop to generate predictor plots
  for(ii in 1:length(names(preds2))){
    
    if(leg.labs[ii] != ""){
      ll = leg.labs[ii]
    } else{
      ll = NULL
    }
    
    BB_strata_buff <- buffer(BB_strata, width = -2000)
    
    pp <- terra::subset(preds2, grep(names(preds2)[[ii]], names(preds2), value = T)) %>%
      mask(BB_strata_buff)
    
    pred_df <- data.frame()
    
    pred_df <- cbind(crds(pp), as.data.frame(pp))
    
    colnames(pred_df)[3] <- "value"
    
    
    # Map
    ggplot2::ggplot() +
      #ggplot2::geom_sf(data = survey.sf, fill = "grey95")+
      ggplot2::geom_tile(data = pred_df, aes(x = x, y = y, fill = value))+
      ggplot2::geom_sf(data = st_as_sf(BB_strata), 
                       fill = NA, 
                       color = "black",
                       linewidth = 1)+
      ggplot2::geom_sf(data = region_layers$akland, 
                       fill = "grey70", 
                       color = "black")+
      geom_shadowtext((data.frame(x = -163.6, y = 59.2, lab = names(pp)) %>%
                         sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
                         sf::st_transform(crs = map.crs) %>%
                         cbind(sf::st_coordinates(.), as.data.frame(.))), 
                      mapping = aes(x = X, y = Y, label = lab), color = "black", 
                      bg.color = "lightgrey", size = 2, fontface = "bold") +
      coord_sf(xlim = plot.boundary$x,
               ylim = plot.boundary$y)+
      scale_fill_viridis_c(name = ll)+
      ggplot2::theme_bw() +
      ggplot2:: theme(
        panel.border = ggplot2::element_rect(color = "black", fill = NA),
        panel.background = ggplot2::element_rect(fill = NA, color = "black"),
        #legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
        legend.key.height = unit(0.12, "cm"),
        legend.key.width = unit(0.08, 'cm'),
        panel.grid.major = element_blank(),
        axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 7.5),
        legend.text = ggplot2::element_text(size = 5), legend.title = element_text(size = 5.5),
        legend.position = c(0.87, 0.2), plot.title = element_blank(),
        plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> pred_list[[ii]]
    
  }

  #Arrange plots
  ggarrange(pred_list[[1]], pred_list[[2]], pred_list[[3]], pred_list[[4]],
            pred_list[[5]], pred_list[[6]], pred_list[[7]], pred_list[[8]],
            pred_list[[9]], pred_list[[10]], pred_list[[11]], pred_list[[12]],
            pred_list[[13]], pred_list[[14]], pred_list[[15]], pred_list[[16]], pred_list[[17]],
            pred_list[[18]], pred_list[[19]], pred_list[[20]], pred_list[[21]], pred_list[[22]], 
            pred_list[[23]],
            nrow = 6, ncol = 4)
  
  ggsave("./Figures/Fall.LM.predictor_suite.png", height=12.5, width=13.5, units="in")
  

### PREDICTING WITH MIN/MAX TEMP AND ICE RASTERS --------------------------------------------------
  # Subset preds by hot and cold stanzas and set names
  predict_yr <- 2023 # this is just a dummy variable to generate the plots. Not actually only data from 2023
  
  preds <- rast("./Data/Fall_lm.preds.mm.tif")
  
  preds2 <- subset(preds, grep(predict_yr, names(preds)))
  preds3 <- preds2[[c(-1:-6, -8, -10, -12:-19)]]
  
  preds.hot <- subset(preds3, grep("BT min |currentE min |currentN min |SAP Count min |ice max ", 
                                   names(preds3), invert = TRUE))
  
  names(preds.hot)<- c("Sed", "Depth", "Slope", "Tidemax",
                       "Jan_Feb_BT", "Mar_Apr_BT", "May_Jun_BT", "Jul_Aug_BT", "Sep_Oct_BT", "Nov_Dec_BT", 
                       "Jan_Feb_currentE", "Mar_Apr_currentE", "May_Jun_currentE", "Jul_Aug_currentE", "Sep_Oct_currentE", "Nov_Dec_currentE",
                       "Jan_Feb_currentN", "Mar_Apr_currentN", "May_Jun_currentN", "Jul_Aug_currentN", "Sep_Oct_currentN", "Nov_Dec_currentN", 
                       "Jan_Feb_Ice", "Mar_Apr_Ice", "SAP_Count")
  
  preds.cold <- subset(preds3, grep("BT max |currentE max |currentN max |SAP Count max |ice min ", 
                                    names(preds3), invert = TRUE))

  names(preds.cold)<- c("Sed", "Depth", "Slope", "Tidemax",
                        "Jan_Feb_BT", "Mar_Apr_BT", "May_Jun_BT", "Jul_Aug_BT", "Sep_Oct_BT", "Nov_Dec_BT", 
                        "Jan_Feb_currentE", "Mar_Apr_currentE", "May_Jun_currentE", "Jul_Aug_currentE", "Sep_Oct_currentE", "Nov_Dec_currentE",
                        "Jan_Feb_currentN", "Mar_Apr_currentN", "May_Jun_currentN", "Jul_Aug_currentN", "Sep_Oct_currentN", "Nov_Dec_currentN", 
                        "Jan_Feb_Ice", "Mar_Apr_Ice", "SAP_Count")

# Write function to generate spatial predictions by for hot/cold stanza
  hotcold <- function(preds, resp_data, model_b, model_p, predict_yr, temp){
  # Create data frame for wind (since no spatial layer) to feed models by year
  data2 <- resp_data %>%
    filter(year == predict_yr) %>%
    mutate(fishery = as.factor(fishery))
  
  
  # Generate spatial model predictions for Bristol Bay management area extent using raster predictors
    #PA
    BB_strata_buff <- buffer(BB_strata, width = -1000)
    
    spatpred_b <- suppressWarnings(predict(preds, # raster stack 
                                           model_b, # fitted model
                                           n.trees=model_b$gbm.call$best.trees, # see help
                                           #factors = f,
                                           type="response", 
                                           ext = panel_extent) %>%
                                     #mask(BB_strata, touches = FALSE) %>%
                                     mask(preds$Jan_Feb_Ice) %>%
                                     mask(BB_strata_buff, touches = FALSE) %>%
                                     na.omit())
    
    #Abund
    spatpred_p <- suppressWarnings(predict(preds, # raster stack 
                                           model_p, # fitted model
                                           n.trees=model_p$gbm.call$best.trees, # see help
                                           #factors = ff,
                                           type="response", 
                                           ext = panel_extent) %>%
                                     #mask(BB_strata, touches = FALSE)  %>%
                                     mask(preds$Jan_Feb_Ice) %>%
                                     mask(BB_strata_buff, touches = FALSE) %>%
                                     na.omit())
    
    # Multiply binomial raster by abundance raster (part of delta model framework)
    spatpred <- spatpred_b * spatpred_p
  
  
  # Create data frame of spatial predictions (for pretty mapping)
  spatpred_df <- cbind(crds(spatpred$lyr1), as.data.frame(spatpred$lyr1)) %>%
    dplyr::rename("catch_pp" = lyr1) %>%
    mutate(year = predict_yr) 
  
  # Avg habitat 
  spatpred_df2 <- spatpred_df %>%
    group_by(x, y) %>%
    reframe(catch_pp = mean(catch_pp))
  
  # Find plotting breaks
  quantiles = c(.05, .25, .5, .75)
  quants<-sort(unique(c(0,quantiles,1)))
  
  threshold <- 0.0513
  
  sample <- stats::na.omit(spatpred_df2$catch_pp)
  sample[sample <= threshold] <- NA
  perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
  perc.breaks[1]<-0
  perc.breaks[length(perc.breaks)]<-Inf
  
  
  # Make raster again
  spatpred_df2 %>%
    rast() %>%
    raster() -> spatpred
  
  # Set crs
  crs(spatpred) <- map.crs
  
  # Cut the prediction map by the breaks
  perc.map <- raster::cut(spatpred, breaks = perc.breaks)
  
  # set up the factor maps
  perc.vals <- raster::getValues(perc.map)
  perc.vals[perc.vals == 1] <- NA
  
  # convert the raster to polygons
  percpoly0 <- stars::st_as_stars(perc.map)
  percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
  percpoly2 <- percpoly[percpoly$layer != 1, ]
  
  # we'll need a new outline
  perc.dummy.raster <- raster::raster(perc.map)
  perc.vals2 <- is.na(perc.vals) == F
  perc.dummy.raster <- raster::setValues(perc.dummy.raster, values = perc.vals2)
  
  percdummy0 <- stars::st_as_stars(perc.dummy.raster)
  percdummy <- sf::st_cast(sf::st_as_sf(percdummy0, merge = TRUE))
  percdummy2 <- sf::st_transform(percdummy, sf::st_crs(map.crs))
  
  # Dropping the smallest areas
  percdummy.poly <- sf::st_cast(percdummy2, "POLYGON")
  areas <- sf::st_area(percdummy.poly)
  
  outside <- order(areas, decreasing = T)[1]
  toosmall <- which(as.numeric(areas) < 10^8)
  
  perc.x <- percdummy2$layer[-c(outside, toosmall)]
  perc.y <- percdummy2$geometry[-c(outside, toosmall)]
  percdummy3 <- sf::st_sf(perc.x, perc.y) %>%
    vect() %>%
    crop(BB_strata)
  
  # Set up plot boundary
  plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                      x = c(-167.5, -158)) # plot boundary unprojected
  
  plot.boundary <- plot.boundary.untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y) # plot boundary projected
  
  year_untrans <- data.frame(lab = paste0("'", temp,"'"),
                             x = -158.3, y = 55.3)
  
  year_lab <- year_untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    cbind(st_coordinates(.)) %>%
    as.data.frame()
  
  # Set up title
  title = paste0("Fall Red King Crab Legal Male Encounter Probability")
  
  # Map
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = percpoly2, ggplot2::aes(fill = as.factor(layer)), color = NA) +
    ggplot2::geom_sf(data = st_as_sf(percdummy3),fill=NA, size = .3) +
    ggplot2::geom_sf(data = st_as_sf(area512),
                     fill = NA,
                     color = "purple",
                     linewidth = 0.75)+
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 1)+
    ggplot2::geom_sf(data = region_layers$akland, 
                     fill = "grey70", 
                     color = "black")+
    ggplot2::geom_sf(data = st_as_sf(RKCSA),
                     fill = NA,
                     color = "red",
                     linewidth = 0.75)+
    ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                     fill = NA,
                     color = "red",
                     linewidth = 0.75)+
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    geom_text(data = year_lab, aes(x=X, y=Y, label= lab), fontface = "bold", size=4) +
    viridis::scale_fill_viridis(discrete = T, name = "Percentiles", labels = c("95%", "75%", "50%", "25%")) +
    ggplot2::theme_bw() +
    ggplot2:: theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      #legend.position = legend.pos,
      panel.grid.major = element_blank(),
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11),
      legend.position = "bottom", plot.title = element_text(size = 18),
      plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> perc_rast
  
  ggsave(plot = perc_rast, paste0("./Figures/Fall.LM.percrast.", temp, ".png"), width = 7, height = 6, units = "in")
  
  return(list(spatpred_df, spatpred_df2, perc_rast))
}

# Run function for hot stanza
  hotcold(preds.hot, lm_df, lm.F.modelb, lm.F.modelp, predict_yr, "Warm") -> out.hot

# Calculate RKC abundance proportion inside and outside of closure areas in hot years
  out.hot[[2]] %>%
    st_as_sf(., coords = c("x", "y"), crs = map.crs)  %>%
    vect() -> sp.dat
  
  sum(sp.dat$catch_pp) -> total.hot
  
  sp.dat %>%
    mask(., area512) -> a512.hot
  
  sp.dat %>%
    mask(., RKCSA_sub) -> RKCSA.hot
  
  sum(a512.hot$catch_pp)/total.hot -> a512.perc.hot
  sum(RKCSA.hot$catch_pp)/total.hot -> RKCSA.perc.hot

 # Run function for cold stanza
  hotcold(preds.cold, lm_df, lm.F.modelb, lm.F.modelp, predict_yr, "Cold") -> out.cold

 # Calculate RKC abundance proportion inside and outside of closure areas in cold years
  out.cold[[2]] %>%
    st_as_sf(., coords = c("x", "y"), crs = map.crs)  %>%
    vect() -> sp.dat
  
  sum(sp.dat$catch_pp) -> total.cold
  
  sp.dat %>%
    mask(., area512) -> a512.cold
  
  sp.dat %>%
    mask(., RKCSA_sub) -> RKCSA.cold
  
  sum(a512.cold$catch_pp)/total.cold -> a512.perc.cold
  sum(RKCSA.cold$catch_pp)/total.cold -> RKCSA.perc.cold

  # Arrange plots
  ggarrange(out.hot[[3]], out.cold[[3]], heights = c(1,1), nrow= 1, ncol = 2,
            common.legend = TRUE, legend = "bottom") -> yr.perc_rast_hotv.cold.stanza
  
  ggsave(plot = yr.perc_rast_hotv.cold.stanza, "./Figures/F.LM.yrrast.hotvcold.stanza.png", height=3, width=6, units="in")

