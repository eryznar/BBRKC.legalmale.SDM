#PURPOSE:
# To fit and evaluate boosted regression trees to model Bristol Bay red king crab
# legal male distribution in the fall

#AUTHOR:
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)

### LOAD PROCESSING PARAMETERS --------------------------------------------------
source("./Scripts/load.libs.params.R")

read.csv("./Output/lm_F_train.csv") %>%
  dplyr::select(!X) -> train
read.csv("./Output/lm_F_test.csv") %>%
  dplyr::select(!X) -> test


# ### GENERATE TRAINING/TESTING DATA ----------------------------------------------------------------------------------------
# # Read in spatial covariates
# Fall_lm.preds <- rast("./Data/Fall_lm.predsxy.tif")
# 
# # Read in response data
# lm_df <- read.csv("./Data/legalmale_direct.fish.csv") %>%
#   dplyr::select(!X)
# 
# # Filter response data by season, transform into spatial vectors, mask by BB management area
# lm_df %>%
#   filter(season == "F") %>%
#   sf::st_as_sf(coords = c(x = "longitude", y = "latitude"), crs = sf::st_crs(4326)) %>%
#   sf::st_transform(crs = map.crs) %>%
#   vect(.) %>%
#   mask(., BB_strata) -> data_vect
# 
# preds <- subset(Fall_lm.preds, grep("GFbycatch.F", names(Fall_lm.preds), invert = TRUE))
# 
# # Generate for loop to extract covariate values at bycatch presence points by year, store as list of dfs,
# # assign labels
# pst <- c(1997:2019, 2021:2022)
# predict_yr <- c(1998:2019, 2021:2023)
# 
# datalist = list()
# datalist = vector("list", length = length(predict_yr))
# 
# for (ii in 1:length(predict_yr)){
#   
#   preds_past <- subset(preds, grep("Sep/Oct |Nov/Dec ", names(preds)))
#   preds_pres <- subset(preds, grep("Sep/Oct |Nov/Dec ", names(preds), invert = TRUE))
#   
#   
#   dat_pres <- terra::extract(subset(preds_pres, grep(predict_yr[ii], names(preds_pres), value = T)), 
#                              subset(data_vect, data_vect$year == predict_yr[ii]))
#   
#   dat_past <- terra::extract(subset(preds_past, grep(pst[ii], names(preds_past), value = T)), 
#                              subset(data_vect, data_vect$year == predict_yr[ii]))
#   
#   dat <- cbind(dat_past[,-1], dat_pres[,-1],
#                crds(subset(data_vect, data_vect$year == predict_yr[ii])),
#                subset(data_vect[,c("year", "catch_pp", "fishery")], data_vect$year == predict_yr[ii]))
#   
#   labs <- c("Sep_Oct_BT", "Nov_Dec_BT", "Sep_Oct_currentE", "Nov_Dec_currentE", "Sep_Oct_currentN",
#             "Nov_Dec_currentN", 
#             "Jan_Feb_BT", "Mar_Apr_BT", "May_Jun_BT", "Jul_Aug_BT", 
#             "Jan_Feb_Ice", "Mar_Apr_Ice", 
#             "Sed", "SAP_Count", "Depth", "Slope", 
#             "Jan_Feb_currentE", "Mar_Apr_currentE", "May_Jun_currentE", "Jul_Aug_currentE",
#             "Jan_Feb_currentN", "Mar_Apr_currentN", "May_Jun_currentN", "Jul_Aug_currentN",
#             "Tidemax") #raster labs
#   
#   
#   names(dat) <- c(labs, "x", "y", "predict_year", "catch_pp", "fishery")
#   
#   datalist[[ii]] <- dat
#   
# }
# 
# # Bind data by year
# dplyr::bind_rows(datalist) %>%
#   na.omit() %>%
#   mutate(#geographic_position = Latitude * Longitude,
#     #Current = CurrentN*CurrentE,
#     #CurrentSD = CurrentSDE*CurrentSDN,
#     PA = ifelse(catch_pp == 0, 0, 1)) -> sdmData # bind presence dfs
# 
# 
# # Check for variable collinearity, drop collinear vars
# non_covs <- c("x", "y", "predict_year", "catch_pp", "PA", "fishery")
# 
# sdmData[-which(colnames(sdmData) %in% non_covs)] %>%
#   as.data.frame() %>%
#   na.omit() %>%
#   vifstep(th=5) -> vif_drop 
# 
# dropvars <- vif_drop@excluded 
# 
# sdmData <- sdmData[-which(colnames(sdmData) %in% dropvars)]
# 
# # Create training/testing datasets, randomly drawing 10 times
# iter <- seq(1, 10, by =1)
# 
# trainlist = list()
# trainlist = vector("list", length = length(iter))
# 
# testlist = list()
# testlist = vector("list", length = length(iter))
# 
# set.seed(1)
# 
# for (ii in 1:length(iter)){
#   group1 <- kfold(sdmData, k = 5)
#   trainlist[[ii]] <- sdmData[group1 != 1,]
#   testlist[[ii]]  <- sdmData[group1 == 1, ]
#   
#   trainlist[[ii]]$iter <- iter[ii]
#   testlist[[ii]]$iter <- iter[ii]
# }
# 
# train <- dplyr::bind_rows(trainlist)
# test <- dplyr::bind_rows(testlist)
# 
# 
# 
# ### CREATE FUNCTION TO FIT MODELS ON REPLICATE TRAINING/TESTING SPLITS --------------
# model_iter <- function(train, test, seas, sptemp, iteration){
#   # Filter by iteration
#   train2 <- train %>% 
#     filter(iter %in% iteration) %>%
#     mutate(catch_pp = as.integer(catch_pp),
#            fishery = as.factor(fishery))
#   
#   test2 <- test %>% 
#     filter(iter %in% iteration) %>%
#     mutate(catch_pp = as.integer(catch_pp),
#            fishery = as.factor(fishery))
#   
#   # Specify non-predictor variables
#   if(sptemp == "space"){
#     non_covs <- c("predict_year", "catch_pp", "PA", "iter", "fishery")
#   }else{
#     non_covs <- c("catch_pp", "PA", "iter", "fishery")
#   }
#  
#     
#   # Fit binomial model with spatial covariates
#   model_b <- gbm.step(data = train2, 
#                       gbm.x = which(!colnames(train2) %in% non_covs), # columns of predictors
#                       gbm.y = which(colnames(train2) == "PA"), # column of response
#                       family = "bernoulli", # for PA data
#                       tree.complexity = 5, # model interactions, etc?
#                       learning.rate = 0.05, # influence of each tree 
#                       bag.fraction = 0.5) 
#   
#   saveRDS(model_b, paste0("./Models/lm.modelb.", seas,".", iteration, sptemp, ".rda"))
#   
#   ntreesb <- model_b$n.trees
#   
#   model_p <- gbm.step(data = train2 %>% filter(catch_pp >0), 
#                       gbm.x = which(!colnames(train2) %in% non_covs), # columns of predictors 
#                       gbm.y = which(colnames(train2) == "catch_pp"), # column of response
#                       family = "poisson", # for count data
#                       tree.complexity = 5, # model interactions, etc?
#                       learning.rate = 0.05, # influence of each tree 
#                       bag.fraction = 0.5)
#   
#   saveRDS(model_p, paste0("./Models/lm.modelp.", seas,".", iteration, sptemp, ".rda"))
#   
#   ntreesp <- model_p$n.trees
#   
#   # Calculate AUC for binomial
#   pred <- predict.gbm(model_b, # fitted model to predict
#                       test2, # data to predict to
#                       n.trees=model_b$gbm.call$best.trees, # see help
#                       type="response") # predict probabilities
#   
#   obs <- test2$PA 
#   
#   d <- cbind(obs, pred)
#   pres <- d[d[,1]==1, 2]
#   abs <- d[d[,1]==0, 2]
#   e <- dismo::evaluate(p=pres, a=abs)
#   AUC <- e@auc
#   
#   
#   # Calculate RMSE and spearman's for positive model
#   pred <- predict.gbm(model_p, # fitted model to predict
#                       test2 %>% filter(catch_pp>0), # data to predict to
#                       n.trees=model_p$gbm.call$best.trees, # see help
#                       type="response") # predict probabilities
#   
#   obs <- test2$catch_pp[which(test2$catch_pp >0)]
#   
#   RMSE <- sqrt(mean((obs-pred)^2))
#   
#   d <- as.data.frame(cbind(obs, pred))
#   
#   ggplot(d, aes(x = obs, y = pred))+
#     geom_point()+
#     geom_smooth(method = "lm") + 
#     theme_classic()+
#     xlab("Observed vals")+
#     ylab("Predicted vals") -> pred_obs_plot
#   
#   cor.test(d$obs, d$pred, method = "spearman") -> cor_out
#   
#   rho <- cor_out$estimate
#   p <- cor_out$p.value
#   
#   # Bind metrics into a data frame
#   eval_df <- data.frame(ntreesb = ntreesb, ntreesp = ntreesp, AUC = AUC, RMSE = RMSE, rho = rho, p = p, 
#                         iter = iteration, sptemp = sptemp)
#   
#   write.csv(eval_df, paste0("./Models/model_iteration_", sptemp, ".csv"))
#   
#   return(list(eval_df))
# }

# ### FUN FUNCTION TO FIT MODELS ------------------------------------------------------
# 
# 
# # Run model iteration function
# 1:10 %>%
#   map_df(~model_iter(train, test, "F", "space", .x)) -> model_out
# 
# write.csv(model_out, "./Models/model_iteration_space.csv")
# 
# readRDS("./Models/")
# 
# 1:10 %>%
#   map_df(~model_iter(train, test, "F", "spacetime", .x)) -> model_out
# 
# write.csv(model_out, "./Models/model_iteration_spacetime.csv")
# read.csv("./Models/model_iteration_spacetime.csv")

### COMPARE MODEL RESIDUALS OF BEST NON-SPATIOTEMPORAL MODEL and SPATIOTEMPORAL MODELS--------------------------
  # Load best models for each scenario
  best.b <- readRDS("./Models/lm.modelb.F.8.rda") # non-spatiotemporal
  best.p <- readRDS("./Models/lm.modelp.F.8.rda") # non-spatiotemporal
  
  b.space <- readRDS("./Models/lm.modelb.F.6space.rda") # spatial
  p.space <- readRDS("./Models/lm.modelp.F.6space.rda") # spatial
  
  b.spacetime <- readRDS("./Models/lm.modelb.F.8spacetime.rda") # spatiotemporal
  p.spacetime <- readRDS("./Models/lm.modelp.F.8spacetime.rda") # spatiotemporal
  
  # Set up plot boundary
  plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                      x = c(-167.5, -158)) # plot boundary unprojected
  
  plot.boundary <- plot.boundary.untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y) # plot boundary projected
  
  # Residual patterning (non-spatiotemporal)
  residuals(best.b, type = "mle-mvn") -> resid.b
  residuals(best.p, type = "mle-mvn") -> resid.p
  
  cbind(train %>% filter(iter == 6), resid = resid.b) -> rr
  
  # visualize residuals
  ggplot(rr) + 
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 1)+
    geom_point(aes(y = y, x = x, color = resid), size = 1, alpha = 0.5) +
    scale_color_gradient2(midpoint = 0) +
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 1)+
    ggplot2::geom_sf(data = region_layers$akland,
                     fill = "grey70",
                     color = "black")+
    facet_wrap(~predict_year)+
    labs(y = "Latitude",
         x = "Longitude") +
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    ggtitle("Env (AUC=0.95)")+
    theme_bw() + 
    theme(axis.title = element_text(size = 10),
          title = element_text(size = 10),
          axis.text = element_text(size = 8),
          legend.title =element_text(size = 10)) -> plot.1
  
  # Residual patterning (spatial)
  residuals(b.space, type = "mle-mvn") -> resid.b.space
  residuals(p.space, type = "mle-mvn") -> resid.p.space
  
  cbind(train %>% filter(iter == 6), resid = resid.b.space) -> rr
  
  # visualize residuals
  ggplot(rr) + 
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 1)+
    geom_point(aes(y = y, x = x, color = resid), size = 1, alpha = 0.5) +
    scale_color_gradient2(midpoint = 0) +
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 1)+
    ggplot2::geom_sf(data = region_layers$akland,
                     fill = "grey70",
                     color = "black")+
    facet_wrap(~predict_year)+
    labs(y = "Latitude",
         x = "Longitude") +
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    ggtitle("Env + space (AUC=0.95)")+
    theme_bw() + 
    theme(axis.title = element_text(size = 10),
          title = element_text(size = 10),
          axis.text = element_text(size = 8),
          legend.title =element_text(size = 10)) -> plot.2
  
  
  # Residual patterning (spatiotemporal)
  residuals(b.spacetime, type = "mle-mvn") -> resid.b.spacetime
  residuals(p.spacetime, type = "mle-mvn") -> resid.p.spacetime
  
  cbind(train %>% filter(iter == 8), resid = resid.b.spacetime) -> rr
  
  # visualize residuals
  ggplot(rr) + 
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 1)+
    geom_point(aes(y = y, x = x, color = resid), size = 1, alpha = 0.5) +
    scale_color_gradient2(midpoint = 0) +
    facet_wrap(~predict_year)+
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 1)+
    ggplot2::geom_sf(data = region_layers$akland,
                     fill = "grey70",
                     color = "black")+
    facet_wrap(~predict_year)+
    labs(y = "Latitude",
         x = "Longitude") +
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    ggtitle("Env + space + time (AUC=0.95)")+
    theme_bw() + 
    theme(axis.title = element_text(size = 10),
          title = element_text(size = 10),
          axis.text = element_text(size = 8),
          legend.title =element_text(size = 10)) -> plot.3
  
# cowplot::plot_grid(plot.1, plot.2, plot.3, ncol = 1)
#  
# ggsave("./Figures/residual_plots.png", width = 6.5, height = 10, units = "in") 

## Function to calculate and plot residuals -------------------------------------
eval_spatresid <- function(model_b, model_p, type, train, test, preds, iteration, observ, predict_yr){
  # Filter by iteration
  train2 <- train %>% 
    filter(iter %in% iteration) %>%
    mutate(catch_pp = as.integer(catch_pp),
           fishery = as.factor(fishery))
  
  
  test2 <- test %>% 
    filter(iter %in% iteration) %>%
    mutate(catch_pp = as.integer(catch_pp),
           fishery = as.factor(fishery))
  

  # Calculate residuals from delta model ----
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
              "x", "y",
              "Jan_Feb_BT", "Mar_Apr_BT", "May_Jun_BT", "Jul_Aug_BT", 
              "Jan_Feb_Ice", "Mar_Apr_Ice", 
              "Sed", "SAP_Count", "Depth", "Slope", 
              "Jan_Feb_currentE", "Mar_Apr_currentE", "May_Jun_currentE", "Jul_Aug_currentE",
              "Jan_Feb_currentN", "Mar_Apr_currentN", "May_Jun_currentN", "Jul_Aug_currentN",
              "Tidemax") #raster labs
    
    
    
    # Set names of preds
    names(preds2) <- labs
    
    # # Create data frame for wind (since no spatial layer) to feed models by year
    # data2 <- resp_data %>%
    #   filter(year == predict_yr[ii]) %>%
    #   mutate(fishery = as.factor(fishery), year = as.factor(year))
    
    yr <- predict_yr[ii]
    
    yr <- as.data.frame(yr)
    names(yr) <- "predict_year"
    
    # Generate spatial model predictions for Bristol Bay management area extent using raster predictors
    #PA
    BB_strata_buff <- buffer(BB_strata, width = -1000)
    
    if(type == "spacetime"){
      spatpred_b <- suppressWarnings(terra::predict(preds2, # raster stack 
                                                    model_b, # fitted model
                                                    n.trees=model_b$gbm.call$best.trees, # see help
                                                    const = yr,
                                                    type="response", 
                                                    ext = ext(preds2$Jan_Feb_Ice)) %>%
                                       #mask(BB_strata, touches = FALSE) %>%
                                       mask(preds2$Jan_Feb_Ice) %>%
                                       mask(BB_strata_buff, touches = FALSE) %>%
                                       na.omit())
      
      
      #Abund
      spatpred_p <- suppressWarnings(terra::predict(preds2, # raster stack 
                                             model_p, # fitted model
                                             n.trees=model_p$gbm.call$best.trees, # see help
                                             const = yr,
                                             type="response", 
                                             ext = ext(preds2$Jan_Feb_Ice)) %>%
                                       #mask(BB_strata, touches = FALSE)  %>%
                                       mask(preds2$Jan_Feb_Ice) %>%
                                       mask(BB_strata_buff, touches = FALSE) %>%
                                       na.omit())
      
    }else{
      spatpred_b <- suppressWarnings(terra::predict(preds2, # raster stack 
                                                    model_b, # fitted model
                                                    n.trees=model_b$gbm.call$best.trees, # see help
                                                    type="response", 
                                                    ext = ext(preds2$Jan_Feb_Ice)) %>%
                                       #mask(BB_strata, touches = FALSE) %>%
                                       mask(preds2$Jan_Feb_Ice) %>%
                                       mask(BB_strata_buff, touches = FALSE) %>%
                                       na.omit())
      
      
      #Abund
      spatpred_p <- suppressWarnings(terra::predict(preds2, # raster stack 
                                             model_p, # fitted model
                                             n.trees=model_p$gbm.call$best.trees, # see help
                                             type="response", 
                                             ext = ext(preds2$Jan_Feb_Ice)) %>%
                                       #mask(BB_strata, touches = FALSE)  %>%
                                       mask(preds2$Jan_Feb_Ice) %>%
                                       mask(BB_strata_buff, touches = FALSE) %>%
                                       na.omit())
    }
    
    # Multiply binomial raster by abundance raster (part of delta model framework)
    spatpred <- spatpred_b * spatpred_p
    
    thres <- mean(na.omit(values(spatpred_b$lyr1)))
    # 
    # # # Create data frame of spatial predictions (for pretty mapping)
    # out <- cbind(crds(spatpred$lyr1), as.data.frame(spatpred$lyr1)) %>%
    #   dplyr::rename("catch_pp" = lyr1) %>%
    #   mutate(year = predict_yr[ii])
    # 
    # Observations
    obs2 <- terra::subset(observ, grep(predict_yr[ii], names(observ)))
    
    # Calculate residuals
    res.rast <- spatpred-obs2 #(predicted-observed)
    
    res <- cbind(crds(res.rast$lyr1), as.data.frame(res.rast$lyr1)) %>%
      dplyr::rename("residuals" = lyr1) %>%
      mutate(year = predict_yr[ii])
    
    # Add to prediction df
    c(spatpred, res.rast) -> tt
    names(tt) <- c("predicted", "residuals")
    
    out <- cbind(crds(tt$predicted), as.data.frame(tt)) %>%
      mutate(year = predict_yr[1])


    rbind(spatpred_df, out) -> spatpred_df
  }
  
  
  return(spatpred_df = spatpred_df)
}

observ <- rast("./Data/obs_rasters.tif")

# c(1998:2019, 2021:2023) %>%
# purrr::map_df(~eval_spatresid(b.spacetime, p.spacetime, "spacetime", train, test, preds, 8, observ, .x)) -> out.st
# 
# c(1998:2019, 2021:2023) %>%
#   purrr::map_df(~eval_spatresid(best.b, best.p, "env", train, test, preds, 8, observ, .x)) -> out.env

write.csv(out.st, "./Output/resid.df.spacetime.csv")
write.csv(out.env, "./Output/resid.df.env.csv")

# Load data
out.st <- read.csv("./Output/resid.df.spacetime.csv")
out.env <- read.csv("./Output/resid.df.env.csv")

# Function to evaluate models
eval_models <- function(model_b, model_p, train, test, iteration){
  # Filter by iteration
  train2 <- train %>% 
    filter(iter %in% iteration) %>%
    mutate(catch_pp = as.integer(catch_pp),
           fishery = as.factor(fishery))
  
  
  test2 <- test %>% 
    filter(iter %in% iteration) %>%
    mutate(catch_pp = as.integer(catch_pp),
           fishery = as.factor(fishery))
  
  # Calculate AUC for binomial ----
  pred <- predict.gbm(model_b, # fitted model to predict
                      test2, # data to predict to
                      n.trees=model_b$gbm.call$best.trees, # see help
                      type="response") # predict probabilities
  
  obs <- test2$PA 
  
  d <- cbind(obs, pred)
  pres <- d[d[,1]==1, 2]
  abs <- d[d[,1]==0, 2]
  e <- dismo::evaluate(p=pres, a=abs)
  AUC <- e@auc
  
  
  # Calculate RMSE and spearman's for positive model ----
  pred <- predict.gbm(model_p, # fitted model to predict
                      test2 %>% filter(catch_pp>0), # data to predict to
                      n.trees=model_p$gbm.call$best.trees, # see help
                      type="response") # predict probabilities
  
  obs <- test2$catch_pp[which(test2$catch_pp >0)]
  
  RMSE <- sqrt(mean((obs-pred)^2))
  
  d <- as.data.frame(cbind(obs, pred))
  
  cor.test(d$obs, d$pred, method = "spearman") -> cor_out
  
  rho <- cor_out$estimate
  p <- cor_out$p.value
  
  # Calculate PDE for positive ----
  PDE <- (model_p$self.statistics$null-model_p$self.statistics$resid)/model_p$self.statistics$null
  
  # Bind eval metrics ----
  eval_df <- data.frame(AUC = AUC, RMSE = RMSE, rho = rho, PDE = PDE)
  
  return(eval_df = eval_df)
}
 
# Evaluate predictive performance - ENV
eval_models(best.b, best.p, train, test, 8) -> eval.env

eval_models(b.spacetime, p.spacetime, train, test, 8) -> eval.st

# Set up plot boundary
plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                    x = c(-167.5, -158)) # plot boundary unprojected

plot.boundary <- plot.boundary.untrans %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(crs = map.crs) %>%
  sf::st_coordinates() %>%
  as.data.frame() %>%
  dplyr::rename(x = X, y = Y) # plot boundary projected

# Residual patterning (non-spatiotemporal)
rr <- na.omit(out.env)

# visualize residuals
ggplot(rr) + 
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   alpha = 0.4,
                   color = "black",
                   linewidth = 1)+
  geom_point(aes(y = y, x = x, color = residuals),  size = 0.05) +
  scale_color_gradient2(midpoint = 0, high = "#6B6100", low = "#4F53B7", mid = "#F1F1F1") +
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   color = "black",
                   linewidth = 1)+
  ggplot2::geom_sf(data = region_layers$akland,
                   fill = "grey70",
                   color = "black")+
  facet_wrap(~year)+
  labs(y = "Latitude",
       x = "Longitude") +
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y)+
  ggtitle("Env (AUC=0.95, Spearman's rho=0.46, PDE=0.36)")+
  theme_bw() + 
  theme(axis.title = element_text(size = 10),
        title = element_text(size = 10),
        axis.text = element_text(size = 6),
        legend.title =element_text(size = 10)) -> plot.1

# visualize residuals
ggplot(rr) + 
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   alpha = 0.4,
                   color = "black",
                   linewidth = 1)+
  geom_point(aes(y = y, x = x, color = residuals), size = 0.05) +
  scale_color_gradient2(midpoint = 0, high = "#6B6100", low = "#4F53B7", mid = "#F1F1F1") +
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   color = "black",
                   linewidth = 1)+
  ggplot2::geom_sf(data = region_layers$akland,
                   fill = "grey70",
                   color = "black")+
  #facet_wrap(~year)+
  labs(y = "Latitude",
       x = "Longitude") +
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y,
           expand = TRUE)+
  ggtitle("Env (AUC=0.95, Spearman's rho=0.46, PDE=0.36)")+
  theme_bw() + 
  theme(axis.title = element_text(size = 10),
        title = element_text(size = 10),
        axis.text = element_text(size = 6),
        legend.title =element_text(size = 10)) -> plot.2


# Residual patterning (non-spatiotemporal)
rr <- na.omit(out.st[,-c(1:6)])
colnames(rr) <- c("x", "y", "predicted", "residuals", "year")


# visualize residuals
ggplot(rr) + 
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   alpha = 0.4,
                   color = "black",
                   linewidth = 1)+
  geom_point(aes(y = y, x = x, color = residuals), size = 0.05) +
  scale_color_gradient2(midpoint = 0, high = "#6B6100", low = "#4F53B7", mid = "#F1F1F1") +
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   color = "black",
                   linewidth = 1)+
  ggplot2::geom_sf(data = region_layers$akland,
                   fill = "grey70",
                   color = "black")+
  facet_wrap(~year)+
  labs(y = "Latitude",
       x = "Longitude") +
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y)+
  ggtitle("Env+Space+Time (AUC=0.95, Spearman's rho=0.46, PDE=0.35)")+
  theme_bw() + 
  theme(axis.title = element_text(size = 10),
        title = element_text(size = 10),
        axis.text = element_text(size = 6),
        legend.title =element_text(size = 10)) -> plot.3

# visualize residuals
ggplot(rr) + 
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   alpha = 0.4,
                   color = "black",
                   linewidth = 1)+
  geom_point(aes(y = y, x = x, color = residuals),size = 0.05) +
  scale_color_gradient2(midpoint = 0, high = "#6B6100", low = "#4F53B7", mid = "#F1F1F1") +
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   color = "black",
                   linewidth = 1)+
  ggplot2::geom_sf(data = region_layers$akland,
                   fill = "grey70",
                   color = "black")+
  #facet_wrap(~year)+
  labs(y = "Latitude",
       x = "Longitude") +
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y)+
  ggtitle("Env+Space+Time (AUC=0.95, Spearman's rho=0.46, PDE=0.35)")+
  theme_bw() + 
  theme(axis.title = element_text(size = 10),
        title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.title =element_text(size = 10)) -> plot.4



ggsave(plot = plot.1, "./Figures/residual_plotENVxYEAR.v1.png", width = 8.5, height = 8.5, units = "in") 
ggsave(plot = plot.3, "./Figures/residual_plotENV+SPACETIMExYEAR.v1.png", width = 8.5, height = 8.5, units = "in") 
ggsave(plot = plot.2, "./Figures/residual_plotENV.png", width = 8.5, height = 8.5, units = "in") 
ggsave(plot = plot.4, "./Figures/residual_plotENV+SPACETIME.png", width = 8.5, height = 8.5, units = "in") 

# non spacetime
rr <- na.omit(out.env)

rr %>%
  mutate(x_bin = cut_number(x, n = 50, labels = FALSE), y_bin = cut_number(y, n = 50, labels = FALSE)) %>%
  group_by(x_bin, y_bin, year) %>%
  reframe(mean.resid = mean(residuals))  -> pp

ggplot(pp) + 
  geom_tile(aes(y = y_bin, x = x_bin, fill = mean.resid))+
  facet_wrap(~year) +
  scale_fill_gradient2(midpoint = 0, high = "#6B6100", low = "#4F53B7", mid = "#F1F1F1", name = "Residuals")+
  ggtitle("Env (AUC=0.95, Spearman's rho=0.46, PDE=0.36)")+
  theme_bw() + 
  xlab("Longitude bin")+
  ylab("Latitude bin")+
  theme(axis.title = element_text(size = 10),
        title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.position = "bottom") -> plot.1


# spacetime
rr <- na.omit(out.st[,-c(1:6)])
colnames(rr) <- c("x", "y", "predicted", "residuals", "year")


rr %>%
  mutate(x_bin = cut_number(x, n = 50, labels = FALSE), y_bin = cut_number(y, n = 50, labels = FALSE)) %>%
  group_by(x_bin, y_bin, year) %>%
  reframe(mean.resid = mean(residuals))  -> pp

ggplot(pp) + 
  geom_tile(aes(y = y_bin, x = x_bin, fill = mean.resid))+
  facet_wrap(~year) +
  scale_fill_gradient2(midpoint = 0, high = "#6B6100", low = "#4F53B7", mid = "#F1F1F1", name = "Residuals")+
  ggtitle("Env+Space+Time (AUC=0.95, Spearman's rho=0.46, PDE=0.35)")+
  theme_bw() + 
  xlab("Longitude bin")+
  ylab("Latitude bin")+
  theme(axis.title = element_text(size = 10),
        title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.title =element_text(size = 10),
        legend.position = "bottom") -> plot.3


ggsave(plot = plot.1, "./Figures/residual_plotENVxYEAR.v2.png", width = 8.5, height = 11, units = "in") 
ggsave(plot = plot.3, "./Figures/residual_plotENV+SPACETIMExYEAR.v2.png", width = 8.5, height = 11, units = "in") 

ggarrange(plot.1, plot.3, common.legend = TRUE, nrow= 2, legend = "bottom")
ggsave("./Figures/resid.combined.png", width = 8.5, height = 11)
