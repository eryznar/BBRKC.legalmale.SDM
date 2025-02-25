#PURPOSE:
# To fit and evaluate boosted regression trees to model Bristol Bay red king crab
# legal male distribution in the fall

#AUTHOR:
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)

### LOAD PROCESSING PARAMETERS --------------------------------------------------
source("./Scripts/load.libs.params.R")

### GENERATE TRAINING/TESTING DATA ----------------------------------------------------------------------------------------
# Read in spatial covariates
Fall_lm.preds <- rast("./Data/Fall_lm.preds.tif")

# Read in response data
lm_df <- read.csv("./Data/legalmale_direct.fish.csv") %>%
  dplyr::select(!X)

# Write processing function to extract covariates at catch locations, create training/testing data

# Filter response data by season, transform into spatial vectors, mask by BB management area
lm_df %>%
  filter(season == "F") %>%
  sf::st_as_sf(coords = c(x = "longitude", y = "latitude"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(crs = map.crs) %>%
  terra::vect() %>%
  terra::mask(BB_strata) -> data_vect

preds <- subset(Fall_lm.preds, grep("GFbycatch.F", names(Fall_lm.preds), invert = TRUE))

# Generate for loop to extract covariate values at bycatch presence points by year, store as list of dfs,
# assign labels
pst <- c(1997:2019, 2021:2022)
predict_yr <- c(1998:2019, 2021:2023)

datalist = list()
datalist = vector("list", length = length(predict_yr))

for (ii in 1:length(predict_yr)){
  
  preds_past <- subset(preds, grep("Sep/Oct |Nov/Dec ", names(preds)))
  preds_pres <- subset(preds, grep("Sep/Oct |Nov/Dec ", names(preds), invert = TRUE))
  
  
  dat_pres <- terra::extract(subset(preds_pres, grep(predict_yr[ii], names(preds_pres), value = T)), 
                             subset(data_vect, data_vect$year == predict_yr[ii]))
  
  dat_past <- terra::extract(subset(preds_past, grep(pst[ii], names(preds_past), value = T)), 
                             subset(data_vect, data_vect$year == predict_yr[ii]))
  
  dat <- cbind(dat_past[,-1], dat_pres[,-1],
               crds(subset(data_vect, data_vect$year == predict_yr[ii])),
               subset(data_vect[,c("year", "catch_pp", "fishery")], data_vect$year == predict_yr[ii]))
  
  labs <- c("Sep_Oct_BT", "Nov_Dec_BT", "Sep_Oct_currentE", "Nov_Dec_currentE", "Sep_Oct_currentN",
            "Nov_Dec_currentN", 
            "Jan_Feb_BT", "Mar_Apr_BT", "May_Jun_BT", "Jul_Aug_BT", 
            "Jan_Feb_Ice", "Mar_Apr_Ice", 
            "Sed", "SAP_Count", "Depth", "Slope", 
            "Jan_Feb_currentE", "Mar_Apr_currentE", "May_Jun_currentE", "Jul_Aug_currentE",
            "Jan_Feb_currentN", "Mar_Apr_currentN", "May_Jun_currentN", "Jul_Aug_currentN",
            "Tidemax") #raster labs
  
  
  names(dat) <- c(labs, "x", "y", "predict_year", "catch_pp", "fishery")
  
  datalist[[ii]] <- dat
  
}

# Bind data by year
dplyr::bind_rows(datalist) %>%
  na.omit() %>%
  mutate(#geographic_position = Latitude * Longitude,
    #Current = CurrentN*CurrentE,
    #CurrentSD = CurrentSDE*CurrentSDN,
    PA = ifelse(catch_pp == 0, 0, 1)) -> sdmData # bind presence dfs


# Check for variable collinearity, drop collinear vars
non_covs <- c("x", "y", "predict_year", "catch_pp", "PA", "fishery")

sdmData[-which(colnames(sdmData) %in% non_covs)] %>%
  as.data.frame() %>%
  na.omit() %>%
  vifstep(th=5) -> vif_drop 

dropvars <- vif_drop@excluded 

#dropvars <- dropvars[-which(dropvars %in% c("Jul_Aug_SST_min", "Sep_Oct_SST_min", "Jul_Aug_SST_max", "bottom_temp_max", "Sed"))]

sdmData <- sdmData[-which(colnames(sdmData) %in% dropvars)]

# Create training/testing datasets, randomly drawing 10 times
iter <- seq(1, 10, by =1)

trainlist = list()
trainlist = vector("list", length = length(iter))

testlist = list()
testlist = vector("list", length = length(iter))

set.seed(1)

for (ii in 1:length(iter)){
  group1 <- kfold(sdmData, k = 5)
  trainlist[[ii]] <- sdmData[group1 != 1,]
  testlist[[ii]]  <- sdmData[group1 == 1, ]
  
  trainlist[[ii]]$iter <- iter[ii]
  testlist[[ii]]$iter <- iter[ii]
}

train <- dplyr::bind_rows(trainlist)
test <- dplyr::bind_rows(testlist)


write.csv(train, "./Output/New models/lm_F_train.csv")
write.csv(test, "./Output/New models/lm_F_test.csv")


### CREATE FUNCTION TO FIT MODELS ON REPLICATE TRAINING/TESTING SPLITS --------------
model_iter <- function(train, test, seas, iteration){
  # Filter by iteration
  train2 <- train %>% 
    filter(iter %in% iteration) %>%
    mutate(catch_pp = as.integer(catch_pp),
           fishery = as.factor(fishery))
  
  test2 <- test %>% 
    filter(iter %in% iteration) %>%
    mutate(catch_pp = as.integer(catch_pp),
           fishery = as.factor(fishery))
  
  # Specify non-predictor variables
  non_covs <- c("x", "y", "predict_year", "catch_pp", "PA", "iter", "fishery")
  
  # Fit binomial model
  model_b <- gbm.step(data = train2, 
                      gbm.x = which(!colnames(train2) %in% non_covs), # columns of predictors
                      gbm.y = which(colnames(train2) == "PA"), # column of response
                      family = "bernoulli", # for PA data
                      tree.complexity = 5, # model interactions, etc?
                      learning.rate = 0.05, # influence of each tree 
                      bag.fraction = 0.5) 
  
  saveRDS(model_b, paste0("./Models/New models/lm.modelb.", seas,".", iteration, ".rda"))
  
  ntreesb <- model_b$n.trees
  
  model_p <- gbm.step(data = train2 %>% filter(catch_pp >0), 
                      gbm.x = which(!colnames(train2) %in% non_covs), # columns of predictors 
                      gbm.y = which(colnames(train2) == "catch_pp"), # column of response
                      family = "poisson", # for count data
                      tree.complexity = 5, # model interactions, etc?
                      learning.rate = 0.05, # influence of each tree 
                      bag.fraction = 0.5)
  
  saveRDS(model_p, paste0("./Models/New models/lm.modelp.", seas,".", iteration, ".rda"))
  
  ntreesp <- model_p$n.trees
  
  # Calculate AUC for binomial
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
  
  
  # Calculate RMSE and spearman's for positive model
  pred <- predict.gbm(model_p, # fitted model to predict
                      test2 %>% filter(catch_pp>0), # data to predict to
                      n.trees=model_p$gbm.call$best.trees, # see help
                      type="response") # predict probabilities
  
  obs <- test2$catch_pp[which(test2$catch_pp >0)]
  
  RMSE <- sqrt(mean((obs-pred)^2))
  
  d <- as.data.frame(cbind(obs, pred))
  
  ggplot(d, aes(x = obs, y = pred))+
    geom_point()+
    geom_smooth(method = "lm") + 
    theme_classic()+
    xlab("Observed vals")+
    ylab("Predicted vals") -> pred_obs_plot
  
  cor.test(d$obs, d$pred, method = "spearman") -> cor_out
  
  rho <- cor_out$estimate
  p <- cor_out$p.value
  
  # Bind metrics into a data frame
  eval_df <- data.frame(ntreesb = ntreesb, ntreesp = ntreesp, AUC = AUC, RMSE = RMSE, rho = rho, p = p, 
                        iter = iteration)
  
  write.csv(eval_df, "./Models/New models/model_iteration.csv")
  
  return(list(eval_df))
}

### RUN FUNCTION TO FIT MODELS ------------------------------------------------------
read.csv("./Output/New models/lm_F_train.csv") %>%
  dplyr::select(!X) -> train
read.csv("./Output/New models/lm_F_test.csv") %>%
  dplyr::select(!X) -> test

1:10 %>%
  map_df(~model_iter(train, test, "F", .x)) -> model_out

saveRDS(model_out, "./Models/New models/eval_out.rda")

readRDS("./Models/New models/eval_out.rda") # iteration 6!
