
source("./Scripts/load.libs.params.R")


Fall_lm.preds <- rast("./Data/Fall_lm.preds.tif")

btemp <- subset(Fall_lm.preds, grep("mean.temp ", names(Fall_lm.preds)))
head(as.data.frame(btemp[[1]]))


### Compute anomalies -----
JF <- subset(btemp, grep("Jan/Feb ", names(btemp))) %>%
      cbind(as.data.frame(.), crds(.)) %>%
        dplyr::select(28:56) %>%
      mutate(cell.mean = rowMeans(.[,-c(28:29)])) %>%
      pivot_longer(., cols = 1:27, names_to = "time", values_to = "sst") %>%
      mutate(period = substr(time, 1, 7),
             year = substr(time, 19, 22),
             anom = sst - cell.mean) %>%
      dplyr::select(!time) 
  
MA <- subset(btemp, grep("Mar/Apr ", names(btemp))) %>%
      cbind(as.data.frame(.), crds(.)) %>%
      dplyr::select(28:56) %>%
      mutate(cell.mean = rowMeans(.[,-c(28:29)])) %>%
      pivot_longer(., cols = 1:27, names_to = "time", values_to = "sst") %>%
      mutate(period = substr(time, 1, 7),
             year = substr(time, 19, 22),
             anom = sst - cell.mean) %>%
      dplyr::select(!time) 


MJ <- subset(btemp, grep("May/Jun ", names(btemp))) %>%
      cbind(as.data.frame(.), crds(.)) %>%
      dplyr::select(28:56) %>%
      mutate(cell.mean = rowMeans(.[,-c(28:29)])) %>%
      pivot_longer(., cols = 1:27, names_to = "time", values_to = "sst") %>%
      mutate(period = substr(time, 1, 7),
             year = substr(time, 19, 22),
             anom = sst - cell.mean) %>%
      dplyr::select(!time) 

JA <- subset(btemp, grep("Jul/Aug ", names(btemp))) %>%
      cbind(as.data.frame(.), crds(.)) %>%
      dplyr::select(28:56) %>%
      mutate(cell.mean = rowMeans(.[,-c(28:29)])) %>%
      pivot_longer(., cols = 1:27, names_to = "time", values_to = "sst") %>%
      mutate(period = substr(time, 1, 7),
             year = substr(time, 19, 22),
             anom = sst - cell.mean) %>%
      dplyr::select(!time) 

SO <- subset(btemp, grep("Sep/Oct ", names(btemp))) %>%
      cbind(as.data.frame(.), crds(.)) %>%
      dplyr::select(27:54) %>%
      mutate(cell.mean = rowMeans(.[,-c(27:28)])) %>%
      pivot_longer(., cols = 1:26, names_to = "time", values_to = "sst") %>%
      mutate(period = substr(time, 1, 7),
             year = substr(time, 19, 22),
             anom = sst - cell.mean) %>%
      dplyr::select(!time) 

ND <- subset(btemp, grep("Nov/Dec ", names(btemp))) %>%
      cbind(as.data.frame(.), crds(.)) %>%
      dplyr::select(27:54) %>%
      mutate(cell.mean = rowMeans(.[,-c(27:28)])) %>%
      pivot_longer(., cols = 1:26, names_to = "time", values_to = "sst") %>%
      mutate(period = substr(time, 1, 7),
             year = substr(time, 19, 22),
             anom = sst - cell.mean) %>%
      dplyr::select(!time) 


# Bind all dfs
df <- rbind(JF, MA, MJ, JA, SO, ND)

df <- rbind(JF, JA, SO)


df %>%
  group_by(year) %>%
  reframe(mean.anom = mean(anom)) %>%
  mutate(cc = ifelse(mean.anom < 0, "Cold", "Warm")) -> anom.df

# Plot
ggplot(anom.df, aes(year, mean.anom, color = cc)) +
  geom_line(group = 1, color = "black")+
  geom_point(size = 2)+
  scale_color_manual(values = c("cyan3", "brown1"), labels = c("cold", "warm"))+
  theme_bw()

### Read in spatial model predictions for directed fishery -----
spat.df <- read.csv("./Output/spatial_predictions.csv")

westBLZ <- westBLZ - RKCSA_sub - area512

# Join closure areas
bind_rows(list(st_as_sf(RKCSA_sub) %>% mutate(AREA = "Red King Crab Savings Area"), 
               st_as_sf(area512) %>% mutate(AREA = "NMFS Statistical Area 512"),
               st_as_sf(westBLZ) %>% mutate(AREA = "Bycatch Limitation Zone 1 (west of 162°)"),
               st_as_sf(NBBTCA) %>% mutate(AREA = "Nearshore Bristol Bay trawl closure area"))) -> closure_areas

# Join with spatial predictions df
spat.df %>%
  st_as_sf(., coords = c("x", "y"), crs = map.crs) %>% # transform into spatial object
  na.omit(.) %>% 
  st_join(., closure_areas) %>%
  cbind(., st_coordinates(.)) %>% # transform from spatial object to data frame
  as.data.frame(.) %>%
  dplyr::select(!c(X, FID, geometry)) %>%
  mutate(AREA = case_when((is.na(AREA) == TRUE) ~ "Other Bristol Bay",
                          TRUE ~ AREA)) %>%
  filter(AREA %in% c("Bycatch Limitation Zone 1 (west of 162°)", "Red King Crab Savings Area", "NMFS Statistical Area 512"))-> area.df

# calculate proportion in each area
area.df %>%
  group_by(year) %>%
  mutate(total_catch = sum(catch_pp)) %>%
  ungroup() %>%
  group_by(year, AREA) %>%
  reframe(total_catch = unique(total_catch),
          area_catch = sum(catch_pp),
          prop_catch = area_catch/total_catch) -> prop.df

# join with sst anom df
right_join(prop.df, anom.df %>% mutate(year = as.numeric(year)), by = "year") %>%
  na.omit() %>%
  mutate(cc = case_when((year == 2017) ~ "Cold",
                        TRUE ~ cc)) -> prop.sst.df

prop.sst.df %>%
  filter(year == 2019) %>%
  mutate(year = 2020, total_catch = NA, area_catch = NA, prop_catch = NA, cc = "Warm") -> dummy

rbind(prop.sst.df, dummy) -> prop.sst.df

# Plot
# ggplot()+
#   geom_bar_pattern(prop.sst.df, mapping = aes(as.factor(year), prop_catch, color = cc, pattern = AREA,
#                                               pattern_fill = cc,
#                                               pattern_color = cc),
#                    stat = "identity", position = "stack", pattern_density = 0.2, pattern_spacing = 0.025, fill = "white",
#                    pattern_key_scale_factor = 0.6)+
#   scale_pattern_manual(values = c(`Area 512` = "stripe",
#                                   `Red king crab savings area` = "none"), name = "")+
#   theme_bw()+
#   scale_color_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
#   scale_pattern_fill_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
#   scale_pattern_color_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
#   scale_x_discrete(breaks = seq(min(prop.sst.df$year), max(prop.sst.df$year), by = 5))+
#   xlab("Year")+
#   ylab("Proportion of predicted catch")+
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 12),
#         legend.key = element_rect(colour = "black", fill = "black"))

ggplot()+
  geom_point(prop.sst.df, mapping = aes(as.factor(year), prop_catch, color = cc, shape = AREA), size = 2)+
  scale_color_manual(values = c(cold = "#A1A6C8", warm = "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
  theme_bw()

ggplot()+
  geom_bar(prop.sst.df, mapping = aes(as.factor(year), prop_catch, fill = cc), color = "lightgrey",
                   stat = "identity")+
  theme_bw()+
  facet_wrap(~AREA)+
  #scale_color_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
  scale_fill_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
  scale_x_discrete(breaks = seq(min(prop.sst.df$year), max(prop.sst.df$year), by = 5))+
  xlab("Year")+
  ggtitle("Predicted data")
  ylab("Proportion of predicted catch")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) -> pred.plot
  
  
### Read in observations for directed fishery ----
# Load data
read.csv("./Output/lm_F_train.csv") %>%
  dplyr::select(!X) -> train
read.csv("./Output/lm_F_test.csv") %>%
  dplyr::select(!X) -> test

rbind(train, test) %>%
  filter(iter == lm_iter) -> dat

### Join observations with spatial predictions df ----
dat %>%
  st_as_sf(., coords = c("x", "y"), crs = map.crs) %>% # transform into spatial object
  na.omit(.) %>% 
  st_join(., closure_areas) %>%
  cbind(., st_coordinates(.)) %>% # transform from spatial object to data frame
  as.data.frame(.) %>%
  dplyr::select(!c(FID, geometry)) %>%
  mutate(AREA = case_when((is.na(AREA) == TRUE) ~ "Other Bristol Bay",
                          TRUE ~ AREA)) %>%
  filter(AREA %in% c("Bycatch Limitation Zone 1 (west of 162°)", "Red King Crab Savings Area", "NMFS Statistical Area 512"))-> area.df2

# calculate proportion in each area
area.df2 %>%
  group_by(predict_year) %>%
  mutate(total_catch = sum(catch_pp)) %>%
  ungroup() %>%
  group_by(predict_year, AREA) %>%
  reframe(total_catch = unique(total_catch),
          area_catch = sum(catch_pp),
          prop_catch = area_catch/total_catch) %>%
  dplyr::rename(year = predict_year) -> prop.df2

# join with sst anom df
right_join(prop.df2, anom.df %>% mutate(year = as.numeric(year)), by = "year") %>%
  na.omit() %>%
  mutate(cc = case_when((year == 2017) ~ "Cold",
                        TRUE ~ cc)) -> prop.sst.df2

prop.sst.df2 %>%
  filter(year == 2019) %>%
  mutate(year = 2020, total_catch = NA, area_catch = NA, prop_catch = NA, cc = "Warm") -> dummy

rbind(prop.sst.df2, dummy) -> prop.sst.df2

rbind(prop.sst.df %>% mutate(type = "Predicted"), prop.sst.df2 %>% mutate(type = "Observed")) -> dat2

# Plot
ggplot()+
  geom_bar(dat2 %>% filter(type == "Observed"), mapping = aes(as.factor(year), prop_catch, fill = cc), color = "lightgrey",
           stat = "identity")+
  theme_bw()+
  facet_wrap(~factor(AREA, levels = c("Bycatch Limitation Zone 1 (west of 162°)", "Red King Crab Savings Area", 
                                            "NMFS Statistical Area 512", "Nearshore Bristol Bay trawl closure area")), 
              labeller = label_wrap_gen(width = 28, multi_line = TRUE), nrow = 3)+
  #scale_color_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
  scale_fill_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
  scale_x_discrete(breaks = seq(min(prop.sst.df$year), max(prop.sst.df$year), by = 5))+
  xlab("Year")+
  ylab("Predicted abundance proportion")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none", 
        legend.direction = "horizontal") -> plot.out

ggsave(plot = plot.out, "./Figures/warmVcold_catchprop.png", width = 6, height = 8, units = "in")


## Top 5 warmest and coldest year (via anomalies) 
dat2 %>%
  filter(cc == "Warm", type == "Predicted", year!=2020) %>%
  dplyr::select(year, mean.anom) %>%
  distinct() %>%
  top_n(5, mean.anom) %>%
  pull(year) -> warm.yrs

dat2 %>%
  filter(cc == "Cold", type == "Predicted", year!=2020) %>%
  dplyr::select(year, mean.anom) %>%
  distinct() %>%
  top_n(-5, mean.anom) %>%
  pull(year) -> cold.yrs

top.yrs <- c(warm.yrs, cold.yrs)

dat2 %>%
  filter(year %in% top.yrs, type == "Predicted") -> dat3

#Plot
ggplot()+
  geom_bar(dat3, mapping = aes(factor(year, levels = c("1999", "2008", "2009", "2010", "2012",
                                                       "2003", "2015", "2016", "2018", "2019")), 
                               prop_catch, fill = cc, group = cc), color = "lightgrey",
           stat = "identity", position = "dodge")+
  theme_bw()+
  facet_wrap(~factor(AREA, levels = c("Bycatch Limitation Zone 1 (west of 162°)", "Red King Crab Savings Area", 
                                      "NMFS Statistical Area 512", "Nearshore Bristol Bay trawl closure area")), 
         labeller = label_wrap_gen(width = 28, multi_line = TRUE))+
  #scale_color_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
  scale_fill_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
  xlab("Year")+
  ylab("Proportion of catch")+
  theme(axis.text.x = element_text(size = 10.5, angle = 45, vjust = 1, hjust=1),
        axis.text.y = element_text(size = 10.5),
        axis.title = element_text(size = 10.5),
        strip.text = element_text(size = 10.5),
        legend.text = element_text(size = 10.5),
        legend.position = "bottom", 
        legend.direction = "horizontal") -> plot.out

ggsave(plot = plot.out, "./Figures/predicted_top5warmVcold_catchprop.png", width = 8.5, height = 6, units = "in")

dat3 %>%
  group_by(AREA, cc) %>%
  reframe(prop_catch2 = mean(prop_catch),
          N = n(),
          std = sd(prop_catch),
          Neff = coda::effectiveSize(prop_catch),
          se = std/sqrt(Neff)) -> hh

ggplot()+
  geom_bar(hh, mapping = aes(cc, prop_catch2, fill = cc), color = "lightgrey",
           stat = "identity", position = "dodge")+
  theme_bw()+
  geom_errorbar(hh, mapping = aes(x = cc, ymin = prop_catch2 - se, ymax = prop_catch2 + se), width = 0)+
  facet_wrap(~factor(AREA, levels = c("Bycatch Limitation Zone 1 (west of 162°)", "Red King Crab Savings Area", 
                                      "NMFS Statistical Area 512", "Nearshore Bristol Bay trawl closure area")), 
             labeller = label_wrap_gen(width = 28, multi_line = TRUE))+
  #scale_color_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
  scale_fill_manual(values = c("#A1A6C8", "#CA9CA4"), labels = c("Cold", "Warm"), name = "")+
  xlab("Year")+
  ylab("Predicted abundance proportion")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none", 
        legend.direction = "horizontal") -> plot.out3

ggsave(plot = plot.out3, "./Figures/predicted_top5warmVcold_AVGcatchprop.png", width = 8.5, height = 6, units = "in")



### Correlations between observed and predicted ----
dat2 %>%
  filter(type == "Predicted") %>%
  dplyr::select(year, total_catch) %>%
  distinct() -> pp

dat2 %>%
  filter(type == "Observed") %>%
  dplyr::select(year, total_catch) %>%
  distinct() -> oo

# Predictions on test data
# Read in training/testing data
lm_F_train <- read.csv("./Output/lm_F_train.csv") %>%
  filter(iter == lm_iter)
test <- read.csv("./Output/lm_F_test.csv") %>%
  filter(iter == lm_iter)

# Best models 
lm.F.modelb <- readRDS("./Models/lm.modelb.F.8.rda")
lm.F.modelp <- readRDS("./Models/lm.modelp.F.8.rda")

# Predict
pred.b <- predict.gbm(lm.F.modelb, # fitted model to predict
                      test, # data to predict to
                      n.trees=lm.F.modelb$gbm.call$best.trees, # see help
                      type="response") # predict probabilities

pred.p <- predict.gbm(lm.F.modelp, # fitted model to predict
                      test %>% filter(catch_pp>0), # data to predict to
                      n.trees=lm.F.modelp$gbm.call$best.trees, # see help
                      type="response") # predict probabilities

# control for overdispersion
pred.p2 <- pred.p * theta 

pred.dat <- data.frame(test %>% filter(catch_pp>0), pred.p = pred.p2)
# Bind predictions
pred.b2 <- cbind(test, pred.b)

pred.dat <- right_join(pred.dat, pred.b2) %>%
  mutate(pred = pred.p * pred.b)

#Pyper and Peterman code written by Franz Mueter, last updated 2000
# Required function for cor.test.PP:
N.effective <- function(mat) {
  # written by Franz Mueter   Last modified: 23 October 2000
  # function to compute effective sample size for pairwise correlations among 
  # autocorrelated variables. 
  # Based on Pyper & Peterman (1998). CJFAS 55:2127-2140.  Eq. (1)
  # where summation was done over j = 1, ..., N/5
  # 
  # mat a matrix of variables, one variable per column
  #
  # function to compute simple estimates of autocorrelation up to lag N/5: 
  # (Eq. 7 in Pyper & Peterman)
  ar.fun <- function(x, max.lag = ceiling(sum(!is.na(x))/5)) {
    res <- rep(NA, max.lag)
    n <- length(x)
    for(i in 1.:max.lag) {
      x.bar <- mean(x, na.rm = T)
      res[i] <- ((n/(n - i)) * sum((x[1:(n - i)] - x.bar) * (x[(i + 1):n] - x.bar), na.rm
                                   = T))/sum((x - x.bar)^2, na.rm = T)
    }
    res
  }
  AR <- apply(mat, 2., ar.fun)
  k <- ncol(mat)
  if(is.matrix(AR)) {
    AR1 <- vector("list", k)
    for(i in 1:k) AR1[[i]] <- AR[, i]
    AR <- AR1  
  }
  N <- t(!is.na(mat)) %*% (!is.na(mat))
  N.lags <- ceiling(N/5.)
  N.eff <- matrix(0., k, k)
  # constrain effective N to smaller than or equal to actual N:
  for(i in 1.:k) {
    for(j in 1.:i) {
      lags <- 1.:N.lags[i, j]
      Nij <- N[i, j]
      N.eff[i, j] <- round((1./Nij + 2./Nij * sum((Nij - lags)/Nij * AR[[i]][lags] * AR[[
        j]][lags]))^(-1.))
    }
  }
  j <- N.eff > N
  N.eff[j] <- N[j]
  N.eff + t(N.eff) - diag(diag(N.eff))
}

# Function to test for significant correlation between two time series
# in the presence of autocorrelation in one or both of the series:
cor.test.PP <- function(x, y) {
  # Function to test for significant correlations between x and y, which may be autocorrelated, 
  # using modified Chelton method after Pyper & Peterman (1998)
  # Eqn. 3 with N*-2 degrees of freedom
  N.eff <- N.effective(cbind(x, y))[1, 2]
  r <- cor(x, y, use="pair")
  fun <- function(alpha, N, r) {
    t2 <- qt(1 - alpha/2, N - 2)^2
    sqrt(t2/(t2 + N - 2)) - abs(r)
  }
  p.value <- uniroot(fun, c(1e-015, 0.9999), N = N.eff, r = r)$root
  cat("Two-sided test\n\n")
  c(correlation = r, P.value = p.value, N.eff= N.eff)
}

cor.test.PP(scale(pred.dat$pred), scale(pred.dat$catch_pp)) -> cor.out

plot(scale(pp$total_catch), type = "l", col = "blue")
lines(scale(oo$total_catch), col = "red")

cor.out
