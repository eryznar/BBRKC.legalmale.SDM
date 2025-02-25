
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
               st_as_sf(NBBTCA) %>% mutate(AREA = "Nearshore Bristol Bay trawl closure area"))) %>%
  mutate(KM2 = as.numeric(st_area(.)/1e6)) -> closure_areas

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
  filter(AREA %in% c("Bycatch Limitation Zone 1 (west of 162°)", "Red King Crab Savings Area", "NMFS Statistical Area 512", "Other Bristol Bay"))-> area.df

# calculate proportion and densityin each area
area.df %>%
  group_by(year) %>%
  mutate(total_catch = sum(catch_pp),
         total_area = sum(KM2)) %>%
  ungroup() %>%
  group_by(year, AREA) %>%
  reframe(total_catch = unique(total_catch),
          total_area = unique(total_area),
          area_catch = sum(catch_pp),
          prop_catch = area_catch/total_catch,
          area_density = area_catch/KM2) %>%
  distinct()-> prop.df

# join with sst anom df
right_join(prop.df, anom.df %>% mutate(year = as.numeric(year)), by = "year") %>%
  #na.omit() %>%
  mutate(cc = case_when((year == 2017) ~ "Cold",
                        TRUE ~ cc)) %>%
  filter(is.na(AREA) == FALSE) -> prop.sst.df

prop.sst.df %>%
  filter(year == 2019) %>%
  mutate(year = 2020, total_catch = NA, total_area = NA,
         area_catch = NA, prop_catch = NA, area_density = NA, cc = "Warm") -> dummy

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
read.csv("./Output/New models/lm_F_train.csv") %>%
  dplyr::select(!X) -> train
read.csv("./Output/New models/lm_F_test.csv") %>%
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
  filter(AREA %in% c("Bycatch Limitation Zone 1 (west of 162°)", "Red King Crab Savings Area", "NMFS Statistical Area 512", "Other Bristol Bay"))-> area.df2

# calculate proportion in each area
area.df2 %>%
  group_by(predict_year) %>%
  mutate(total_catch = sum(catch_pp),
         total_area = sum(KM2)) %>%
  ungroup() %>%
  group_by(predict_year, AREA) %>%
  reframe(total_catch = unique(total_catch),
          total_area = unique(total_area),
          area_catch = sum(catch_pp),
          prop_catch = area_catch/total_catch,
          area_density = area_catch/KM2) %>%
  dplyr::rename(year = predict_year) %>%
  distinct() %>%
  filter(is.na(AREA) == FALSE) -> prop.df2



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
  filter(year %in% top.yrs, type == "Predicted", AREA != "Other Bristol Bay") -> dat3

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

dat2 %>%
  filter(year %in% top.yrs, type == "Predicted") %>%
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
  ylab("Predicted proportion of abundance")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "none", 
        legend.direction = "horizontal") -> plot.out3

ggsave(plot = plot.out3, "./Figures/predicted_top5warmVcold_AVGcatchprop.png", width = 8.5, height = 6, units = "in")

# DOES AVERAGING SPATIAL DISTRIBUTIONS ACROSS ANOMALY YEARS PRODUCE THE SAME RESULTS AS IN FIG 6? ----
### Read in spatial model predictions for directed fishery
spat.df <- read.csv("./Output/spatial_predictions.csv")

# Join with spatial predictions df
right_join(spat.df, anom.df %>% mutate(year = as.numeric(year)), by = "year") %>%
  na.omit() %>%
  mutate(cc = case_when((year == 2017) ~ "Cold",
                        TRUE ~ cc)) %>%
  filter(year %in% top.yrs) -> sst.df

# calculate mean catch_pp
sst.df %>%
  group_by(cc, x, y) %>%
  reframe(mean_catchpp = mean(catch_pp)) -> mean.df

# Calculate and plot quantile distributions
clim <- unique(mean.df$cc)

poly.list <- list()

for(ii in 1:length(clim)){
  mean.df %>%
    filter(cc == clim[ii]) -> plot.df
  
  # Find plotting breaks
  quantiles = c(.05, .25, .5, .75)
  quants<-sort(unique(c(0,quantiles,1)))
  
  threshold <- 0.0513
  
  sample <- stats::na.omit(plot.df$mean_catchpp)
  sample[sample <= threshold] <- NA
  perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
  perc.breaks[1]<-0
  perc.breaks[length(perc.breaks)]<-Inf
  
  
  # Make raster again
  plot.df %>%
    dplyr::select(!cc) %>%
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
  
  polys <- list(percpoly2, percdummy3)
  
  poly.list <- c(poly.list, polys)
  
}


# Set up plot boundary
plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                    x = c(-167.5, -158)) # plot boundary unprojected

plot.boundary <- plot.boundary.untrans %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(crs = map.crs) %>%
  sf::st_coordinates() %>%
  as.data.frame() %>%
  dplyr::rename(x = X, y = Y) # plot boundary projected

year_untrans <- data.frame(lab = c("Warm"),
                           x = -158.3, y = 55.3)

year_lab <- year_untrans %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
  sf::st_transform(crs = map.crs) %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame()

# Map
ggplot2::ggplot() +
  ggplot2::geom_sf(data = poly.list[[1]], ggplot2::aes(fill = as.factor(layer)), color = NA) +
  ggplot2::geom_sf(data = st_as_sf(poly.list[[2]]),fill=NA, size = .3) +
  ggplot2::geom_sf(data = st_as_sf(area512),
                   fill = NA,
                   color = "purple",
                   linewidth = 0.75)+
  ggplot2::geom_sf(data = st_as_sf(westBLZ),
                   fill = NA,
                   color = "blue",
                   linewidth = 1)+
  ggplot2::geom_sf(data = st_as_sf(RKCSA),
                   fill = NA,
                   color = "red",
                   linewidth = 0.75)+
  ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                   fill = NA,
                   color = "red",
                   linewidth = 0.75)+
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   color = "black",
                   linewidth = 1)+
  ggplot2::geom_sf(data = region_layers$akland, 
                   fill = "grey70", 
                   color = "black")+
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y)+
  geom_text(data = data.frame(lab = c("Cold"),
                              x = -158.3, y = 55.3) %>%
              sf::st_as_sf(., coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
              sf::st_transform(., crs = map.crs) %>%
              cbind(st_coordinates(.)) %>%
              as.data.frame(), 
              aes(x=X, y=Y, label= lab), fontface = "bold", size=4) +
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
    plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> perc_rast_cold

# Map
ggplot2::ggplot() +
  ggplot2::geom_sf(data = poly.list[[3]], ggplot2::aes(fill = as.factor(layer)), color = NA) +
  ggplot2::geom_sf(data = st_as_sf(poly.list[[4]]),fill=NA, size = .3) +
  ggplot2::geom_sf(data = st_as_sf(area512),
                   fill = NA,
                   color = "purple",
                   linewidth = 0.75)+
  ggplot2::geom_sf(data = st_as_sf(westBLZ),
                   fill = NA,
                   color = "blue",
                   linewidth = 1)+
  ggplot2::geom_sf(data = st_as_sf(RKCSA),
                   fill = NA,
                   color = "red",
                   linewidth = 0.75)+
  ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                   fill = NA,
                   color = "red",
                   linewidth = 0.75)+
  ggplot2::geom_sf(data = st_as_sf(BB_strata),
                   fill = NA,
                   color = "black",
                   linewidth = 1)+
  ggplot2::geom_sf(data = region_layers$akland, 
                   fill = "grey70", 
                   color = "black")+
  coord_sf(xlim = plot.boundary$x,
           ylim = plot.boundary$y)+
  geom_text(data = data.frame(lab = c("Warm"),
                              x = -158.3, y = 55.3) %>%
              sf::st_as_sf(., coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
              sf::st_transform(., crs = map.crs) %>%
              cbind(st_coordinates(.)) %>%
              as.data.frame(), 
            aes(x=X, y=Y, label= lab), fontface = "bold", size=4) +
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
    legend.position = "bottom", plot.title = element_text(size = 11),
    plot.background = ggplot2::element_rect(fill = "white", color = "white")) -> perc_rast_warm


# Arrange plots
ggarrange(perc_rast_warm, perc_rast_cold, heights = c(1,1), nrow= 1, ncol = 2,
          common.legend = TRUE, legend = "bottom") -> warm.v.coldrast

ggsave(plot = warm.v.coldrast, "./Figures/anomaly_warmVcold_preds.png", height=3, width=6, units="in")

