
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


### Bycatch vs. directed fishery overlap-----
pred.df <- rbind(read.csv("./Output/spatial_predictions.csv") %>% 
                   mutate(model = "Directed fishery") %>%
                   dplyr::select(!X) %>%
                   rename(value = catch_pp),
                 read.csv("./Data/predicted_LMbycatch.csv") %>% 
                   mutate(model = "Bycatch") %>%
                   dplyr::select(!X) %>%
                   rename(value = mean_count))

# Spatial maps by period
pp <- c("1998:2005", "2006:2014", "2015:2023")


pred.df %>%
  filter(year > 1997) %>%
  mutate(phase = case_when((year %in% 1998:2005) ~ "1998:2005",
                           (year %in% 2006:2014) ~ "2006:2014",
                           (year %in% 2015:2023) ~ "2015:2023")) -> pred.df2

# Set up list to store maps
spatpred_list <- list()



# Calculate encounter percentiles
for(ii in 1:length(pp)){
 # BYCATCH
   bycatch <-  pred.df2 %>%
    filter(phase %in% pp[ii], model == "Bycatch", year > 1997) %>%
    group_by(x, y) %>%
    reframe(value = mean(value))
  
  # Find plotting breaks
  quantiles = c(.05, .25, .5, .75)
  quants<-sort(unique(c(0,quantiles,1)))
  
  threshold <- 0.0513
  
  sample <- stats::na.omit(bycatch$value)
  sample[sample <= threshold] <- NA
  perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
  perc.breaks[1]<-0
  perc.breaks[length(perc.breaks)]<-Inf
  
  # Make raster again
  bycatch %>%
    rast() %>%
    raster() -> bycatch.rast
  
  # Set crs
  crs(bycatch.rast) <- map.crs
  
  # Cut the prediction map by the breaks
  perc.map <- raster::cut(bycatch.rast, breaks = perc.breaks)
  
  # set up the factor maps
  perc.vals <- raster::getValues(perc.map)
  perc.vals[perc.vals == 1] <- NA
  
  # convert the raster to polygons
  percpoly0 <- stars::st_as_stars(perc.map)
  percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
  bycatch.poly <- percpoly[percpoly$layer %in% 4:5, ]
  
# DIRECTED FISHERY
  direct.fish <-  pred.df2 %>%
    filter(phase %in% pp[ii], model == "Directed fishery", year > 1997) %>%
    group_by(x, y) %>%
    reframe(value = mean(value))
  
  # Mask to areas where bycatch can occur
  direct.fish %>%
    st_as_sf(., coords = c("x", "y"), crs = map.crs) -> qq
  
  mask(vect(qq), NBBTCA, inverse = TRUE) -> kk
  mask(kk, RKCSA, inverse = TRUE) -> tt
  
  # change back to df
  cbind(crds(tt), as.data.frame(tt)) -> direct.fish
  
  # Find plotting breaks
  quantiles = c(.05, .25, .5, .75)
  quants<-sort(unique(c(0,quantiles,1)))
  
  threshold <- 0.0513
  
  sample <- stats::na.omit(direct.fish$value)
  sample[sample <= threshold] <- NA
  perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
  perc.breaks[1]<-0
  perc.breaks[length(perc.breaks)]<-Inf
  
  # Make raster again
  direct.fish%>%
    rast() %>%
    raster() -> direct.fishrast
  
  # Set crs
  crs(direct.fishrast) <- map.crs
  
  # Cut the prediction map by the breaks
  perc.map <- raster::cut(direct.fishrast, breaks = perc.breaks)
  
  # set up the factor maps
  perc.vals <- raster::getValues(perc.map)
  perc.vals[perc.vals == 1] <- NA
  
  # convert the raster to polygons
  percpoly0 <- stars::st_as_stars(perc.map)
  percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
  directfish.poly <- percpoly[percpoly$layer %in% 4:5, ]
  
  # Set up plot boundary
  plot.boundary.untrans <- data.frame(y = c(54.25, 59.25),
                                      x = c(-167.5, -158)) # plot boundary unprojected
  
  plot.boundary <- plot.boundary.untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    sf::st_coordinates() %>%
    as.data.frame() %>%
    dplyr::rename(x = X, y = Y) # plot boundary projected
  
  
  # Set up year label size
  size = 3.5
  lw = 1
  
  if(pp[ii] == "1998:2005"){
    labs = "1998-\n2005"
  }else if(pp[ii] == "2006:2014"){
    labs = "2006-\n2014"
  }else{
    labs = "2015-\n2023"
  }
  
  
  year_untrans <- data.frame(lab = labs, x = -158.3, y = 55.3)
  
  year_lab <- year_untrans %>%
    sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(4326)) %>%
    sf::st_transform(crs = map.crs) %>%
    cbind(st_coordinates(.)) %>%
    as.data.frame()
  
  
  # Map
  ggplot2::ggplot() +
    #ggplot2::geom_sf(data = survey.sf, fill = "grey95")+
    # ggplot2::geom_sf(data = bycatch.poly, ggplot2::aes(fill = as.factor(layer)), col = NA, alpha = 0.5) +
    # ggplot2::geom_sf(data = directfish.poly, ggplot2::aes(fill = as.factor(layer)), col = NA, alpha = 0.5) +
    ggplot2::geom_sf(data = bycatch.poly, mapping = aes(fill = as.factor(1)), col = NA, alpha = 0.5) +
    ggplot2::geom_sf(data = directfish.poly, mapping = aes(fill = as.factor(2)), col = NA, alpha = 0.5) +
    scale_fill_manual(values = c("darkgoldenrod1", "darkcyan"), name = "", labels = c("Bycatch", "Directed fishery"))+
    #ggplot2::geom_sf(data = st_as_sf(percdummy3),fill=NA, size = .3) +
    # ggplot2::geom_sf(data = st_as_sf(NBBTCA),
    #                  fill = NA,
    #                  color = "black",
    #                  linewidth = lw)+
    ggplot2::geom_sf(data = st_as_sf(westBLZ),
                     fill = NA,
                     color = "blue",
                     linewidth = lw)+
    ggplot2::geom_sf(data = st_as_sf(area512),
                     fill = NA,
                     color = "purple",
                     linewidth = lw)+
    ggplot2::geom_sf(data = st_as_sf(RKCSA_sub),
                     fill = NA,
                     color = "red",
                     linewidth = lw)+
    ggplot2::geom_sf(data = st_as_sf(RKCSA),
                     fill = NA,
                     color = "red",
                     linewidth = lw)+
    ggplot2::geom_sf(data = st_as_sf(BB_strata),
                     fill = NA,
                     color = "black",
                     linewidth = 1.75)+
    ggplot2::geom_sf(data = region_layers$akland,
                     fill = "grey70",
                     color = "black")+
    geom_text(data = year_lab, aes(x=X, y=Y, label= lab), fontface = "bold", size=size) +
    #ggtitle(title)+
    
    coord_sf(xlim = plot.boundary$x,
             ylim = plot.boundary$y)+
    #viridis::scale_fill_viridis(option = "cividis", discrete = T, name = "Percentiles", labels = c("95%", "75%", "50%", "25%")) +
    ggplot2::theme_bw() +
    ggplot2:: theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      #legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      #legend.position = legend.pos,
      panel.grid.major = element_blank(),
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 11), legend.title = ggplot2::element_text(size = 11),
      legend.position = "bottom", plot.title = element_text(size = 18)) -> spatpred_list[[ii]]
  
}

# Arrange plots
ggarrange(spatpred_list[[1]],
          spatpred_list[[2]],
          spatpred_list[[3]],
          nrow=3, ncol=1, common.legend = TRUE, legend = "bottom") -> phase.perc_rast

ggsave(plot = phase.perc_rast, "./Figures/predicted.bycatchVdirectedfish.png", width = 6, height = 8, units = "in")

# Percent overlap by year
df <- pred.df %>% filter(model == "Directed fishery", year == 1998)

# Set up list to store maps
overlap.df <- data.frame()

yy <- c(1998:2019, 2021:2023)

# Calculate encounter percentiles
for(ii in 1:length(yy)){
  # BYCATCH
  bycatch <-  pred.df2 %>%
    filter(year %in% yy[ii], model == "Bycatch") %>%
    group_by(x, y) %>%
    reframe(value = mean(value))
  
  # Find plotting breaks
  quantiles = c(.05, .25, .5, .75)
  quants<-sort(unique(c(0,quantiles,1)))
  
  threshold <- 0.0513
  
  sample <- stats::na.omit(bycatch$value)
  sample[sample <= threshold] <- NA
  perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
  perc.breaks[1]<-0
  perc.breaks[length(perc.breaks)]<-Inf
  
  # Make raster again
  bycatch %>%
    rast() %>%
    raster() -> bycatch.rast
  
  # Set crs
  crs(bycatch.rast) <- map.crs
  
  # Cut the prediction map by the breaks
  perc.map <- raster::cut(bycatch.rast, breaks = perc.breaks)
  
  # set up the factor maps
  perc.vals <- raster::getValues(perc.map)
  perc.vals[perc.vals == 1] <- NA
  
  # convert the raster to polygons
  percpoly0 <- stars::st_as_stars(perc.map)
  percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
  bycatch.poly <- percpoly[percpoly$layer %in% 4:5, ]
  
  st_union(bycatch.poly) -> byc.poly
  
  byc.area <- st_area(byc.poly)
  
  # DIRECTED FISHERY
  direct.fish <-  pred.df2 %>%
    filter(year %in% yy[ii], model == "Directed fishery") %>%
    group_by(x, y) %>%
    reframe(value = mean(value))
  
  # Mask to areas where bycatch can occur
  direct.fish %>%
    st_as_sf(., coords = c("x", "y"), crs = map.crs) -> qq
  
  mask(vect(qq), NBBTCA, inverse = TRUE) -> kk
  mask(kk, RKCSA, inverse = TRUE) -> tt
  
  # change back to df
  cbind(crds(tt), as.data.frame(tt)) -> direct.fish
  
  # Find plotting breaks
  quantiles = c(.05, .25, .5, .75)
  quants<-sort(unique(c(0,quantiles,1)))
  
  threshold <- 0.0513
  
  sample <- stats::na.omit(direct.fish$value)
  sample[sample <= threshold] <- NA
  perc.breaks <- stats::quantile(sample, probs = quants, na.rm = TRUE, names = FALSE)
  perc.breaks[1]<-0
  perc.breaks[length(perc.breaks)]<-Inf
  
  # Make raster again
  direct.fish%>%
    rast() %>%
    raster() -> direct.fishrast
  
  # Set crs
  crs(direct.fishrast) <- map.crs
  
  # Cut the prediction map by the breaks
  perc.map <- raster::cut(direct.fishrast, breaks = perc.breaks)
  
  # set up the factor maps
  perc.vals <- raster::getValues(perc.map)
  perc.vals[perc.vals == 1] <- NA
  
  # convert the raster to polygons
  percpoly0 <- stars::st_as_stars(perc.map)
  percpoly <- sf::st_as_sf(percpoly0,merge = TRUE)
  directfish.poly <- percpoly[percpoly$layer %in% 4:5, ]
 
  st_union(directfish.poly) -> df.poly
  
  df.area <- st_area(df.poly)
  
  # Calculate area of overlap
  st_intersection(df.poly, byc.poly) -> intersect
  
  intersect.area <- st_area(intersect) 
  
  # Join data
  df <- data.frame(year = yy[ii], 
                           bycatch.area = byc.area, 
                           directfish.area = df.area, 
                           overlap.area = intersect.area) %>%
                mutate(bycatch.area = as.numeric(bycatch.area/1000), # convert to km2
                       directfish.area = as.numeric(directfish.area/1000),
                       overlap.area = as.numeric(overlap.area/1000),
                       fishery.area = sum(bycatch.area, directfish.area) - overlap.area,
                       total.overlap = (overlap.area/fishery.area)*100,
                       bycatch.overlap = (overlap.area/bycatch.area)*100,
                       df.overlap = (overlap.area/directfish.area)*100)
  
  overlap.df <- rbind(overlap.df, df)
  
}

dummy <- overlap.df %>% 
          filter(year == 1998) %>%
          mutate(year = 2020,
                 bycatch.area = NA,
                 directfish.area = NA,
                 overlap.area = NA,
                 fishery.area = NA,
                 total.overlap = NA,
                 bycatch.overlap = NA,
                 df.overlap = NA)

overlap.df <- rbind(dummy, overlap.df)

# Plot
ggplot()+
  geom_line(overlap.df, mapping = aes(year, total.overlap), color = "cadetblue", linewidth = 1)+
  geom_point(overlap.df, mapping = aes(year, total.overlap), color = "cadetblue", size = 2)+
  theme_bw()+
  ylab("Percent overlap")+
  xlab("Year")

ggplot() +
  geom_line(overlap.df, mapping = aes(year, scale(bycatch.area), color = "1"), linewidth = 1)+
  geom_point(overlap.df, mapping = aes(year, scale(bycatch.area), color = "1"), size = 2)+
  geom_line(overlap.df, mapping = aes(year, scale(directfish.area), color = "2"), linewidth = 1)+
  geom_point(overlap.df, mapping = aes(year, scale(directfish.area), color = "2"), size = 2)+
  theme_bw()+
  scale_color_manual(values = c("darkgoldenrod", "steelblue"), labels = c("Bycatch", "Directed fishery"))+
  ylab("Predicted 50th percentile area (km^2)")+
  xlab("Year")
  
  cor.test(scale(overlap.df$bycatch.area), scale(overlap.df$directfish.area))
  
  
  # Spatial correlation plots
  head(pred.df)
  
  pred.df %>%
    filter(year > 1997) %>%
    mutate(phase = case_when((year %in% 1998:2005) ~ "1998:2005",
                             (year %in% 2006:2014) ~ "2006:2014",
                             (year %in% 2015:2023) ~ "2015:2023")) -> pred.df2
  
  pred.df2 %>% 
    group_by(x, y, year, model) %>%
    reframe(value = mean(value)) -> cor.dat
  
 
  corr.df <- data.frame()
  yrs <- c(1998:2019, 2021:2023)
  
  for(ii in 1:length(yrs)){
    cor.dat %>%
      filter(year == yrs[ii]) -> dat
    
    dat %>%
      filter(model == "Directed fishery") %>%
      dplyr::select(x, y, value) %>%
      rast(.) %>%
      raster(.) -> df
    
    dat %>% 
      filter(model == "Bycatch") %>%
      dplyr::select(x, y, value) %>%
      rast(.) %>%
      raster(.) -> byc
    
    
    resample(byc, df) -> byc
    
    corr.map <- SpatialPack::modified.ttest(scale(as.data.frame(byc)$value), 
                                             scale(as.data.frame(df)$value),
                                             coords = coordinates(byc))
    
    ss <- data.frame(year = yrs[ii], 
                     corr = round(corr.map$corr, 2),
                     p = round(corr.map$p.value, 2))
    
    corr.df <- rbind(corr.df, ss)
  }
  

  
  resample(early.bycatch, early.df) -> early.bycatch
  
  as.data.frame(early.bycatch)
  
  
  corr.map1 <- SpatialPack::modified.ttest(scale(as.data.frame(early.bycatch)$value), 
                                           scale(as.data.frame(early.df)$value),
                                           coords = coordinates(early.bycatch))
  

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
