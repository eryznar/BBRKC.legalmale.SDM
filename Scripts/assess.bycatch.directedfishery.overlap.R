source("./Scripts/load.libs.params.R")



### Bycatch vs. directed fishery overlap-----
pred.df <- rbind(read.csv(paste0(dir, "spatial_predictions.csv")) %>% 
                   mutate(model = "Directed fishery") %>%
                   dplyr::select(!X) %>%
                   rename(value = catch_pp),
                 read.csv(paste0(dir, "predicted_LMbycatch.csv")) %>% 
                   mutate(model = "Bycatch") %>%
                   dplyr::select(!X) %>%
                   rename(value = mean_count))

### Weighted coordinates ---
# mask df data to where trawl fishing can occur
pred.df %>%
  st_as_sf(., coords = c("x", "y"), crs = map.crs) -> qq

mask(vect(qq), NBBTCA, inverse = TRUE) -> kk
mask(kk, RKCSA, inverse = TRUE) -> tt

# change back to df
cbind(crds(tt), as.data.frame(tt)) -> pred.df2

# Set up list to store maps
wtdcrds.df <- data.frame()

yy <- c(1998:2019, 2021:2023)

# Calculate encounter percentiles
for(ii in 1:length(yy)){
  # BYCATCH
  bycatch <-  pred.df %>%
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
  
  #byc.vals <- mask(bycatch.rast, bycatch.poly)
  byc.vals <- bycatch.rast
  
  top.byc <- cbind(coordinates(byc.vals), as.data.frame(byc.vals))
  
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
  
  # Resample to bycatch raster resolution (coarser)
  direct.fishrast <- resample(direct.fishrast, bycatch.rast)
  
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
  
  #df.vals <- mask(direct.fishrast, directfish.poly)
  df.vals <- direct.fishrast
  
  top.df <- cbind(coordinates(df.vals), as.data.frame(df.vals))
  
  thres <- mean(na.omit(top.df$value))
  
  
  # Join data and calc wtd crds
  rbind(top.byc %>% mutate(model = "Bycatch"), top.df %>% mutate(model = "Direct.fish")) %>%
    filter(is.na(value) == FALSE) %>%
   mutate(year = yy[ii]) %>%
   group_by(year, model) %>%
   mutate(total = sum(value),
          norm_val = value/total) %>%
   ungroup() %>%
   group_by(year, model) %>%
   reframe(N = n(),
           mean_weighted_lat = sum(y*norm_val)/sum(norm_val),
           mean_weighted_lon = sum(x*norm_val)/sum(norm_val)) -> wtd.crds
  
  wtdcrds.df <- rbind(wtdcrds.df, wtd.crds)
  
}

wtdcrds.df %>%
  group_by(model) %>%
  mutate(scale_lat = scale(mean_weighted_lat)[,1],
         scale_lon = scale(mean_weighted_lon)[,1]) %>%
  ungroup() -> wtd.crds2

# Plot
ggplot(wtd.crds2, aes(year, scale_lat, color = model)) +
  geom_line(linewidth = 1)+
  # geom_line(wtd_crds2, mapping = aes(year, bycatch_wtd_lat, color = as.factor(1)), linewidth = 1)+
  # geom_line(wtd_crds2, mapping = aes(year, directfish_wtd_lat, color = as.factor(2)), linewidth = 1) + 
  theme_bw()+
  scale_x_continuous(breaks = seq(1998, 2023, by = 2))+
  scale_color_manual(values = c("darkturquoise", "darkgoldenrod"), labels = c("Bycatch",
                                                                              "Directed fishery"),
                     name = "")+
  ylab("Weighted latitude")

wtd.crds2 %>%
  dplyr::select(year, model, scale_lat) %>%
  pivot_wider(names_from= model, values_from = scale_lat) -> jj


cor.test.PP(jj$Bycatch, jj$Direct.fish) -> out

ggplot(wtd.crds2, aes(scale_lon, scale_lat, color = model))+
  #geom_point()+
  theme_bw()+
  geom_text(aes(label = year))+
  xlab("Bycatch weighted latitude")+
  ylab("Directed fishery weighted latitude")
  


# Spatial maps by period ---
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
    mutate(value = rescale(value)) %>% # rescaling for warren's I and Schoener's D
    rast(.) %>%
    raster(.) -> df
  
  dat %>% 
    filter(model == "Bycatch") %>%
    dplyr::select(x, y, value) %>%
    mutate(value = rescale(value)) %>%
    rast(.) %>%
    raster(.) -> byc
  
  
  resample(byc, df) -> byc
  
  over <- ENMTools::raster.overlap(scale(rast(df)), scale(rast(byc)))
  
  
  n_bins <- 1024 #arbitrary
  
  bins <- seq(floor(min(na.omit(c(values(byc), values(df))))), ceiling(max(c(na.omit(values(byc), values(df))))), length = n_bins)
  
  # Calculate hellinger's distance
  hd<- hell_dist(scale(na.omit(values(byc))), scale(na.omit(values(df))), from=bins[1], to=bins[n_bins])
  #hd <- hellinger(all.BB, pref.samp, lower = -Inf, upper = Inf)
  
  
  ss <- data.frame(year = yrs[ii], 
                   corr = round(over$rank.cor, 2),
                   I = round(over$I, 2),
                   SD = round(over$D, 2),
                   HD = hd)
  
  corr.df <- rbind(corr.df, ss)
}

corr.df %>%
  filter(year == 2021) %>%
  mutate(year = 2020,
         I = NA,
         corr = NA,
         HD = NA,
         SD = NA) -> dummy

rbind(dummy, corr.df) -> corr.df

ggplot()+
  geom_line(corr.df, mapping = aes(year, corr, color = "1"), linewidth = 1)+
  geom_point(corr.df, mapping = aes(year, corr, color = "1"), size = 1.5)+
  # geom_line(corr.df, mapping = aes(year, I, color = "2"), linewidth = 1)+
  # geom_point(corr.df, mapping = aes(year, I, color = "2"), size = 1.5)+
  geom_line(corr.df, mapping = aes(year, SD, color = "3"), linewidth = 1)+
  geom_point(corr.df, mapping = aes(year, SD, color = "3"), size = 1.5)+
  theme_bw()+
  scale_color_manual(values = c("darkgoldenrod", "#009B95"), 
                     labels = c("Spearman's rho", "Schoener's D"),
                     name = "")+
  ylab("Value")+
  xlab("")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12)) -> p1


# Plot
ggplot()+
  geom_line(overlap.df, mapping = aes(year, total.overlap), color = "black", linewidth = 1)+
  geom_point(overlap.df, mapping = aes(year, total.overlap), color = "black", size = 2)+
  theme_bw()+
  ylab("Percent overlap")+
  xlab("Year")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12)) -> p2

p1 + p2 + plot_layout(nrow = 2)

ggsave("./Figures/bycatchVdirectedfish_corrTS.png", width = 6, height = 6, units = "in")


