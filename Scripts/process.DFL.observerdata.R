#PURPOSE: To process rkc data for use in BBRKC legal male fall distribution models

### LOAD LIBS/PARAMS ----
source("./Scripts/load.libs.params.R")


### PROCESS DF OBSERVER DATA -----------------------------------------------------------------
# Load 1995-1999 and 2000-2023 potsum data to get zero catch pots
# msr_pot(Y/N) denotes a measure pot vs. count pot, but both have total legal crab enumerated. Measure pot has 
# Detailed biological data collected (i.e. size, shell condition, etc.)
rbind(read.csv("./Data/RKC-1990-2021_potsum.csv"),
      read.csv("./Data/RKC-2022_23_potsum.csv")) %>%
  mutate(source = ifelse(description %in% c("BERING SEA BAIRDI CRAB", "BERING SEA BAIRDI", "EASTERN BERING SEA TANNER CRAB",
                                            "WESTERN BERING SEA TANNER CRAB", "PRIBILOF RED KING CRAB", 
                                            "CDQ PRIBILOF RED AND BLUE KING CRAB", "BERING SEA TANNER CRAB (B.S. WEST OF 166W"),
                         "Observer - bycatch", "Observer - directed fishery"),
         fishery = ifelse(source == "Observer - bycatch", "Tanner crab", "Red king crab")) %>%
  mutate(sampdate = lubridate::mdy(sampdate), 
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate),
         season = case_when(month %in% c(12, 1, 2) ~ "W",
                            month %in% c(3, 4, 5) ~ "Sp",
                            month %in% c(6, 7, 8) ~ "Su",
                            month %in% c(9, 10, 11) ~ "F")) %>%
  filter(latitude != -9 , longitude != -9) %>%
  group_by(year, month, day, season, latitude, longitude, source, fishery, adfg) %>%
  dplyr::reframe(catch_pp = sum(tot_legal),
                 cpue = catch_pp/sum(soaktime)) -> potsum

### PROCESS DF LOGBOOK DATA ---------------------------------------------------------------------------------
#Processed following Zacher et al. 2018 (need to remove strings >20km??)
dfl_0516 <- read.csv("./Data/dfl_0516.adfg.csv") # with adfg vessel code added
dfl_17<- read.csv("./Data/2017.18_BBRKC_DFLS.csv") # new data from Leah
dfl_18 <- read.csv("./Data/dfl_18.adfg.csv")
dfl_1920 <- read.csv("./Data/2019.20_BBRKC_DFLs.csv")
dfl_1920adfg <- read.csv("./Data/dfl_1920.adfg.csv") %>%
  dplyr::select(c(cif_form, adfg, packetnum, interview_date, haul_date, set_time, begin_lat, begin_lon, pots, crab_count))
dfl_1920 <- right_join(dfl_1920adfg, dfl_1920) %>%
  filter(is.na(adfg) == FALSE)
dfl_2021 <- read.csv("./Data/2020.21_BBRKC_DFLs.csv")  # new data from Leah
dfl_2021adfg <- read.csv("./Data/dfl_2021.adfg.csv") %>%
  dplyr::select(!c(X, crab_count, soaktime, statarea, frac, cpue, soakhrs))
dfl_2021 <- right_join(dfl_2021adfg, dfl_2021) #adding adfg code
dfl_2324 <- read.csv("./Data/dfl_2324.adfg.csv")
dfl_2425 <- read.csv("./Data/2024.25_BBRKC_DFLs.csv")

# Extract time informaiton, calculate catch per pot (ADD IN SOAK TIME??)
# 2005-2016
dfl_0516 %>%
  dplyr::filter(Pots_Hauled > 5 & Pots_Hauled <= 100) %>% #filtering for pots >5, and <=100
  mutate(lat = (Start_Latitude + End_Latitude)/2, #converting to mid-lat/long
         lon = (Start_Longitude + End_Longitude)/2,
         sampdate = lubridate::ymd(DATE_HAUL), 
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate),
         catch = RKC_Haul.Number,
         catch_pp = catch/Pots_Hauled) %>%
  rename(pots_total = Pots_Hauled) %>%
  dplyr::select(sampdate, year, month, day, lat, lon, catch, catch_pp, adfg, pots_total) -> dfl_0516CLEAN2

# 2017
dfl_17 %>%
  mutate(num = row_number()) %>%
 group_by(adfg, set_date, set_time, haul_date, haul_time, Soak, pots_total,
            crab_count, cpue, BegLat_DD, BeginLong_DD, EndLat_DD, EndLong_DD) %>%
  mutate(N = n()) %>%
  ungroup() -> dups

dups %>% 
  filter(N>1) %>%
  arrange(EndLat_DD, cpue) %>%
  mutate(num = row_number()) %>%
  filter(num %% 2 == 0) %>%
  rbind(., dups %>% filter(N==1))%>%
  dplyr::select(!c(num, N)) -> dfl_17CLEAN

write.csv(dfl_17CLEAN, "./Data/Clean DFL data/dfl_2017-18.csv")

dfl_17CLEAN %>%
  dplyr::filter(pots_total > 5 & pots_total <= 100) %>% #filtering for pots >5, and <=100
  mutate(lat = (BegLat_DD + EndLat_DD)/2, #converting to mid-lat/long
         lon = (BeginLong_DD + EndLong_DD)/2,
         sampdate = lubridate::mdy(haul_date), 
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate),
         catch = crab_count,
         catch_pp = catch/pots_total) %>%
  dplyr::select(sampdate, year, month, day, lat, lon, catch, catch_pp, adfg, pots_total) -> dfl_1718CLEAN2

# 2018-2019
dfl_18 %>%
  mutate(num = row_number()) %>%
  group_by(adfg, set_date, set_time, haul_date, haul_time, Soak, Pots_Total,
           crab_count, CPUE, begin_lat, begin_lon, end_lat, end_lon) %>%
  mutate(N = n()) %>%
  ungroup() -> dups

dups %>% 
  filter(N>1) %>%
  arrange(end_lat, CPUE) %>%
  mutate(num = row_number()) %>%
  filter(num %% 2 == 0) %>%
  rbind(., dups %>% filter(N==1))%>%
  dplyr::select(!c(num, N)) -> dfl_18CLEAN

write.csv(dfl_18CLEAN, "./Data/Clean DFL data/dfl_2018-19.csv")


dfl_18CLEAN %>%
  dplyr::filter(Pots_Total > 5 & Pots_Total <= 100) %>% #filtering for pots >5, and <=100
  mutate(lat = (begin_lat + end_lat)/2, #converting to mid-lat/long
         lon = (begin_lon + end_lon)/2,
         sampdate = lubridate::ymd(haul_date), 
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate),
         catch = crab_count,
         catch_pp = catch/Pots_Total) %>%
  rename(pots_total = Pots_Total) %>%
  filter(is.na(sampdate) == FALSE) %>%
  dplyr::select(sampdate, year, month, day, lat, lon, catch, catch_pp, adfg, pots_total) -> dfl_1819CLEAN2

# 2019-2020
dfl_1920 %>%
  mutate(num = row_number()) %>%
  group_by(adfg, packetnum, set_date, set_time, haul_date, haul_time, pots, crab_count, lost_pots,
          begin_lat, begin_lon, end_lat, end_lon) %>%
  mutate(N = n()) %>%
  ungroup() -> dups

dups %>% 
  filter(N>1) %>%
  arrange(end_lat, crab_count) %>%
  mutate(num = row_number()) %>%
  filter(num %% 2 == 0) %>%
  rbind(., dups %>% filter(N==1))%>%
  dplyr::select(!c(num, N)) -> dfl_1920CLEAN

write.csv(dfl_1920CLEAN, "./Data/Clean DFL data/dfl_2019-20.csv")

 
dfl_1920CLEAN %>%
  mutate(lat = (begin_lat + end_lat)/2, #converting to mid-lat/long
         lon = (begin_lon + end_lon)/2,
         sampdate = lubridate::mdy(haul_date),
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate),
         catch = crab_count,
         catch_pp = catch/(pots-lost_pots),
         pots_total = pots-lost_pots) %>%
  dplyr::select(sampdate, year, month, day, lat, lon, catch, catch_pp, adfg, pots_total) -> dfl_1920CLEAN2

# 2020-2021 
dfl_2021 %>%
  mutate(num = row_number()) %>%
  group_by(packetnum, adfg, interview_date, trip_week_start, trip_week_end, daysfished, set_date, set_time, pots_total,
          haul_date, crab_count, cpue, begin_lat, end_lat, begin_lon, end_lon) %>%
  mutate(N = n()) %>%
  ungroup() -> dups

dups %>% 
  filter(N>1) %>%
  arrange(end_lat, cpue) %>%
  mutate(num = row_number()) %>%
  filter(num %% 2 == 0) %>%
  rbind(., dups %>% filter(N==1)) %>%
  dplyr::select(!c(num, N)) -> dfl_2021CLEAN

write.csv(dfl_2021CLEAN, "./Data/Clean DFL data/dfl_2020-21.csv")

  
dfl_2021CLEAN %>%
  dplyr::filter(pots_total > 5 & pots_total <= 100) %>% #filtering for pots >5, and <=100
  mutate(lat = (begin_lat + end_lat)/2, #converting to mid-lat/long
         lon = (begin_lon + end_lon)/2,
         sampdate = lubridate::mdy(haul_date),
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate),
         catch = crab_count,
         catch_pp = catch/pots_total) %>%
  dplyr::select(sampdate, year, month, day, lat, lon, catch, catch_pp, adfg, pots_total) -> dfl_2021CLEAN2

# 2023-2024
dfl_2324 %>%
  mutate(num = row_number()) %>%
  group_by(packetnum, adfg, Interview_Date, Trip_Wk_Start, Trip_Wk_End, Days_Fished, Set_Date, Set_Time, Pot_Count,
           Lost_Pots, Haul_Date, Crab_Count, Begin_Lat., End_Lat., Begin_Lon., End_Lon.) %>%
  mutate(N = n()) %>%
  ungroup() -> dups

dups %>% 
  filter(N>1) %>%
  arrange(End_Lon., Crab_Count)  %>%
  mutate(num = row_number()) %>%
  filter(num %% 2 == 0) %>%
  rbind(., dups %>% filter(N==1)) %>%
  dplyr::select(!c(num, N)) -> dfl_2324CLEAN

write.csv(dfl_2324CLEAN, "./Data/Clean DFL data/dfl_2023-24.csv")


dfl_2324CLEAN %>%
  mutate(pots_total = Pot_Count-Lost_Pots) %>%
  rename(pots = Pot_Count, crab_count = Crab_Count) %>%
  dplyr::filter(pots > 5 & pots <= 100) %>% #filtering for pots >5, and <=100
  mutate(lat = (Begin_Lat. + End_Lat.)/2, #converting to mid-lat/long
         lon = (Begin_Lon. + End_Lon.)/2,
         sampdate = lubridate::mdy(Haul_Date),
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate),
         catch = crab_count,
         catch_pp = catch/pots_total) %>%
  dplyr::select(sampdate, year, month, day, lat, lon, catch, catch_pp, adfg, pots_total) -> dfl_2324CLEAN2

# 2024-2025
dfl_2425 %>%
  mutate(num = row_number()) %>%
  group_by(ADFG, Packet, Interview_Date, Trip_Wk_Start, Trip_Wk_End, Days_Fished, Set_Date, Set_Time, Pots_Hauled,
           Lost_Pots, Haul_Date, Haul_Time, Crab_Count, Begin_Lat., End_Lat., Begin_Lon., End_Lon., CPUE) %>%
  mutate(N = n()) %>%
  ungroup() -> dups

dups %>% 
  filter(N>1) %>%
  arrange(End_Lon., Crab_Count)  %>%
  mutate(num = row_number()) %>%
  filter(num %% 2 == 0) %>%
  rbind(., dups %>% filter(N==1)) %>%
  dplyr::select(!c(num, N)) -> dfl_2425CLEAN

write.csv(dfl_2425CLEAN, "./Data/Clean DFL data/dfl_2024-25.csv")


dfl_2425CLEAN %>%
  mutate(pots_total = Pot_Count-Lost_Pots) %>%
  rename(pots = Pot_Count, crab_count = Crab_Count) %>%
  dplyr::filter(pots > 5 & pots <= 100) %>% #filtering for pots >5, and <=100
  mutate(lat = (Begin_Lat. + End_Lat.)/2, #converting to mid-lat/long
         lon = (Begin_Lon. + End_Lon.)/2,
         sampdate = lubridate::mdy(Haul_Date),
         year = lubridate::year(sampdate),
         month = lubridate::month(sampdate),
         day = lubridate::day(sampdate),
         catch = crab_count,
         catch_pp = catch/pots_total) %>%
  dplyr::select(sampdate, year, month, day, lat, lon, catch, catch_pp, pots_total) -> dfl_2425CLEAN2

# Join data frames, specify seasons
#rbind(dfl_0516, dfl_17, dfl_18, dfl_1920, dfl_20) %>%
rbind(dfl_0516CLEAN2, dfl_1718CLEAN2, dfl_1819CLEAN2, dfl_1920CLEAN2, dfl_2021CLEAN2, dfl_2324CLEAN2) %>%
  mutate(season = case_when(month %in% c(12, 1, 2) ~ "W",
                            month %in% c(3, 4, 5) ~ "Sp",
                            month %in% c(6, 7, 8) ~ "Su",
                            month %in% c(9, 10, 11) ~ "F")) %>%
  rename(latitude = lat, longitude = lon) %>%
  na.omit(catch_pp) -> dfls.1

# Get unique dfl adfg codes by year, month, day fishing
dfls.1 %>%
  dplyr::select(year, month, day, adfg) %>%
  distinct() -> dfls.adfg

# Specify potsum rownumbers for filtering later
potsum %>%
  mutate(rownum = 1:nrow(potsum)) -> potsum2

# Join potsum and dfls to match duplicate entries by year, month, day, and adfg code
right_join(potsum2, dfls.adfg, by = c("year", "month", "day", "adfg")) %>%
  filter(is.na(latitude) == FALSE) -> dfl.dups

# Filter potsum observer data to exclude rownumbers for the duplicate entries with DFLs
potsum2 %>%
  filter(!rownum %in% dfl.dups$rownum)%>%
  group_by(year, season, latitude, longitude, source, fishery)  %>%
  dplyr::reframe(catch_pp = sum(catch_pp),
                 cpue = sum(cpue)) -> potsum.clean
 
# Summarise dfls by year, season, coordinatees
dfls.1 %>%
  group_by(year, season, latitude, longitude) %>%
  reframe(catch_pp = sum(catch_pp)) -> dfls.clean

### JOIN DF LOGBOOK AND OBSERVER DATA ----------------------------------------------------------------------
rbind(potsum.clean %>%
        dplyr::select(!cpue),
      dfls.clean %>%
        mutate(source = "Logbook",
               fishery = "Red king crab")) %>%
  filter(-150>longitude & -170<longitude, latitude<60 & latitude>54) -> lm_df

lm_df %>%
  sf::st_as_sf(coords = c(x = "longitude", y = "latitude"), crs = sf::st_crs(in.crs)) %>%
  sf::st_transform(crs = map.crs) %>%
  terra::vect() %>%
  terra::mask(., BB_strata) -> data_vect

cbind(crds(data_vect), as.data.frame(data_vect)) %>%
  sf::st_as_sf(coords = c(x = "x", y = "y"), crs = sf::st_crs(map.crs)) %>%
  sf::st_transform(crs = in.crs) -> data_vect2

cbind(st_coordinates(data_vect2), as.data.frame(data_vect2)) %>%
  dplyr::select(year, season, Y, X, catch_pp, source, fishery) %>%
  rename(longitude = X, latitude = Y) %>%
  filter(is.na(catch_pp) == FALSE) -> lm_df2


write.csv(lm_df2, "./Data/legalmale_direct.fish.csv") 
