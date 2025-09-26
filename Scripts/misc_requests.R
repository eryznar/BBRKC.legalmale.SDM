# PURPOSE:
# To respond to miscellaneous requests from reviewers

# AUTHOR: 
# Emily Ryznar - NOAA/NMFS RACE Shellfish Assessment Program (emily.ryznar@noaa.gov)

# LOAD LIBS/PARAMS ------
source("./Scripts/load.libs.params.R")

# Calculating sampling distribution by year and data type
# Read in response data
lm_df <- read.csv("./Data/legalmale_direct.fish.csv")

lm_df %>%
  group_by(year, source) %>%
  reframe(total_obs = n()) -> jj

jj %>%
  pivot_wider(., names_from = source, values_from = total_obs) -> pp

pp %>%
  group_by(year) %>%
  mutate(total = sum(na.omit(c(`Observer - directed fishery`, `Observer - bycatch`, Logbook)))) -> hh

View(hh)

# Evaluating poisson dist
lm_df %>%
  filter(year > 1996) -> dd

hist(dd$catch_pp, breaks = 50, main = "", xlab = "Catch per pot")
ggsave(plot,file="graph1.pdf")

png("./Figures/respdata_hist.png", height = 4, width = 5, units = "in", res = 600)

install.packages("MASS")
library(MASS)

rbind(lm_F_train, lm_F_test) -> dat
dat$catch_pp <- as.integer(dat$catch_pp)

colnames(dat)

nb_model <- glm.nb(catch_pp ~ Sep_Oct_BT + Nov_Dec_BT + Sep_Oct_currentE + Nov_Dec_currentE +Sep_Oct_currentN+
                   Nov_Dec_currentN + Jan_Feb_BT + Jul_Aug_BT + Jan_Feb_Ice + Mar_Apr_Ice + Sed+ SAP_Count + Depth+            
                   Slope+ Jan_Feb_currentE + Mar_Apr_currentE + May_Jun_currentE + Jan_Feb_currentN + Mar_Apr_currentN+
                   May_Jun_currentN+ Jul_Aug_currentN + Tidemax, data = dat)
summary(nb_model)
nb_model$theta
