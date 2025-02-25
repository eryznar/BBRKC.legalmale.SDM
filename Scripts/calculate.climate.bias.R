### LOAD PROCESSING PARAMETERS 
source("./Scripts/load.libs.params.R")


# Read in spatial covariates
preds <- rast("./Data/Fall_lm.preds.tif")

# Read in response data
lm_df <- read.csv("./Data/legalmale_direct.fish.csv")

# Read in training/testing data
train <- read.csv("./Output/New models/lm_F_train.csv") %>%
  filter(iter == lm_iter)
test <- read.csv("./Output/New models/lm_F_test.csv") %>%
  filter(iter == lm_iter)


# Function to calculate Cohen's D and Hellinger's distance between Bristol Bay and area sampled by fishery
#gsub('[0-9]+', '', names(preds)) -> names(preds)

# Temp
names(train)[2:23] <- c("Sep/Oct mean.temp ", "Nov/Dec mean.temp ", "Sep/Oct mean.uEast ", "Nov/Dec mean.uEast ",
                  "Sep/Oct mean.vNorth ", "Nov/Dec mean.vNorth ", "Jan/Feb mean.temp ", "Jul/Aug mean.temp ",
                  "Mar_Apr Ice ", "Jan_Feb Ice ", "Sed ", "SAP Count ", "Depth ", "Slope ", 
                  "Jan/Feb mean.uEast ", "Mar/Apr mean.uEast ", "May/Jun mean.uEast ", "Jan/Feb mean.vNorth ",
                  "Mar/Apr mean.vNorth ", "May/Jun mean.vNorth ", "Jul/Aug mean.vNorth ", "Tidemax ")

#subset(preds, names(preds) %in% names(train)[2:23]) -> preds2

nn <- names(train)[2:23]
#nn<- nn[c(1:2, 7:14)]

norm_vec <- function(x) sqrt(sum(x^2))

hell_dist <- function (p, q, from, to, n = 1024) {
  P <- density(p, kernel = "gaussian", from = from, to = to, n = n)
  p <- P$y
  p <- p / sum(p)
  Q <- density(q, kernel = "gaussian", from = from, to = to, n = n)
  q <- Q$y
  q <- q / sum(q)
  hd <- norm_vec(sqrt(p) - sqrt(q)) / sqrt(2)
  hd
}

eval.out <- data.frame()

for(ii in 1:length(nn)){
  
  # subset preds by variable
  pp <- subset(preds, grep(nn[ii], names(preds)))
  dd <- as.data.frame(pp) %>%
    mutate(index = 1:nrow(.))
  
  pivot_longer(dd, !index, names_to = nn[ii], values_to = "value") -> kk
  
  kk %>% pull(value) -> all.BB
  
  # subset variable in training data
  pref.samp <- train %>%
                  dplyr::select(nn[ii]) %>%
                  pull()
  
  n_bins <- 1024 #arbitrary
  
  bins <- seq(floor(min(c(all.BB, pref.samp))), ceiling(max(c(all.BB, pref.samp))), length = n_bins)
  
  # Calculate cohen's D
  cd <- (cohen.d(all.BB, pref.samp))$estimate
  
  # Calculate hellinger's distance
  hd<- hell_dist(all.BB, pref.samp, from=bins[1], to=bins[n_bins])
  #hd <- hellinger(all.BB, pref.samp, lower = -Inf, upper = Inf)
  
  eval.out <- rbind(eval.out, data.frame(var = nn[ii],
                                         cd = cd,
                                         hd = hd))
  
}
