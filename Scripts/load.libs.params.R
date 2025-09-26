### LOAD PACKAGES -----------------------------------------------------------------------------------------
library(sf)
library(ggmap)
#library(rgdal)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(raster)
library(terra)
library(lubridate)
library(akgfmaps)
library(biomod2)
library(gam)
library(purrr)
library(effsize)
library(dismo)
library(gbm) 
library(pROC)
library(ggrepel)
library(geosphere)
library(usdm)
library(statip)
library(ggpattern)
library(corrplot)


### SET SPATIAL DETAILS ---------------------------------------------------------
crs.latlon <- "epsg:4326" #lat lon crs

map.crs <- "EPSG:3338"

in.crs = "+proj=longlat +datum=NAD83"

dir <- "Y:/KOD_Research/Ryznar/BBRKC legal male SDM/Data/"

# LOAD SPATIAL LAYERS -----------------------------------------------------------
region_layers <- akgfmaps::get_base_layers(select.region = "bs.south", set.crs="auto")

survey_gdb <- "./Data/SAP_layers.gdb"

# Bristol Bay strata multipolygon
BB_strata <- sf::st_read(survey_gdb, layer = "BristolBaySurveyStrata") %>%
              # can also use sf::st_read() to read layers in as sf objects
              vect() %>%
            project(., map.crs)


st_read(paste0(dir, "Closure areas/RKCSA_sub.shp")) %>%
  vect() -> RKCSA_sub

st_read(paste0(dir,"Closure areas/RKCSA.shp")) %>%
  vect() -> RKCSA

st_read(paste0(dir, "Closure areas/area512.shp")) %>%
  vect() -> area512
st_read(paste0(dir, "Closure areas/fullBLZ.shp")) %>%
  vect() -> fullBLZ

st_read(paste0(dir, "Closure areas/BLZwestof162.shp")) %>%
  vect() -> westBLZ

st_read(paste0(dir, "Closure areas/NBBTCA.shp")) %>%
  vect() -> NBBTCA


# OTHER PARAMETERS --------------------------------------------------------------
lm_iter <- 6
thres <- 0.6626053

readRDS("./Models/nb_model.rda") -> nb_model
theta <- nb_model$theta # dispersion parameter from negative binomial model


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

