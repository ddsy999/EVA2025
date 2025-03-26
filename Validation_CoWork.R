# Library
library(evmix)
library(tidyverse)
library(zoo)
library(POT)
library(extRemes)
library(texmex)

workingDir = "/Users/choijisoo/Documents/Github/"
# Read Source Code 
source(paste0(workingDir,"EVA2025/eva-reich2014-MCMC.R"))

run1_CO <- readRDS(paste0(workingDir,"EVA2025/Validation/CO/run1011.001.rds"))
run2_CO <- readRDS(paste0(workingDir,"EVA2025/Validation/CO/run1231.013.rds"))
run3_CO <- readRDS(paste0(workingDir,"EVA2025/Validation/CO/run1251.011.rds"))
run4_CO <- readRDS(paste0(workingDir,"EVA2025/Validation/CO/run1251.019.rds"))

run1_NC <- readRDS(paste0(workingDir,"EVA2025/Validation/NC/run1031.002.rds"))
run2_NC <- readRDS(paste0(workingDir,"EVA2025/Validation/NC/run1111.006.rds"))
run3_NC <- readRDS(paste0(workingDir,"EVA2025/Validation/NC/run1231.011.rds"))
run4_NC <- readRDS(paste0(workingDir,"EVA2025/Validation/NC/run1231.020.rds"))

month_vec <- month(as.Date(rep(c(1:365), time = 165)))
wint_vec <- c(1:60225)[month_vec %in% c(1,2,3,4,11,12)] 
summ_vec <- c(1:60225)[month_vec %in% c(5:10)]


y_NC <- array(0, dim = c(60225, 5, 5, 4))
for(t in 1:60225){
  y_NC[t, , ,1] <- run1_NC[,,t]
  y_NC[t, , ,2] <- run2_NC[,,t]
  y_NC[t, , ,3] <- run3_NC[,,t]
  y_NC[t, , ,4] <- run4_NC[,,t]
}
y_NC_wint <- y_NC[wint_vec,,,]
y_NC_summ <- y_NC[summ_vec,,,]


y_CO <- array(0, dim = c(60225, 5, 5, 4))
for(t in 1:60225){
  y_CO[t, , ,1] <- run1_CO[,,t]
  y_CO[t, , ,2] <- run2_CO[,,t]
  y_CO[t, , ,3] <- run3_CO[,,t]
  y_CO[t, , ,4] <- run4_CO[,,t]
}
y_CO_wint <- y_CO[wint_vec,,,]
y_CO_summ <- y_CO[summ_vec,,,]




# 아래 부분은 각자 주석해제 해서 사용 

#### JS ( 돌려야 하는 것 1 )
# thresh <- quantile(y_CO, prob = .9)
# mcmc_CO_wint <- prelim2_MCMC(y_CO_wint, thresh = thresh, iters = 5000, burn = 1000)
# saveRDS(mcmc_CO_wint, "MCMCsamples_winter_dataCO.RDS")


#### YW ( 돌려야 하는 것 2 )
# thresh <- quantile(y_NC, prob = .9)
# mcmc_NC_summ <- prelim2_MCMC(y_NC_summ, thresh = thresh, iters = 5000, burn = 1000)
# saveRDS(mcmc_NC_summ, "MCMCsamples_summer_dataNC.RDS")

#### JH ( 돌려야 하는 것 3 )
# thresh <- quantile(y_NC, prob = .9)
# mcmc_NC_wint <- prelim2_MCMC(y_NC_wint, thresh = thresh, iters = 5000, burn = 1000)
# saveRDS(mcmc_NC_wint, "MCMCsamples_winter_dataNC.RDS")

