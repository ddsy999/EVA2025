library(evmix)
library(tidyverse)
library(zoo)
library(POT)
library(extRemes)
library(texmex)
library(stabledist) #positive stable disribution

rps<-function(n,alpha){
  #### PS(alpha) generation as given by Stephenson(2003)
  unif <- runif(n)*pi
  stdexp.ps <- rexp(n,1)
  logs <- (1-alpha)/alpha * log(sin((1-alpha)*unif)) + 
    log(sin(alpha*unif)) - (1-alpha)/alpha * log(stdexp.ps) - 
    1/alpha * log(sin(unif))
  return(exp(logs))}

qGPD<-function(tau,B,thresh){
  
  if(is.matrix(B)){
    prob<-B[,1]
    sig<-B[,2]
    xi<-B[,3]
  }
  
  if(!is.matrix(B)){
    prob<-B[1]
    sig<-B[2]
    xi<-B[3]
  }
  
  tau2<-1-(tau-prob)/(1-prob)
  Y<-thresh+ifelse(tau<prob,0,sig*(tau2^(-xi)-1)/xi)
  Y}


workingDir = "C:/Users/ddsy9/Documents/GitHub/EVA2025/"
source(paste0(workingDir,"eva-reich2014-MCMC.R"))
MCMCResult_winter_dataNC <- readRDS(paste0(workingDir,"MCMC_result/MCMCsamples_winter_dataNC.RDS"))
MCMCResult_winter_dataCO <- readRDS(paste0(workingDir,"MCMC_result/MCMCsamples_winter_dataCO.RDS"))
MCMCResult_summer_dataNC <- readRDS(paste0(workingDir,"MCMC_result/MCMCsamples_summer_dataNC.RDS"))
MCMCResult_summer_dataCO <- readRDS(paste0(workingDir,"MCMC_result/MCMCsamples_summer_dataCO.RDS"))
# MCMCResult_winter_dataOri <- readRDS(paste0(workingDir,"MCMC_result/MCMCsamples_winter_dataOriginal.RDS"))
MCMCResult_summer_dataOri <- readRDS(paste0(workingDir,"MCMC_result/MCMCsamples_summer_dataOriginal.RDS"))


run1_CO <- readRDS(paste0(workingDir,"Validation/CO/run1011.001.rds"))
run2_CO <- readRDS(paste0(workingDir,"Validation/CO/run1231.013.rds"))
run3_CO <- readRDS(paste0(workingDir,"Validation/CO/run1251.011.rds"))
run4_CO <- readRDS(paste0(workingDir,"Validation/CO/run1251.019.rds"))

run1_NC <- readRDS(paste0(workingDir,"Validation/NC/run1031.002.rds"))
run2_NC <- readRDS(paste0(workingDir,"Validation/NC/run1111.006.rds"))
run3_NC <- readRDS(paste0(workingDir,"Validation/NC/run1231.011.rds"))
run4_NC <- readRDS(paste0(workingDir,"Validation/NC/run1231.020.rds"))

load(paste0(workingDir,"data/run1.RData"))
load(paste0(workingDir,"data/run2.RData"))
load(paste0(workingDir,"data/run3.RData"))
load(paste0(workingDir,"data/run4.RData"))



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


y_Ori <- array(0, dim = c(60225, 5, 5, 4))
for(t in 1:60225){
  y_Ori[t, , ,1] <- run1[,,t]
  y_Ori[t, , ,2] <- run2[,,t]
  y_Ori[t, , ,3] <- run3[,,t]
  y_Ori[t, , ,4] <- run4[,,t]
}
y_Ori_wint <- y_CO[wint_vec,,,]
y_Ori_summ <- y_CO[summ_vec,,,]



################## Data Change ###########################
sourceData = y_Ori_summ 
samplePriorDistn = MCMCResult_summer_dataOri 
sourceData = sourceData[[1]]
samplePriorDistn = samplePriorDistn[[1]]
###########################################################


####### Setting ################
M = nrow(samplePriorDistn)
# L = 25
# nd = 25
L = 5
nd = 5


####### Treshold
# u = numeric(25)
u_univarite <- quantile(sourceData, prob = .9)
u = rep(u_univarite,25) %>% as.vector()

####### dw2
if(nd == 5){ 
  dw2<-as.matrix(dist(1:nd,diag=T,upper=T))^2
}else if(nd == 25){
  nd_sqrt <- sqrt(nd)
  dw2<-as.matrix(dist(expand.grid(1:nd_sqrt, 1:nd_sqrt),diag=T,upper=T))^2
} 

## Y만드는 함수
X2Y <- function(X, u, pi, sigma, xi){
  U = exp(-1/X)
  if(U<= pi){result = u
  }else{result = qgpd(p=(1-U)/(1-pi),sigma=sigma , xi=xi, u = u)
  }
  return(result)
}


Y_SampleResult = NULL
maxIter = 60225*100

prelim_prob1 = numeric(1)
prelim_prob2 = numeric(1)
prelim_prob3 = numeric(1)
prelim_prob3_buffer = 0

final_prob1 = numeric(1)
final_prob2 = numeric(1)
final_prob3 = numeric(1)
final_prob3_buffer = 0

######## Iter 
for( iter in 1:maxIter){
  #####################################  
  ######## sampling prior Distn
  # set.seed(Sys.time())
  # rm(.Random.seed,envir=globalenv())
  indx = sample(1:M,1)
  samplePrior <- samplePriorDistn[indx,]
  alpha = expit(samplePrior["Alpha"])
  xi = samplePrior["Shape"]
  gamma = exp(samplePrior["BW"])
  pis = expit(samplePrior[1:25])
  sigmas = exp(samplePrior[26:50])

  
  ######## B를 만든다 
  # Bstable <- rstable(L, alpha = alpha, beta = 1, pm = 1) #beta = 1 to set positive stable dist
  Bstable <-   rps(L, alpha = alpha) #beta = 1 to set positive stable dist

  ######## A를 만든다
  FAC<-fac2FAC(make.fac(dw2,gamma)) # gaussian kernels  K(j)s
  A<-a2A(FAC,Bstable,alpha) 
  
  ######## X를 만든다 
  X <- numeric(nd) 
  for( j in 1:nd){
    X[j] = rGEV(1,A[j]^alpha, alpha*A[j]^alpha, alpha)
  }
  
  ######## Y를 만든다
  Y <- numeric(nd)
  for (j in 1:nd){
    Y[j] = X2Y(X[j], u[j], pis[j], sigmas[j], xi)
  }
  
  
  #############################
  # Prelim 
  if(nd==5){
    if(sum(Y)>85) {prelim_prob1 = prelim_prob1 + 1}
    if(sort(Y,decreasing=T)[3]>4.3) {prelim_prob2 = prelim_prob2 + 1}
    if(sort(Y,decreasing=T)[3]>2.5) {
      if(prelim_prob3_buffer==1){
      prelim_prob3 = prelim_prob3 + 1
      }
      prelim_prob3_buffer = 1 
    }else{
      prelim_prob3_buffer = 0
    }
  }
  # Final 
  if(nd==25){
    if(sum(Y)>1.7) {final_prob1 = final_prob1 + 1}
    if(sort(Y,decreasing=T)[6]>5.7) {final_prob2 = final_prob2 + 1}
    if(sort(Y,decreasing=T)[3]>5) {
      if(final_prob3_buffer==1){
        final_prob3 = final_prob3 + 1
      }
      final_prob3_buffer = 1 
    }else{
      final_prob3_buffer = 0
    }
  }
  
}

prelim_prob1/maxIter*100
prelim_prob2/maxIter*100
prelim_prob3/maxIter*100



