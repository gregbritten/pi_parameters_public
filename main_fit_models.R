library(randomForest)
library(lme4)
source('r/main_process_nc.R')

ndays <- 100
FITS=RSQ <- list()

for(i in 1:ndays){
  print(i)
  d <- D_nc[[i]] %>% filter(complete.cases(sst,par,depth,pico,region),depth<=mld)
  
  fit_g_Ek    <- randomForest(Ek    ~ sst + par + depth + pico + region, data=d, importance=TRUE) 
  fit_g_PBmax <- randomForest(PBmax ~ sst + par + depth + pico + region, data=d, importance=TRUE) 
  fit_g_alpha <- randomForest(alpha ~ sst + par + depth + pico + region, data=d, importance=TRUE) 
  
  FITS[["reg"]][["Ek"]][[i]]    <- fit_g_Ek
  FITS[["reg"]][["PBmax"]][[i]] <- fit_g_PBmax
  FITS[["reg"]][["alpha"]][[i]] <- fit_g_alpha
  RSQ[["reg"]][["Ek"]][[i]]     <- cor(fit_g_Ek$predicted,d$Ek)^2
  RSQ[["reg"]][["PBmax"]][[i]]  <- cor(fit_g_PBmax$predicted,d$PBmax)^2
  RSQ[["reg"]][["alpha"]][[i]]  <- cor(fit_g_Ek$predicted,d$alpha)^2
  
  fit_g_Ek_noreg    <- randomForest(Ek    ~ sst + par + depth + pico, data=d, importance=TRUE) 
  fit_g_PBmax_noreg <- randomForest(PBmax ~ sst + par + depth + pico, data=d, importance=TRUE) 
  fit_g_alpha_noreg <- randomForest(alpha ~ sst + par + depth + pico, data=d, importance=TRUE) 
  
  FITS[["noreg"]][["Ek"]][[i]]    <- fit_g_Ek_noreg
  FITS[["noreg"]][["PBmax"]][[i]] <- fit_g_PBmax_noreg
  FITS[["noreg"]][["alpha"]][[i]] <- fit_g_alpha_noreg
  RSQ[["noreg"]][["Ek"]][[i]]     <- cor(fit_g_Ek_noreg$predicted,d$Ek)^2
  RSQ[["noreg"]][["PBmax"]][[i]]  <- cor(fit_g_PBmax_noreg$predicted,d$PBmax)^2
  RSQ[["noreg"]][["alpha"]][[i]]  <- cor(fit_g_alpha_noreg$predicted,d$alpha)^2
  
  dd <- D_nc[[i]] %>% filter(complete.cases(sst,par,depth,pico),depth<=mld)
  
  fit_g_Ek_noreg_dd    <- randomForest(Ek    ~ sst + par + depth + pico, data=dd, importance=TRUE) 
  fit_g_PBmax_noreg_dd <- randomForest(PBmax ~ sst + par + depth + pico, data=dd, importance=TRUE) 
  fit_g_alpha_noreg_dd <- randomForest(alpha ~ sst + par + depth + pico, data=dd, importance=TRUE) 
  
  FITS[["noreg_dd"]][["Ek"]][[i]]    <- fit_g_Ek_noreg_dd
  FITS[["noreg_dd"]][["PBmax"]][[i]] <- fit_g_PBmax_noreg_dd
  FITS[["noreg_dd"]][["alpha"]][[i]] <- fit_g_alpha_noreg_dd
  RSQ[["noreg_dd"]][["Ek"]][[i]]     <- cor(fit_g_Ek_noreg_dd$predicted,dd$Ek)^2
  RSQ[["noreg_dd"]][["PBmax"]][[i]]  <- cor(fit_g_PBmax_noreg_dd$predicted,dd$PBmax)^2
  RSQ[["noreg_dd"]][["alpha"]][[i]]  <- cor(fit_g_alpha_noreg_dd$predicted,dd$alpha)^2
  
  d <- D_nc[[i]] %>% filter(complete.cases(sst,par,depth,pico,region),depth<=40)
  fit_g_Ek_lm    <- lmer(Ek    ~ sst + par + depth + pico + (1|region), data=d) 
  fit_g_PBmax_lm <- lmer(PBmax ~ sst + par + depth + pico + (1|region), data=d) 
  fit_g_alpha_lm <- lmer(alpha ~ sst + par + depth + pico + (1|region), data=d) 
  
  FITS[["lm"]][["Ek"]][[i]]    <- fit_g_Ek_lm
  FITS[["lm"]][["PBmax"]][[i]] <- fit_g_PBmax_lm
  FITS[["lm"]][["alpha"]][[i]] <- fit_g_alpha_lm
  RSQ[["lm"]][["Ek"]][[i]]     <- cor(predict(fit_g_Ek_lm),d$Ek)^2
  RSQ[["lm"]][["PBmax"]][[i]]  <- cor(predict(fit_g_PBmax_lm),d$PBmax)^2
  RSQ[["lm"]][["alpha"]][[i]]  <- cor(predict(fit_g_alpha_lm),d$alpha)^2
}

saveRDS(file='results/FITS.rds',FITS)
saveRDS(file='results/RSQ.rds',RSQ)
