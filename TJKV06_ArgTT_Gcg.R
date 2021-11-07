suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Amelia"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("future"))

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)
set.seed(12345)

load("./prebrms_ArgTT.RData")


fit_ArgTT_Gcg <- brm_multiple(GcgpM~FullGroup*TimeNom+(TimeNom|AnimalID),
                           iter=8000,
                           warmup=1000,
			   thin=5,
			   family="skew_normal",
                           data = imps_ArgTT_Horm,
                           control = list(adapt_delta = 0.99,
                                          max_treedepth=12),
                           chains = 4,
                           silent=0)

Gcg_summary<-summary(fit_ArgTT_Gcg)

save(Gcg_summary,fit_ArgTT_Gcg,file="./fit_ArgTT_Gcg.RData")
