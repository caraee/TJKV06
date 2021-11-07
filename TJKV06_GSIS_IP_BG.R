suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Amelia"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("future"))

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)
set.seed(12345)

load("./prebrms.RData")

fit_GSIS_BG <- brm_multiple(BG~FullGroup*(WeeksNom/TimeNom)+
                              (WeeksNom|AnimalID),
                            iter=8000,
                            warmup=1000,
			    thin=4,
                            chains=8,
			    family="gamma",
                            prior=set_prior("normal(-0.5,5)", class = "b"),
			    control = list(max_treedepth = 12,
			    adapt_delta=0.9),
			    data = imps_IP,
			    future=T,
			    silent=0)

GSIS_BG_summary<-summary(fit_GSIS_BG)

save(fit_GSIS_BG,GSIS_BG_summary,file="./fit_GSIS_BG.RData")
