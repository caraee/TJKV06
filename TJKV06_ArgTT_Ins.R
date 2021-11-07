suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Amelia"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("future"))

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)
set.seed(12345)

load("./prebrms_ArgTT.RData")

fit_ArgTT_Ins <- brm_multiple(bf(InsnM~FullGroup*TimeNom+(1|AnimalID),
sigma~FullGroup),
                          iter=8000,
			  warmup=1000,
			  thin=5,
                          data = imps_ArgTT_Horm,
			  family="skew_normal",
                          #prior = set_prior("normal(0,1)", class = "b"),
                          control = list(adapt_delta = 0.9,
                                         max_treedepth=12),
                          #chains = 4,
                          silent=0)

save(fit_ArgTT_Ins,file="./fit_ArgTT_Ins.RData")
