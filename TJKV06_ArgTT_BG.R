suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Amelia"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("future"))

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)
set.seed(12345)

load("./prebrms_ArgTT.RData")

fit_ArgTT_BG <- brm_multiple(BG~FullGroup*TimeNom+(TimeNom|AnimalID),
                             data = imps_ArgTT,
                             iter=8000,
			     thin=5,
			     warmup=1000,
			     prior=c(set_prior("normal(10,3)", class = "Intercept"),
                                   set_prior("normal(0,5)", class = "b")),
			    #family = "skew_normal",
                             control = list(adapt_delta=0.95,
			     max_treedepth = 12),
                             future=T)

save(fit_ArgTT_BG,file="./fit_ArgTT_BG.RData")
