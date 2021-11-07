suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Amelia"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("future"))

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)
set.seed(12345)

load("./prebrms_ITT.RData")

fit_ITT_BG <- brm_multiple(BG~FullGroup*TimeNom+(1|AnimalID),
                         data = imps_ITT,
                         iter=8000,
			 warmup=1000,
			 thin=5,
			 family="skew_normal",
			 prior=c(set_prior("normal(10,5)", class = "Intercept"),
                                    set_prior("normal(0,5)", class = "b")),
                         control = list(max_treedepth = 12,
			 adapt_delta=0.99),
			 chains=4,
                         future=T)

ITT_BG_summary<-summary(fit_ITT_BG)

save(ITT_BG_summary,fit_ITT_BG,file="fit_ITT_BG.RData")
