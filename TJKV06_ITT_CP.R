suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Amelia"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("future"))

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)
set.seed(12345)

load("./prebrms_ITT.RData")

fit_ITT_CP <- brm_multiple(bf(CPnM~FullGroup*TimeNom+(1|AnimalID),
sigma~FullGroup), 
                           iter=8000,
                           warmup=2000,
                           thin=5,
                           family="skew_normal",
                           prior=c(set_prior("normal(1,5)", class = "Intercept"),
                                   set_prior("normal(0,5)", class = "b")),
                           control = list(adapt_delta=0.99,
                                          max_treedepth=15),
                           data = imps_ITT_CP,
                           #chains=1,
                           #future=T,
                           silent=0,
                           backend = "cmdstanr",
                           threads = threading(parallel::detectCores()),
                           save_pars = save_pars(all = TRUE))

ITT_CP_summary<-summary(fit_ITT_CP)

save(ITT_CP_summary,fit_ITT_CP,file="fit_ITT_CP2.RData")
