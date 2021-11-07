suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Amelia"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("future"))

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)
set.seed(12345)

load("./prebrms.RData")

fit_GSIS_IP_CP <- brm_multiple(bf(CPnM ~ FullGroup * (WeeksNom/TimeNom) + (1 | AnimalID), 
                                  sigma ~ FullGroup),
                               iter=8000,
                               warmup=1000,
                               thin=5,
                               family = "skew_normal",
                               chains=8,
                               prior=c(set_prior("normal(0,5)", class = "b"),
                               set_prior("normal(0,5)", class = "Intercept")),
                               control = list(max_treedepth = 15,
                                              adapt_delta=0.9),
                               data = imps_IP_CP,
                                silent=0,
            backend = "cmdstanr",
            threads = threading(parallel::detectCores()),
            save_pars = save_pars(all = TRUE))

GSIS_IP_CP_sum<-summary(fit_GSIS_IP_CP)

save(GSIS_IP_CP_sum,fit_GSIS_IP_CP,file="./fit_GSIS_IP_CP.RData")


