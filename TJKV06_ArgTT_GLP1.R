suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("Amelia"))
suppressPackageStartupMessages(library("brms"))
suppressPackageStartupMessages(library("future"))

options(mc.cores = parallel::detectCores())
plan(multisession)
rstan::rstan_options(auto_write = TRUE)
set.seed(12345)

load("./prebrms_ArgTT.RData")

fit_ArgTT_GLP1 <- brm_multiple(bf(GLP1pM~FullGroup*TimeNom+(TimeNom|AnimalID)),
                            iter=8000,
                            warmup=2000,
                            family="skew_normal",
			    thin=5,
                            prior = c(set_prior("normal(0,5)", class = "b"),
                                      set_prior("normal(10,5)", class = "Intercept")),
                            control = list(adapt_delta=0.95),
                            data = imps_ArgTT_Horm,
                            #chains=8,
                            silent=0,
                            backend = "cmdstanr",
                            threads = threading(parallel::detectCores()),
                            save_pars = save_pars(all = TRUE))

GLP1_summary<-summary(fit_ArgTT_GLP1)

save(fit_ArgTT_GLP1,GLP1_summary,file="./fit_ArgTT_GLP1_2.RData")

