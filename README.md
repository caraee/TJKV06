# TJKV06
This repository contains R code for performing multiple imputations and Bayesian
analyses on all experiments related to animal experiments in this study. 

Missing data were imputed using the Amelia package, which has excellent help. See
https://gking.harvard.edu/amelia for examples and help. Priors based on the 
limit of detection of assays or the glucometer were used were appropriate. An
arbitrary high value (147 mM) was chosen for the maximum for blood glucose. 

Imputations were inspected to ensure that the predicted values were reasonable.

The brms package was used as an interface with Stan. An excellent starting point
is https://paul-buerkner.github.io/brms/. Most models were run using a compute 
cluster (https://www.computecanada.ca/) as they are time consuming to run. No
model had issues with convergence or divergence. Note that for multiple
imputations, the Rhat values reported in the model summary do not reflect the
Rhats for the individual imputations - see 
https://cran.r-project.org/web/packages/brms/vignettes/brms_multivariate.html.

Where possible/depending on when I wrote it, code uses tidyverse conventions.
