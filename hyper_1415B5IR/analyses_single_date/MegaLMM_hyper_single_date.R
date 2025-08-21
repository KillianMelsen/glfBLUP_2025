# This script produces the intermediate results required to produce several figures
# of the manuscript supplementary material. If these intermediate results are not
# included in the repository, first run this script, then run the script at
# `hyper_1415B5IR/analyses_single_date/MegaLMM_hyper_single_date_plotting.R` to
# produce the figures.

# Loading libraries:
library(rlist)
library(glfBLUP)
library(ggplot2)
library(grid)

# Setting seed:
set.seed(1997)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading kinship:
load("genotypes/K_hyper.RData")

# Which dataset should we look at:
dataset <- 1

# Loading hyperspectral dataset:
datalist <- list.load(file = sprintf("hyper_1415B5IR/datasets/splines/hyper_dataset_%d.RData", dataset))

# Storing data:
d <- datalist$data

# Subsetting to 10-03-2015:
d <- d[c(1, grep(".*_150310", names(d)), ncol(d))]

# Make training data and store feature/focal trait names:
sec <- names(d[2:(ncol(d) - 1)])
foc <- names(d)[ncol(d)]

# M = 5 ========================================================================
# Some MegaLMM parameters:
M <- 5
burnin <- 10000
posterior <- 100000
thin <- 2

# Data manipulation:
d.train <- droplevels(d[which(!is.na(d$Y)), ])
d.test <- droplevels(d[which(is.na(d$Y)), ])

# Calculating genotypic means (BLUEs):
d.train <- glfBLUP:::genotypeMeans(d.train)
d.test <- glfBLUP:::genotypeMeans(d.test)

# Rescaling now we only have means:
d.train[, 2:ncol(d.train)] <- sapply(d.train[, 2:ncol(d.train)], scale)
d.test[, 2:ncol(d.test)] <- sapply(d.test[, 2:ncol(d.test)], scale)

d <- rbind(d.test, d.train)

# Reordering the columns (G, Y, SEC):
d <- cbind(d["G"], d["Y"], d[sec])
names(d)[3:ncol(d)] <- substr(names(d)[3:ncol(d)], 3, 5)

# MegaLMM config:
run_parameters <- MegaLMM::MegaLMM_control(
  scale_Y = FALSE,
  burn = 0,
  K = M,
  save_current_state = TRUE,
  thin = thin
)

priors = MegaLMM::MegaLMM_priors(
  tot_Y_var = list(V = 0.5, nu = 3),
  tot_F_var = list(V = 18/20, nu = 20),
  Lambda_prior = list(
    sampler = MegaLMM::sample_Lambda_prec_horseshoe,
    prop_0 = 0.1,
    delta = list(shape = 3, scale = 1),
    delta_iterations_factor = 100
  ),
  h2_priors_resids_fun = function(h2s, n) 1,
  h2_priors_factors_fun = function(h2s, n) 1
)

# Creating run ID:
run_ID <- "hyper_1415B5IR/megalmm_states/MegaLMM_hyper_single_date_M5"

# Initializing MegaLMM:
MegaLMM_state = MegaLMM::setup_model_MegaLMM(d[, 2:ncol(d)],
                                             ~ 1 + (1|G),
                                             data = d,
                                             relmat = list(G = K),
                                             run_parameters = run_parameters,
                                             run_ID = run_ID)

maps = MegaLMM::make_Missing_data_map(MegaLMM_state, verbose = T)
MegaLMM_state <- MegaLMM::set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)

MegaLMM_state <- MegaLMM::set_priors_MegaLMM(MegaLMM_state, priors)
MegaLMM_state <- MegaLMM::initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state <- MegaLMM::initialize_MegaLMM(MegaLMM_state, verbose = T)

MegaLMM_state$Posterior$posteriorSample_params <- c("Lambda")
MegaLMM_state$Posterior$posteriorFunctions <- list(pred = "U_R[,1] + U_F %*% Lambda[,1]")

# Clearing posterior samples if they already exist for some reason:
if (file.exists(run_ID)) {
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
}

# Burn-in, collect samples, reorder factors, clear posterior, repeat 10 times
n_iter <- burnin / 10
n_burn_in <- 10
for (i in 1:n_burn_in) {
  MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state)
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
}

# Clearing the burn-in samples:
MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)

# Collecting the posterior samples in a single go:
n_iter <- posterior
n_sampling <- 1
for (i in 1:n_sampling) {
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
  MegaLMM_state <- MegaLMM::save_posterior_chunk(MegaLMM_state)
}

# Reloading the saved posterior samples:
Lambda_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "Lambda")
pred_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "pred")

# Saving the arrays containing all posterior predictions and loadings
save(Lambda_samples, file = "hyper_1415B5IR/megalmm_hyper_single_date_arrays/LAMBDA_M5.RData")
save(pred_samples, file = "hyper_1415B5IR/megalmm_hyper_single_date_arrays/PREDS_M5.RData")

# Deleting MegaLMM state files:
unlink(run_ID, recursive = TRUE)

# M = 10 =======================================================================
# Setting seed:
set.seed(1997)
rm(Lambda_samples, pred_samples)

# Some MegaLMM parameters:
M <- 10
burnin <- 10000
posterior <- 100000
thin <- 2

# MegaLMM config:
run_parameters <- MegaLMM::MegaLMM_control(
  scale_Y = FALSE,
  burn = 0,
  K = M,
  save_current_state = TRUE,
  thin = thin
)

priors = MegaLMM::MegaLMM_priors(
  tot_Y_var = list(V = 0.5, nu = 3),
  tot_F_var = list(V = 18/20, nu = 20),
  Lambda_prior = list(
    sampler = MegaLMM::sample_Lambda_prec_horseshoe,
    prop_0 = 0.1,
    delta = list(shape = 3, scale = 1),
    delta_iterations_factor = 100
  ),
  h2_priors_resids_fun = function(h2s, n) 1,
  h2_priors_factors_fun = function(h2s, n) 1
)

# Creating run ID:
run_ID <- "hyper_1415B5IR/megalmm_states/MegaLMM_hyper_single_date_M10"

# Initializing MegaLMM:
MegaLMM_state = MegaLMM::setup_model_MegaLMM(d[, 2:ncol(d)],
                                             ~ 1 + (1|G),
                                             data = d,
                                             relmat = list(G = K),
                                             run_parameters = run_parameters,
                                             run_ID = run_ID)

maps = MegaLMM::make_Missing_data_map(MegaLMM_state, verbose = T)
MegaLMM_state <- MegaLMM::set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)

MegaLMM_state <- MegaLMM::set_priors_MegaLMM(MegaLMM_state, priors)
MegaLMM_state <- MegaLMM::initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state <- MegaLMM::initialize_MegaLMM(MegaLMM_state, verbose = T)

MegaLMM_state$Posterior$posteriorSample_params <- c("Lambda")
MegaLMM_state$Posterior$posteriorFunctions <- list(pred = "U_R[,1] + U_F %*% Lambda[,1]")

# Clearing posterior samples if they already exist for some reason:
if (file.exists(run_ID)) {
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
}

# Burn-in, collect samples, reorder factors, clear posterior, repeat 10 times
n_iter <- burnin / 10
n_burn_in <- 10
for (i in 1:n_burn_in) {
  MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state, drop_cor_threshold = 0.6)
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
}

# Clearing the burn-in samples:
MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)

# Collecting the posterior samples in a single go:
n_iter <- posterior
n_sampling <- 1
for (i in 1:n_sampling) {
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
  MegaLMM_state <- MegaLMM::save_posterior_chunk(MegaLMM_state)
}

# Reloading the saved posterior samples:
Lambda_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "Lambda")
pred_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "pred")

# Saving the arrays containing all posterior predictions and loadings
save(Lambda_samples, file = "hyper_1415B5IR/megalmm_hyper_single_date_arrays/LAMBDA_M10.RData")
save(pred_samples, file = "hyper_1415B5IR/megalmm_hyper_single_date_arrays/PREDS_M10.RData")

# Deleting MegaLMM state files:
unlink(run_ID, recursive = TRUE)

# M = 3 ========================================================================
# Setting seed:
set.seed(1997)
rm(Lambda_samples, pred_samples)

# Some MegaLMM parameters:
M <- 3
burnin <- 10000
posterior <- 100000
thin <- 2

# MegaLMM config:
run_parameters <- MegaLMM::MegaLMM_control(
  scale_Y = FALSE,
  burn = 0,
  K = M,
  save_current_state = TRUE,
  thin = thin
)

priors = MegaLMM::MegaLMM_priors(
  tot_Y_var = list(V = 0.5, nu = 3),
  tot_F_var = list(V = 18/20, nu = 20),
  Lambda_prior = list(
    sampler = MegaLMM::sample_Lambda_prec_horseshoe,
    prop_0 = 0.1,
    delta = list(shape = 3, scale = 1),
    delta_iterations_factor = 100
  ),
  h2_priors_resids_fun = function(h2s, n) 1,
  h2_priors_factors_fun = function(h2s, n) 1
)

# Creating run ID:
run_ID <- "hyper_1415B5IR/megalmm_states/MegaLMM_hyper_single_date_M3"

# Initializing MegaLMM:
MegaLMM_state = MegaLMM::setup_model_MegaLMM(d[, 2:ncol(d)],
                                             ~ 1 + (1|G),
                                             data = d,
                                             relmat = list(G = K),
                                             run_parameters = run_parameters,
                                             run_ID = run_ID)

maps = MegaLMM::make_Missing_data_map(MegaLMM_state, verbose = T)
MegaLMM_state <- MegaLMM::set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)

MegaLMM_state <- MegaLMM::set_priors_MegaLMM(MegaLMM_state, priors)
MegaLMM_state <- MegaLMM::initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state <- MegaLMM::initialize_MegaLMM(MegaLMM_state, verbose = T)

MegaLMM_state$Posterior$posteriorSample_params <- c("Lambda")
MegaLMM_state$Posterior$posteriorFunctions <- list(pred = "U_R[,1] + U_F %*% Lambda[,1]")

# Clearing posterior samples if they already exist for some reason:
if (file.exists(run_ID)) {
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
}

# Burn-in, collect samples, reorder factors, clear posterior, repeat 10 times
n_iter <- burnin / 10
n_burn_in <- 10
for (i in 1:n_burn_in) {
  MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state, drop_cor_threshold = 0.6)
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
}

# Clearing the burn-in samples:
MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)

# Collecting the posterior samples in a single go:
n_iter <- posterior
n_sampling <- 1
for (i in 1:n_sampling) {
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
  MegaLMM_state <- MegaLMM::save_posterior_chunk(MegaLMM_state)
}

# Reloading the saved posterior samples:
Lambda_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "Lambda")
pred_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "pred")

# Saving the arrays containing all posterior predictions and loadings
save(Lambda_samples, file = "hyper_1415B5IR/megalmm_hyper_single_date_arrays/LAMBDA_M3.RData")
save(pred_samples, file = "hyper_1415B5IR/megalmm_hyper_single_date_arrays/PREDS_M3.RData")

# Deleting MegaLMM state files:
unlink(run_ID, recursive = TRUE)


