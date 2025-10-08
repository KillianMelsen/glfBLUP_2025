#!/usr/bin/env Rscript
# This script produces a second, independent set of intermediate results for the 
# MegaLMM CV2VEG hyperspectral HEAT analyses.
#
# The intermediate results for MegaLMM are unfortunately not reproducible due to
# monte carlo errors, even if setting a seed. This is specific to the MegaLMM
# software and unfortunately outside of our control. The second independent set
# of intermediate results serves as a check to show that the impact of the monte
# carlo errors on the mean accuracies reported in the manuscript are negligible.

CV <- "CV2"
prep <- "nosplines" # No VEGsplines cause of a single VEG date...

# Loading libraries:
library(rlist)
library(tictoc)
library(MegaLMM)
library(glfBLUP)
source("helper_functions.R")
library(MCMCglmm)
library(coda)
library(ape)

# Setting seed:
set.seed(1997)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading kinship and marker data:
load("genotypes/K_hyper.RData")

# Hyperspectral data:
tic("MegaLMM CV2VEG")

datasets <- 1:250
n.datasets <- length(datasets)

set.seed(1997)
    
# Setting up result storage:
acc <- numeric(n.datasets)

# Running:
for (run in datasets) {
  
  cat(sprintf("Running dataset %d / %d...\n\n", run, n.datasets))
  
  # Loading hyperspectral dataset:
  datalist <- list.load(file = sprintf("hyper_1415HEAT/datasets/%s/hyper_dataset_%d.RData", prep, run))
  
  # Storing data and prediction target:
  dates <- c("150414")
  d <- datalist$data
  select <- which(substr(names(d), 7, 12) %in% dates)
  d <- d[c(1, select, ncol(d))]
  pred.target <- datalist$pred.target
  test.set <- datalist$test.set
  train.set <- datalist$train.set
  d.train <- droplevels(d[which(!is.na(d$Y)), ])
  d.test <- droplevels(d[which(is.na(d$Y)), ])
  
  # Subsetting K (only really happens for the first dataset...):
  K <- K[unique(d$G), unique(d$G)]
  
  ### Redundancy filter the secondary features using training data only ----
  sec <- names(d[2:(ncol(d) - 1)])
  foc <- names(d)[ncol(d)]
  temp <- glfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 0.95, verbose = FALSE)
  d.train <- cbind(temp$data.RF, d.train[foc])
  d.test <- d.test[colnames(d.train)]
  sec.RF <- names(d.train[2:(ncol(d.train) - 1)])

  # Calculating genotypic means (BLUEs):
  d.train <- glfBLUP:::genotypeMeans(d.train)
  d.test <- glfBLUP:::genotypeMeans(d.test)
  
  # Rescaling now we only have means:
  d.train[, 2:ncol(d.train)] <- sapply(d.train[, 2:ncol(d.train)], scale)
  d.test[, 2:ncol(d.test)] <- sapply(d.test[, 2:ncol(d.test)], scale)
  
  d <- rbind(d.test, d.train)
  
  # Reordering the columns (G, Y, SEC):
  d <- cbind(d[, 1], d[, ncol(d)], d[, (2:(ncol(d) - 1))])
  names(d)[1:2] <- c("G", "Y")
  
  # Setting test set secondary features to NA if in CV1:
  if (CV == "CV1") {
    d[which(is.na(d$Y)), 3:ncol(d)] <- NA
  }
  
  ### Model ##############################################################
  tic(run)
  
  # MegaLMM config:
  run_parameters <- MegaLMM::MegaLMM_control(
    scale_Y = FALSE,
    burn = 0,
    K = 20,
    save_current_state = TRUE,
    thin = 2
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
  run_ID <- sprintf("hyper_1415HEAT/megalmm_states/%s_%sVEG_hyper_dataset_%d_RF", prep, CV, run)
  
  # Initializing MegaLMM:
  MegaLMM_state = MegaLMM::setup_model_MegaLMM(d[, 2:ncol(d)],
                                               ~ 1 + (1|G),
                                               data = d,
                                               relmat = list(G = K),
                                               run_parameters = run_parameters,
                                               run_ID = run_ID)
  
  maps = MegaLMM::make_Missing_data_map(MegaLMM_state, verbose = TRUE)
  MegaLMM_state <- MegaLMM::set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)
  
  MegaLMM_state <- MegaLMM::set_priors_MegaLMM(MegaLMM_state, priors)
  MegaLMM_state <- MegaLMM::initialize_variables_MegaLMM(MegaLMM_state)
  MegaLMM_state <- MegaLMM::initialize_MegaLMM(MegaLMM_state, verbose = TRUE)
  
  MegaLMM_state$Posterior$posteriorSample_params <- c("Lambda")
  MegaLMM_state$Posterior$posteriorFunctions <- list(pred = "U_R[,1] + U_F %*% Lambda[,1]")
  
  # Clearing posterior samples if they already exist for some reason:
  if (file.exists(run_ID)) {
    MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
  }
  
  # Burn-in, collect 100 samples, reorder factors, clear posterior, repeat
  # (10 times for a total of 1000 burn-in samples, as in the MegaLMM paper):
  # Don't make trace plots now, we can check whether the chains have become stationary
  # after we finish collecting all posterior samples.
  #
  # NO PRINTING THE PROGRESS IN PARALLEL RUNS!
  n_iter <- 100
  n_burn_in <- 10
  for (i in 1:n_burn_in) {
    MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state)
    MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
    MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, verbose = F)
  }
  
  # Clearing the burn-in samples:
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
  
  # Collecting the posterior samples:
  # 500 iterations with a thinning rate of 2 gives 250 posterior samples.
  # As in the MegaLMM paper.
  # So 1000 burn-in samples and 250 posterior samples.
  n_iter <- 500
  n_sampling <- 1
  for (i in 1:n_sampling) {
    MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, verbose = F)
    MegaLMM_state <- MegaLMM::save_posterior_chunk(MegaLMM_state)
  }
  
  # Reloading the saved posterior samples:
  Lambda_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "Lambda")
  pred_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "pred")
  
  # Calculating means:
  mean_Lambda <- MegaLMM::get_posterior_mean(Lambda_samples)
  mean_pred <- MegaLMM::get_posterior_mean(pred_samples)
  
  # Storing genotype names in the order of mean_pred:
  names <- gsub("([0-9]*)::G", "\\1", rownames(mean_pred))
  mean_pred <- as.numeric(mean_pred)
  names(mean_pred) <- names
  mean_pred.test <- mean_pred[which(names(mean_pred) %in% test.set)]
  mean_pred.test <- mean_pred.test[match(pred.target$G, names(mean_pred.test))]
  
  toc(log = TRUE)
  ########################################################################
  
  #### Runcie & Cheng 2019 correction --------------------------------------
  temp <- estimate_gcor(data = data.frame(ID = pred.target$G,
                                          obs = pred.target$pred.target,
                                          pred = mean_pred.test),
                        Knn = K[pred.target$G, pred.target$G],
                        method = "MCMCglmm",
                        normalize = T)
  acc[run] <- temp["g_cor"]
  
  # Deleting MegaLMM state files:
  unlink(run_ID, recursive = TRUE)
}

# Retrieve computational times:
tictoc.logs <- tic.log(format = FALSE)
tic.clearlog()
comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
    
toc()

results <- data.frame(acc = acc,
                      comptimes = comptimes)

# Making correct CV label:
if (CV == "CV1") {
  lab <- "a"
} else if (CV == "CV2") {
  lab <- "b"
}

# Export results:
write.csv(results, sprintf("MegaLMM_rep_check_scripts/hyper_1415HEAT/results/%s/12%s_hyper_results_megalmm_%sVEG.csv", prep, lab, CV))


