# Loading libraries:
library(rlist)
library(tictoc)
library(MegaLMM)
library(gfBLUP)
library(doParallel)

# Setting seed:
set.seed(1997)

# Setting working directory:
# wd <- "C:/Users/Killian/Desktop/gfblup-methodological-paper"
wd <- "~/gfblup_methodology"
setwd(wd)

# Loading kinship and marker data:
load("K_hyper.RData")



# Part 1: datasets 1-20 for CV1 ====================================================================================================================
set.seed(1997)
CV <- "CV1"
datasets <- 1:10

n.threads <- 10
work <- split(datasets, ceiling(seq_along(datasets) / ceiling(length(datasets) / n.threads)))

cat(sprintf("Running MegaLMM on dataset %d to %d (%s, total %d) in parallel using %d threads\n", datasets[1], datasets[length(datasets)], CV, length(datasets), n.threads))
tic(sprintf("MegaLMM %s, dataset %d to %d", CV, datasets[1], datasets[length(datasets)]))

doParallel::registerDoParallel(cores = n.threads)
invisible(
par.results <- foreach::foreach(worker = 1:length(work), .combine = "c", .packages = c("MegaLMM", "gfBLUP", "rlist", "tictoc")) %dopar% {
  
  par.work <- work[[worker]]
  
  # Setting up result storage:
  acc <- numeric(length(par.work))
  
  # Running:
  for (run in 1:length(par.work)) {
    
    # Loading simulated dataset:
    datalist <- rlist::list.load(file = sprintf("analyses_datasets/hyper/hyper_dataset_%d.RData", par.work[run]))
    
    # Storing data and prediction target:
    d <- datalist$data
    pred.target <- datalist$pred.target
    test.set <- datalist$test.set
    train.set <- datalist$train.set
    d.train <- droplevels(d[which(!is.na(d$Y)), ])
    d.test <- droplevels(d[which(is.na(d$Y)), ])
    
    # Redundancy filtering (0.95):
    temp <- gfBLUP:::redundancyFilter(data = na.omit(d)[, 1:(ncol(d) - 1)], tau = 0.95, verbose = FALSE)
    d.train.RF <- cbind(temp$data.RF, d.train["Y"])
    d.test.RF <- d.test[names(d.train.RF)]
    
    # Calculating genotypic means (BLUEs):
    d.train.RF <- gfBLUP:::genotypeMeans(d.train.RF)
    d.test.RF <- gfBLUP:::genotypeMeans(d.test.RF)
    
    # Rescaling now we only have means:
    d.train.RF[, 2:ncol(d.train.RF)] <- sapply(d.train.RF[, 2:ncol(d.train.RF)], scale)
    d.test.RF[, 2:ncol(d.test.RF)] <- sapply(d.test.RF[, 2:ncol(d.test.RF)], scale)
    
    d.RF <- rbind(d.test.RF, d.train.RF)
    
    # Reordering the columns (G, Y, SEC):
    d.RF <- cbind(d.RF[, 1], d.RF[, ncol(d.RF)], d.RF[, (2:(ncol(d.RF) - 1))])
    names(d.RF)[1:2] <- c("G", "Y")
    
    # Setting test set secondary features to NA if in CV1:
    if (CV == "CV1") {
      d.RF[which(is.na(d.RF$Y)), 3:ncol(d)] <- NA
    }
    
    # Start of timing:
    tic(run)
    
    # MegaLMM config:
    run_parameters <- MegaLMM::MegaLMM_control(
      drop0_tol = 1e-10,
      scale_Y = FALSE,
      h2_divisions = 20,
      h2_step_size = NULL,
      burn = 0,
      K = 20,
      save_current_state = TRUE,
      thin = 2
    )
    
    priors = MegaLMM::MegaLMM_priors(
      tot_Y_var = list(V = 0.5, nu = 10),
      tot_F_var = list(V = 18/20, nu = 100000),
      Lambda_prior = list(
        sampler = MegaLMM::sample_Lambda_prec_horseshoe,
        prop_0 = 0.1,
        delta = list(shape = 3, scale = 1),
        delta_iterations_factor = 100
      ),
      h2_priors_resids_fun = function(h2s, n) 1,
      h2_priors_factors_fun = function(h2s, n) 1
    )
    
    # Doing the redundancy filtering again within the tictoc measurement so it is also included in the comptime.
    # Note that we don't do anything with the output. We only redo it because the scaling should not be included
    # in the comptimes.
    # Redundancy filtering (0.95):
    temp_timing <- gfBLUP:::redundancyFilter(data = na.omit(d)[, 1:(ncol(d) - 1)], tau = 0.95, verbose = TRUE)
    d.train.RF_timing <- cbind(temp_timing$data.RF, d.train["Y"])
    d.test.RF_timing <- d.test[names(d.train.RF_timing)]
    
    run_ID <- sprintf("megalmm_states/hyper/%s_hyper_dataset_%d", CV, par.work[run])
    
    # Initializing MegaLMM:
    MegaLMM_state = MegaLMM::setup_model_MegaLMM(d.RF[, 2:ncol(d.RF)],
                                                 ~ 1 + (1|G),
                                                 data = d.RF,
                                                 relmat = list(G = K),
                                                 run_parameters = run_parameters,
                                                 run_ID = run_ID)
    
    maps = MegaLMM::make_Missing_data_map(MegaLMM_state, verbose = FALSE)
    MegaLMM_state <- MegaLMM::set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)
    
    MegaLMM_state <- MegaLMM::set_priors_MegaLMM(MegaLMM_state, priors)
    MegaLMM_state <- MegaLMM::initialize_variables_MegaLMM(MegaLMM_state)
    MegaLMM_state <- MegaLMM::initialize_MegaLMM(MegaLMM_state, verbose = FALSE)
    
    MegaLMM_state$Posterior$posteriorSample_params <- c("Lambda")
    MegaLMM_state$Posterior$posteriorFunctions <- list(pred = "U_R[,1] + U_F %*% Lambda[,1]")
    
    if (file.exists(run_ID)) {
      MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
    }
    
    # n_threads = optimize_n_threads(MegaLMM_state, seq(1, RcppParallel::defaultNumThreads(), by = 1), times = 2)
    # set_MegaLMM_nthreads(n_threads$optim)
    
    # Burn-in, collect 200 samples, reorder factors, clear posterior, repeat
    # (10 times for a total of 2000 burn-in samples):
    # Don't make trace plots now, we'll check whether the chains have become stationary
    # after we finish collecting all posterior samples.
    #
    # NO PRINTING THE PROGRESS IN PARALLEL RUNS!
    n_iter <- 100
    n_burn_in <- 10
    for (i in 1:n_burn_in) {
      MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state, drop_cor_threshold = 0.6)
      MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
      MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
    }
    
    # Clearing the burn-in samples:
    MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
    
    # Collecting the posterior samples:
    # 2000 iterations with a thinning rate of 2 gives 1000 posterior samples.
    #
    # So 2000 burn-in samples and 1000 posterior samples.
    n_iter <- 2000
    n_sampling <- 1
    for (i in 1:n_sampling) {
      MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
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
    
    acc[run] <- cor(pred.target$pred.target, mean_pred.test)
    
    # Deleting MegaLMM state files:
    unlink(run_ID, recursive = TRUE)
    
  }
  
  # Retrieve computational times:
  tictoc.logs <- tic.log(format = FALSE)
  tic.clearlog()
  comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
  
  # Collect results:
  result <- list(data.frame(acc = acc,
                            comptimes = comptimes))
  
  names(result) <- sprintf("worker_%d", worker)
  return(result)
  
})
doParallel::stopImplicitCluster()
toc()

# Restructuring parallel results for saving:
acc <- numeric()
comptimes <- numeric()
for (j in 1:length(work)) {
  
  acc <- c(acc, par.results[[sprintf("worker_%d", j)]]$acc)
  comptimes <- c(comptimes, par.results[[sprintf("worker_%d", j)]]$comptimes)
  
}

results <- data.frame(acc = acc,
                      comptimes = comptimes)

# Making correct CV label:
if (CV == "CV1") {
  lab <- "a"
} else if (CV == "CV2") {
  lab <- "b"
}

# Export results:
write.csv(results, sprintf("analyses_results/hyper/12%s_hyper_results_megalmm_d%d_d%d.csv", lab, datasets[1], datasets[length(datasets)]))



# rm(list=ls())
# # Setting seed:
# set.seed(1997)
# 
# # Setting working directory:
# # wd <- "C:/Users/Killian/Desktop/gfblup-methodological-paper"
# wd <- "~/gfblup_methodology"
# setwd(wd)
# 
# # Loading kinship and marker data:
# load("K_hyper.RData")
# 
# # Part 2: datasets 1-250 for CV2 ====================================================================================================================
# set.seed(1997)
# CV <- "CV2"
# datasets <- 1:24
# 
# n.threads <- 12
# work <- split(datasets, ceiling(seq_along(datasets) / ceiling(length(datasets) / n.threads)))
# 
# cat(sprintf("Running MegaLMM on dataset %d to %d (%s, total %d) in parallel using %d threads\n", datasets[1], datasets[length(datasets)], CV, length(datasets), n.threads))
# tic(sprintf("MegaLMM %s, dataset %d to %d", CV, datasets[1], datasets[length(datasets)]))
# 
# doParallel::registerDoParallel(cores = n.threads)
# invisible(
#   par.results <- foreach::foreach(worker = 1:length(work), .combine = "c", .packages = c("MegaLMM", "gfBLUP", "rlist", "tictoc")) %dopar% {
#     
#     par.work <- work[[worker]]
#     
#     # Setting up result storage:
#     acc <- numeric(length(par.work))
#     
#     # Running:
#     for (run in 1:length(par.work)) {
#       
#       # Loading simulated dataset:
#       datalist <- rlist::list.load(file = sprintf("analyses_datasets/hyper/hyper_dataset_%d.RData", par.work[run]))
#       
#       # Storing data and prediction target:
#       d <- datalist$data
#       pred.target <- datalist$pred.target
#       test.set <- datalist$test.set
#       train.set <- datalist$train.set
#       d.train <- droplevels(d[which(!is.na(d$Y)), ])
#       d.test <- droplevels(d[which(is.na(d$Y)), ])
#       
#       # Redundancy filtering (0.95):
#       temp <- gfBLUP:::redundancyFilter(data = na.omit(d)[, 1:(ncol(d) - 1)], tau = 0.95, verbose = FALSE)
#       d.train.RF <- cbind(temp$data.RF, d.train["Y"])
#       d.test.RF <- d.test[names(d.train.RF)]
#       
#       # Calculating genotypic means (BLUEs):
#       d.train.RF <- gfBLUP:::genotypeMeans(d.train.RF)
#       d.test.RF <- gfBLUP:::genotypeMeans(d.test.RF)
#       
#       # Rescaling now we only have means:
#       d.train.RF[, 2:ncol(d.train.RF)] <- sapply(d.train.RF[, 2:ncol(d.train.RF)], scale)
#       d.test.RF[, 2:ncol(d.test.RF)] <- sapply(d.test.RF[, 2:ncol(d.test.RF)], scale)
#       
#       d.RF <- rbind(d.test.RF, d.train.RF)
#       
#       # Reordering the columns (G, Y, SEC):
#       d.RF <- cbind(d.RF[, 1], d.RF[, ncol(d.RF)], d.RF[, (2:(ncol(d.RF) - 1))])
#       names(d.RF)[1:2] <- c("G", "Y")
#       
#       # Setting test set secondary features to NA if in CV1:
#       if (CV == "CV1") {
#         d.RF[which(is.na(d.RF$Y)), 3:ncol(d)] <- NA
#       }
#       
#       # Start of timing:
#       tic(run)
#       
#       # MegaLMM config:
#       run_parameters <- MegaLMM::MegaLMM_control(
#         drop0_tol = 1e-10,
#         scale_Y = FALSE,
#         h2_divisions = 20,
#         h2_step_size = NULL,
#         burn = 0,
#         K = 20,
#         save_current_state = TRUE,
#         thin = 2
#       )
#       
#       priors = MegaLMM::MegaLMM_priors(
#         tot_Y_var = list(V = 0.5, nu = 10),
#         tot_F_var = list(V = 18/20, nu = 100000),
#         Lambda_prior = list(
#           sampler = MegaLMM::sample_Lambda_prec_horseshoe,
#           prop_0 = 0.1,
#           delta = list(shape = 3, scale = 1),
#           delta_iterations_factor = 100
#         ),
#         h2_priors_resids_fun = function(h2s, n) 1,
#         h2_priors_factors_fun = function(h2s, n) 1
#       )
#       
#       # Doing the redundancy filtering again within the tictoc measurement so it is also included in the comptime.
#       # Note that we don't do anything with the output. We only redo it because the scaling should not be included
#       # in the comptimes.
#       # Redundancy filtering (0.95):
#       temp_timing <- gfBLUP:::redundancyFilter(data = na.omit(d)[, 1:(ncol(d) - 1)], tau = 0.95, verbose = TRUE)
#       d.train.RF_timing <- cbind(temp_timing$data.RF, d.train["Y"])
#       d.test.RF_timing <- d.test[names(d.train.RF_timing)]
#       
#       run_ID <- sprintf("megalmm_states/hyper/%s_hyper_dataset_%d", CV, par.work[run])
#       
#       # Initializing MegaLMM:
#       MegaLMM_state = MegaLMM::setup_model_MegaLMM(d.RF[, 2:ncol(d.RF)],
#                                                    ~ 1 + (1|G),
#                                                    data = d.RF,
#                                                    relmat = list(G = K),
#                                                    run_parameters = run_parameters,
#                                                    run_ID = run_ID)
#       
#       maps = MegaLMM::make_Missing_data_map(MegaLMM_state, verbose = FALSE)
#       MegaLMM_state <- MegaLMM::set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)
#       
#       MegaLMM_state <- MegaLMM::set_priors_MegaLMM(MegaLMM_state, priors)
#       MegaLMM_state <- MegaLMM::initialize_variables_MegaLMM(MegaLMM_state)
#       MegaLMM_state <- MegaLMM::initialize_MegaLMM(MegaLMM_state, verbose = FALSE)
#       
#       MegaLMM_state$Posterior$posteriorSample_params <- c("Lambda")
#       MegaLMM_state$Posterior$posteriorFunctions <- list(pred = "U_R[,1] + U_F %*% Lambda[,1]")
#       
#       if (file.exists(run_ID)) {
#         MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
#       }
#       
#       # n_threads = optimize_n_threads(MegaLMM_state, seq(1, RcppParallel::defaultNumThreads(), by = 1), times = 2)
#       # set_MegaLMM_nthreads(n_threads$optim)
#       
#       # Burn-in, collect 200 samples, reorder factors, clear posterior, repeat
#       # (10 times for a total of 2000 burn-in samples):
#       # Don't make trace plots now, we'll check whether the chains have become stationary
#       # after we finish collecting all posterior samples.
#       #
#       # NO PRINTING THE PROGRESS IN PARALLEL RUNS!
#       n_iter <- 100
#       n_burn_in <- 10
#       for (i in 1:n_burn_in) {
#         MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state, drop_cor_threshold = 0.6)
#         MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
#         MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
#       }
#       
#       # Clearing the burn-in samples:
#       MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
#       
#       # Collecting the posterior samples:
#       # 2000 iterations with a thinning rate of 2 gives 1000 posterior samples.
#       #
#       # So 2000 burn-in samples and 1000 posterior samples.
#       n_iter <- 2000
#       n_sampling <- 1
#       for (i in 1:n_sampling) {
#         MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
#         MegaLMM_state <- MegaLMM::save_posterior_chunk(MegaLMM_state)
#       }
#       
#       # Reloading the saved posterior samples:
#       Lambda_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "Lambda")
#       pred_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "pred")
#       
#       # Calculating means:
#       mean_Lambda <- MegaLMM::get_posterior_mean(Lambda_samples)
#       mean_pred <- MegaLMM::get_posterior_mean(pred_samples)
#       
#       # Storing genotype names in the order of mean_pred:
#       names <- gsub("([0-9]*)::G", "\\1", rownames(mean_pred))
#       mean_pred <- as.numeric(mean_pred)
#       names(mean_pred) <- names
#       mean_pred.test <- mean_pred[which(names(mean_pred) %in% test.set)]
#       mean_pred.test <- mean_pred.test[match(pred.target$G, names(mean_pred.test))]
#       
#       toc(log = TRUE)
#       
#       acc[run] <- cor(pred.target$pred.target, mean_pred.test)
#       
#       # Deleting MegaLMM state files:
#       unlink(run_ID, recursive = TRUE)
#       
#     }
#     
#     # Retrieve computational times:
#     tictoc.logs <- tic.log(format = FALSE)
#     tic.clearlog()
#     comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
#     
#     # Collect results:
#     result <- list(data.frame(acc = acc,
#                               comptimes = comptimes))
#     
#     names(result) <- sprintf("worker_%d", worker)
#     return(result)
#     
#   })
# doParallel::stopImplicitCluster()
# toc()
# 
# # Restructuring parallel results for saving:
# acc <- numeric()
# comptimes <- numeric()
# for (j in 1:length(work)) {
#   
#   acc <- c(acc, par.results[[sprintf("worker_%d", j)]]$acc)
#   comptimes <- c(comptimes, par.results[[sprintf("worker_%d", j)]]$comptimes)
#   
# }
# 
# results <- data.frame(acc = acc,
#                       comptimes = comptimes)
# 
# # Making correct CV label:
# if (CV == "CV1") {
#   lab <- "a"
# } else if (CV == "CV2") {
#   lab <- "b"
# }
# 
# # Export results:
# write.csv(results, sprintf("analyses_results/hyper/12%s_hyper_results_megalmm_d%d_d%d.csv", lab, datasets[1], datasets[length(datasets)]))

