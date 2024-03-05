CV <- "CV1"


# Loading libraries:
library(rlist)
library(tictoc)
library(MegaLMM)
library(gfBLUP)
library(doParallel)

# Setting seed:
set.seed(1997)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading kinship and marker data:
load("genotypes/K_hyper.RData")

# Hyperspectral data:
tic("MegaLMM CV1")

datasets <- 1:250
n.datasets <- length(datasets)
n.cores <- 10
work <- split(datasets, ceiling(seq_along(datasets) / ceiling(n.datasets / n.cores)))
cl <- parallel::makeCluster(n.cores, outfile = sprintf("logs/MegaLMM_CV1_hyper_datasets_%s_%s.txt", datasets[1], datasets[n.datasets]))
doParallel::registerDoParallel(cl)

invisible(
  par.results <- foreach::foreach(k = 1:length(work), .packages = c("MegaLMM", "gfBLUP", "rlist", "tictoc"), .combine = "c") %dopar% {
    
    par.work <- work[[k]]
    set.seed(1997)
    
    # Setting up result storage:
    acc <- numeric(length(par.work))
    
    # Running:
    for (run in 1:length(par.work)) {
      
      # Loading hyperspectral dataset:
      datalist <- list.load(file = sprintf("hyper/datasets/hyper_dataset_%d.RData", par.work[run]))
      
      # Storing data and prediction target:
      d <- datalist$data
      pred.target <- datalist$pred.target
      test.set <- datalist$test.set
      train.set <- datalist$train.set
      d.train <- droplevels(d[which(!is.na(d$Y)), ])
      d.test <- droplevels(d[which(is.na(d$Y)), ])
      
      # Calculating genotypic means (BLUEs):
      d.train <- gfBLUP:::genotypeMeans(d.train)
      d.test <- gfBLUP:::genotypeMeans(d.test)
      
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
      tic(par.work[run])
      
      # MegaLMM config:
      run_parameters <- MegaLMM::MegaLMM_control(
        drop0_tol = 1e-10,
        scale_Y = FALSE,
        h2_divisions = 20,
        h2_step_size = NULL,
        burn = 0,
        K = 100,
        save_current_state = TRUE,
        thin = 2
      )
      
      priors = MegaLMM::MegaLMM_priors(
        tot_Y_var = list(V = 0.5, nu = 10),
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
      run_ID <- sprintf("hyper/megalmm_states/%s_hyper_dataset_%d", CV, par.work[run])
      
      # Initializing MegaLMM:
      MegaLMM_state = MegaLMM::setup_model_MegaLMM(d[, 2:ncol(d)],
                                                   ~ 1 + (1|G),
                                                   data = d,
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
        MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state, drop_cor_threshold = 0.6)
        MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
        MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
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
      ########################################################################
      
      acc[run] <- cor(pred.target$pred.target, mean_pred.test)
      
      # Deleting MegaLMM state files:
      unlink(run_ID, recursive = TRUE)
    }
    
    # Retrieve computational times:
    tictoc.logs <- tic.log(format = FALSE)
    tic.clearlog()
    comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
    
    # Collect results:
    worker.result <- list(data.frame(acc = acc,
                                     comptimes = comptimes))
    
    names(worker.result) <- sprintf("worker_%d", k)
    return(worker.result)
    
  })
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
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
write.csv(results, sprintf("hyper/results/12%s_hyper_results_megalmm_%s.csv", lab, CV))


