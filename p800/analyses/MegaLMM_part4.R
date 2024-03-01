# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!

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
load("genotypes/K_sim.RData")

# Simulated genetic parameters:
h2.sec <- c("05", "07", "09")
comms <- c("02", "05", "08")
h2.foc <- c("01", "03", "05", "07", "09")
CVs <- c("CV1", "CV2")

combis <- length(h2.foc) * length(CVs)
combi <- 1

par.combis <- data.frame(comm = rep(comms, each = length(h2.sec)),
                         h2s = rep(h2.sec, length(comms)))

part <- "4"

# All simulated data for both CV1 and CV2:
tic("MegaLMM")
for (CV in CVs) {
  for (h2y in h2.foc) {
    
    tic(sprintf("MegaLMM p800 combi %d / %d, (CV = %s, h2y = %s)", combi, combis, CV, h2y))

    cl <- parallel::makeCluster(9, outfile = sprintf("logs/MegaLMM_sim_p800_h2y%s_%s_part%s.txt", h2y, CV, part))
    registerDoParallel(cl)
    
    invisible(
      foreach::foreach(i = 1:nrow(par.combis), .packages = c("rlist", "tictoc")) %dopar% {
        
        comm <- par.combis[i, "comm"]
        h2s <- par.combis[i, "h2s"]
        
        # Number of simulated datasets to load:
        first <- 17
        last <- 20
        n.sim <- length(first:last)
        
        # Setting up result storage:
        acc <- numeric(n.sim)

        # Running (SET SEED IN EACH PARALLEL WORKER!):
        set.seed(1997)
        run <- 1
        for (sim in first:last) {
          
          cat(sprintf("Running MegaLMM on %s p800 dataset %d / %d (part %s, dataset %s), (h2s = %s, comm = %s, h2y = %s), combi %d / %d\n",
                      CV, run, n.sim, part, sim, h2s, comm, h2y, combi, combis))
          
          # Loading simulated dataset:
          datalist <- list.load(file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
          
          # Storing data and prediction target:
          d <- datalist$data.real
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
          
          # Start of timing:
          tic(run)
          
          # MegaLMM config:
          run_parameters <- MegaLMM::MegaLMM_control(
            drop0_tol = 1e-10,
            scale_Y = FALSE,
            h2_divisions = 20,
            h2_step_size = NULL,
            burn = 0,
            K = 16,
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
          
          # Creating run ID:
          run_ID <- sprintf("p800/megalmm_states/%s_p800_dataset_h2y%s_comm%s_h2s%s_%d", CV, h2y, comm, h2s, sim)
          
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
          # (10 times for a total of 1000 burn-in samples):
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
          # 2000 iterations with a thinning rate of 2 gives 1000 posterior samples.
          #
          # So 1000 burn-in samples and 1000 posterior samples.
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
          mean_pred.test <- mean_pred.test[match(names(pred.target), names(mean_pred.test))]
          
          toc(log = TRUE)
          
          acc[run] <- cor(pred.target, mean_pred.test)
          
          # Deleting MegaLMM state files:
          unlink(run_ID, recursive = TRUE)
          
          # Increment run number:
          run <- run + 1
        }
        
        # Retrieve computational times:
        tictoc.logs <- tic.log(format = FALSE)
        tic.clearlog()
        comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
        
        # Collect results:
        results <- data.frame(acc = acc,
                              comptimes = comptimes)
        
        
        # Making correct CV label:
        if (CV == "CV1") {
          lab <- "a"
        } else if (CV == "CV2") {
          lab <- "b"
        }
        
        # Export results:
        write.csv(results, sprintf("p800/results/h2s%s/12%s_p800_results_MegaLMM_%s_h2y%s_comm%s_h2s%s_17to20.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
      }
    )
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    toc()
    combi <- combi + 1
  }
}
toc()


