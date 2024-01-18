# Overall ======================================================================

# Loading libraries:
library(rlist)
library(tictoc)
library(MegaLMM)

# Numbers of MegaLMM factors:
Ms <- c(5, 10, 20, 50)

# Number of burn-in iterations (swithcing factors 10 times):
burnin <- 5000

# Number of sample iterations (thinning rate of 2)
posterior <- 50000
for (M in Ms) {
  
  # Setting seed:
  set.seed(1997)
  
  # Setting working directory:
  setwd("C:/Users/killi/Desktop/LINUX/MEGALMM")
  
  # Loading kinship:
  load("K_hyper.RData")
  
  # Loading simulated dataset:
  datalist <- list.load(file = "hyper_dataset_1.RData")
  
  # Storing data and prediction target:
  d <- datalist$data
  
  # +++++++++++++++++++++++++++++++++++++++++
  # Looking only at a single date:
  select <- c(1, seq(7, 617, 10), ncol(d))
  d <- d[, select]
  # +++++++++++++++++++++++++++++++++++++++++
  
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
  
  tic(sprintf("M%s", M))
  
  # MegaLMM config:
  run_parameters <- MegaLMM::MegaLMM_control(
    drop0_tol = 1e-10,
    scale_Y = FALSE,
    h2_divisions = 20,
    h2_step_size = NULL,
    burn = 0,
    K = M,
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
  run_ID <- sprintf("states/CV2_hyper_dataset_1_M%s", M)
  
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
  
  toc(log = TRUE)
  
  # Saving the arrays containing all posterior predictions and loadings
  save(Lambda_samples, file = sprintf("posterior_arrays/LAMBDA_M%s.RData", M))
  save(pred_samples, file = sprintf("posterior_arrays/PREDS_M%s.RData", M))
  
  # Traceplots:
  traceplot_array(Lambda_samples, facet_dim = 3, n_per_facet = M, name = sprintf("MegaLMM_traceplots/LAMBDA_M%s.pdf", M))
  traceplot_array(pred_samples, facet_dim = 2, name = sprintf("MegaLMM_traceplots/PREDS_M%s.pdf", M))
  
  # Calculating posterior means:
  mean_Lambda <- MegaLMM::get_posterior_mean(Lambda_samples)
  mean_pred <- MegaLMM::get_posterior_mean(pred_samples)
  
  # Plotting +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  loadings <- as.data.frame(t(mean_Lambda))
  names(loadings) <- paste0("F", 1:nrow(mean_Lambda))
  loadings.wl <- loadings[-1,]
  loadings.y <- loadings[1,]
  
  loadings.wl$wavelength <- as.numeric(substr(rownames(loadings.wl), 3, 5))
  
  loadings.wl.long <- tidyr::pivot_longer(loadings.wl, names_to = "Factor", cols = 1:(ncol(loadings.wl) - 1))
  
  save(loadings.y, file = sprintf("mean_posterior_loadings_Y/LAMBDA_Y_M%s.RData", M))
  
  # Find useful factors (i.e., those with absolute Y loadings higher than 0.05):
  useful <- names(loadings.y[1,])[which(abs(loadings.y[1,]) > 0.05)]
  loadings.wl.long$Useful <- ifelse(loadings.wl.long$Factor %in% useful, "Yes", "No")
  
  pal <- RColorBrewer::brewer.pal(12, "Paired")
  
  colors <- character(nrow(mean_Lambda))
  i <- 1
  j <- 1
  for (f in paste0("F", 1:nrow(mean_Lambda))) {
    if (f %in% useful) {
      colors[j] <- pal[i]
      i <- i + 1
    } else {
      colors[j] <- "black"
    }
    j <- j + 1
  }
  
  names(colors) <- paste0("F", 1:nrow(mean_Lambda))
  
  library(ggplot2)
  
  loadings.wl.long$Factor <- factor(loadings.wl.long$Factor, levels = paste0("F", 1:nrow(mean_Lambda)))
  
  ggplot(data = loadings.wl.long, mapping = aes(x = wavelength, y = value, color = Factor)) +
    geom_point(mapping = aes(shape = Useful)) +
    theme_dark() +
    ylim(c(-1, 1)) +
    scale_color_manual(values = colors)
  
  ggsave(filename = sprintf("mean_posterior_loadings_WL/LAMBDA_WL_M%s.jpg", M), height = 20, width = 40, units = "cm")
  
  rm(list = setdiff(ls(), c("burnin", "Ms", "posterior")))
}

tic.log()
