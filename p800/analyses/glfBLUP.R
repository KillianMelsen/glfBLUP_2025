# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!
# Should take about 4.5 hours to run

# Loading libraries:
library(rlist)
library(tictoc)
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

combis <- length(h2.foc)
combi <- 1

par.combis <- data.frame(comm = rep(comms, each = length(h2.sec)),
                         h2s = rep(h2.sec, length(comms)))

# All simulated data for both CV1 and CV2:
tic("glfBLUP")
for (h2y in h2.foc) {
  
  tic(sprintf("glfBLUP CV1/CV2 p800 combi %d / %d, (h2y = %s)", combi, combis, h2y))
  
  cl <- parallel::makeCluster(9, outfile = sprintf("logs/glfBLUP_sim_p800_h2y%s.txt", h2y))
  registerDoParallel(cl)
  
  invisible(
  foreach::foreach(i = 1:nrow(par.combis), .packages = c("rlist", "tictoc")) %dopar% {
    
    comm <- par.combis[i, "comm"]
    h2s <- par.combis[i, "h2s"]
    
    # Number of simulated datasets to load:
    first <- 1
    last <- 100
    n.sim <- length(first:last)
    
    # Setting up result storage:
    CV1.acc <- numeric(n.sim)
    CV2.acc <- numeric(n.sim)
    penG <- rep(NA, n.sim)
    penE <- rep(NA, n.sim)
    subset <- rep(NA, n.sim)
    extra <- vector("list", n.sim)
    
    # Running (SET SEED IN EACH PARALLEL WORKER!):
    set.seed(1997)
    for (sim in first:last) {
      
      cat(sprintf("Running glfBLUP CV1/CV2 on p800 dataset %d / %d, (h2s = %s, comm = %s, h2y = %s), combi %d / %d\n",
                  sim, n.sim, h2s, comm, h2y, combi, combis))
      
      # Loading simulated dataset:
      datalist <- list.load(file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
      
      # Storing data and prediction target:
      d <- datalist$data.real
      pred.target <- datalist$pred.target
      
      ### Model ############################################################
      tic(sim)
      ### 1. Make training data and store feature/focal trait names ------------------------------------------------------------------------
      d.train <- droplevels(na.omit(d))
      sec <- names(d[2:(ncol(d) - 1)])
      foc <- names(d)[ncol(d)]
      
      ### 2. Redundancy filter the secondary features using training data only -------------------------------------------------------------
      temp <- glfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 1, verbose = F)
      d.train.RF <- cbind(temp$data.RF, d.train[foc])
      d.RF <- d[names(d.train.RF)]
      sec.RF <- names(d.RF[2:(ncol(d.RF) - 1)])
      
      ### 3. Regularization ----------------------------------------------------------------------------------------------------------------
      folds <- glfBLUP::createFolds(genos = unique(as.character(d.train.RF$G)))
      tempG <- glfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "genetic", dopar = FALSE)
      tempE <- glfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "residual", dopar = FALSE)
      Rg.RF.reg <- tempG$optCor
      
      ### 4. Fitting factor model ----------------------------------------------------------------------------------------------------------
      # data is only used to determine the sample size for the MP-bound. what specifies that it's a genetic correlation matrix, so the
      # number of training genotypes should be used, and not the number of training individuals (= genotypes * replicates).
      FM.fit <- glfBLUP::factorModel(data = d.train.RF[c("G", sec.RF)], cormat = Rg.RF.reg, what = "genetic")
      
      #### 5. Getting factor scores (also for the test set in CV2!) ------------------------------------------------------------------------
      # Loadings and uniquenesses were estimated on the correlation scale, but should be on the covariance scale for genetic-thomson scores:
      D <- sqrt(diag(tempG$Sg)) # Getting standard deviations
      L.cov <- diag(D) %*% FM.fit$loadings
      PSI.cov <- outer(D, D) * FM.fit$uniquenesses
      
      # CV1 Factor scores:
      CV1.d.RF <- d.RF
      CV1.d.RF[which(is.na(CV1.d.RF$Y)), 2:ncol(CV1.d.RF)] <- NA
      CV1.F.scores <- glfBLUP::factorScores(data = CV1.d.RF[c("G", sec.RF)],
                                           loadings = L.cov,
                                           uniquenesses = PSI.cov,
                                           m = FM.fit$m,
                                           type = "genetic-thomson-repdiv",
                                           Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)
      
      CV1.d.final <- cbind(CV1.F.scores, CV1.d.RF$Y)
      names(CV1.d.final)[ncol(CV1.d.final)] <- "Y"
      names(CV1.d.final)[1] <- "G"
      
      # CV2 Factor scores:
      # First recenter/rescale the training and test secondary data together:
      d.RF[sec.RF] <- sapply(d.RF[sec.RF], scale)
      CV2.F.scores <- glfBLUP::factorScores(data = d.RF[c("G", sec.RF)],
                                           loadings = L.cov,
                                           uniquenesses = PSI.cov,
                                           m = FM.fit$m,
                                           type = "genetic-thomson-repdiv",
                                           Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)
      
      CV2.d.final <- cbind(CV2.F.scores, d.RF$Y)
      names(CV2.d.final)[ncol(CV2.d.final)] <- "Y"
      names(CV2.d.final)[1] <- "G"
      
      #### 6. Selecting the relevant factors -----------------------------------------------------------------------------------------------
      selection <- glfBLUP::factorSelect(CV1.d.final, procedure = "leaps", verbose = F)
      
      #### 7. Multi-trait genomic prediction -----------------------------------------------------------------------------------------------
      CV1.temp <- glfBLUP::glfBLUP(data = CV1.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
      CV2.temp <- glfBLUP::glfBLUP(data = CV2.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
      toc(log = TRUE)
      ######################################################################
      
      CV1.acc[sim] <- cor(pred.target, CV1.temp$preds[match(names(pred.target), names(CV1.temp$preds))])
      CV2.acc[sim] <- cor(pred.target, CV2.temp$preds[match(names(pred.target), names(CV2.temp$preds))])
      penG[sim] <- tempG$optPen
      penE[sim] <- tempE$optPen
      subset[sim] <- paste(selection, collapse = "-")
      extra[[sim]] <- list(loadings = FM.fit$loadings,
                           uniquenesses = FM.fit$uniquenesses,
                           m = FM.fit$m,
                           m.selected = selection,
                           Sg = CV1.temp$Sg,
                           Se = CV1.temp$Se,
                           h2s = CV1.temp$h2s)
      
    }
    
    # Retrieve computational times:
    tictoc.logs <- tic.log(format = FALSE)
    tic.clearlog()
    comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
    
    # Collect results:
    CV1.results <- data.frame(acc = CV1.acc,
                              comptimes = comptimes,
                              penG = penG,
                              penE = penE,
                              subset = subset)
    
    CV2.results <- data.frame(acc = CV2.acc,
                              comptimes = comptimes,
                              penG = penG,
                              penE = penE,
                              subset = subset)
    
    # Export results:
    write.csv(CV1.results, sprintf("p800/results/h2s%s/3a_p800_results_glfblup_CV1_h2y%s_comm%s_h2s%s.csv", h2s, h2y, comm, h2s))
    write.csv(CV2.results, sprintf("p800/results/h2s%s/3b_p800_results_glfblup_CV2_h2y%s_comm%s_h2s%s.csv", h2s, h2y, comm, h2s))
    
    
    list.save(extra, file = sprintf("p800/results/h2s%s/3a_p800_extra_results_glfblup_CV1_h2y%s_comm%s_h2s%s.RData", h2s, h2y, comm, h2s))
    list.save(extra, file = sprintf("p800/results/h2s%s/3b_p800_extra_results_glfblup_CV2_h2y%s_comm%s_h2s%s.RData", h2s, h2y, comm, h2s))
    
  })
  doParallel::stopImplicitCluster()
  parallel::stopCluster(cl)
  toc()
  combi <- combi + 1
}
toc()

