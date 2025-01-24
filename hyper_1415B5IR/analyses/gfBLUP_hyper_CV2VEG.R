# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!
# Runtime: ~ s for 250 datasets.
prep <- "VEGsplines"
# Loading libraries:
library(rlist)
library(tictoc)
library(doParallel)
library(gfBLUPold)
library(gfBLUP)
source("helper_functions/Estimate_gcor_prediction.R")
library(MCMCglmm)
library(coda)
library(ape)

# Setting seed:
set.seed(1997)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading kinship:
load("genotypes/K_hyper.RData")

# Hyperspectral data:
tic("gfBLUP")

n.datasets <- 250
n.cores <- 5
work <- split(1:n.datasets, ceiling(seq_along(1:n.datasets) / ceiling(n.datasets / n.cores)))
cl <- parallel::makeCluster(n.cores, outfile = "logs/gfBLUP_hyper_CV2VEG.txt")
doParallel::registerDoParallel(cl)

invisible(
  par.results <- foreach::foreach(i = 1:length(work), .packages = c("rlist", "tictoc", "gfBLUPold", "gfBLUP"), .combine = "c") %dopar% {
    
    par.work <- work[[i]]
    set.seed(1997)
    
    # Setting up result storage:
    CV1.acc <- CV2.acc <- CV2.RC.acc <- numeric(length(par.work))
    penG <- numeric(length(par.work))
    penE <- numeric(length(par.work))
    subset <- character(length(par.work))
    extra <- vector("list", length(par.work))
    
    # Running:
    for (run in 1:length(par.work)) {
      
      # Loading hyperspectral dataset:
      datalist <- list.load(file = sprintf("hyper/datasets/%s/hyper_dataset_%d.RData", prep, par.work[run]))
      
      # Storing data and prediction target:
      # 9 feb is last day of VEG, 25 feb is heading, 10 march is start of grain filling:
      dates <- c("150110", "150119", "150204", "150209")
      d <- datalist$data
      select <- which(substr(names(d), 7, 12) %in% dates)
      d <- d[c(1, select, ncol(d))]
      pred.target <- datalist$pred.target
      
      ### Model ##############################################################
      tic(run)
      
      ### 1. Make training data and store feature/focal trait names ------------------------------------------------------------------------
      d.train <- droplevels(na.omit(d))
      sec <- names(d[2:(ncol(d) - 1)])
      foc <- names(d)[ncol(d)]
      
      ### 2. Redundancy filter the secondary features using training data only -------------------------------------------------------------
      temp <- gfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 0.95, verbose = FALSE)
      d.train.RF <- cbind(temp$data.RF, d.train[foc])
      d.RF <- d[names(d.train.RF)]
      sec.RF <- names(d.RF[2:(ncol(d.RF) - 1)])
      
      ### 3. Regularization ----------------------------------------------------------------------------------------------------------------
      folds <- gfBLUP::createFolds(genos = unique(as.character(d.train.RF$G)))
      tempG <- gfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "genetic", dopar = FALSE, verbose = FALSE)
      tempE <- gfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "residual", dopar = FALSE, verbose = FALSE)
      Rg.RF.reg <- tempG$optCor
      
      ### 4. Fitting factor model ----------------------------------------------------------------------------------------------------------
      # data is only used to determine the sample size for the MP-bound. what specifies that it's a genetic correlation matrix, so the
      # number of training genotypes should be used, and not the number of training individuals (= genotypes * replicates).
      FM.fit <- gfBLUP::factorModel(data = d.train.RF[c("G", sec.RF)], cormat = Rg.RF.reg, what = "genetic", verbose = FALSE)
      
      #### 5. Getting factor scores (also for the test set in CV2!) ------------------------------------------------------------------------
      # Loadings and uniquenesses were estimated on the correlation scale, but should be on the covariance scale for genetic-thomson scores:
      D <- sqrt(diag(tempG$Sg)) # Getting standard deviations
      L.cov <- diag(D) %*% FM.fit$loadings
      PSI.cov <- outer(D, D) * FM.fit$uniquenesses
      
      # # CV1 Factor scores:
      # CV1.d.RF <- d.RF
      # CV1.d.RF[which(is.na(CV1.d.RF$Y)), 2:ncol(CV1.d.RF)] <- NA
      # CV1.F.scores <- gfBLUP::factorScores(data = CV1.d.RF[c("G", sec.RF)],
      #                                      loadings = L.cov,
      #                                      uniquenesses = PSI.cov,
      #                                      m = FM.fit$m,
      #                                      type = "genetic-thomson-repdiv",
      #                                      Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)
      # 
      # CV1.d.final <- cbind(CV1.F.scores, CV1.d.RF$Y)
      # names(CV1.d.final)[ncol(CV1.d.final)] <- "Y"
      # names(CV1.d.final)[1] <- "G"
      
      # CV2 Factor scores:
      # First recenter/rescale the training and test secondary data together:
      d.RF[sec.RF] <- sapply(d.RF[sec.RF], scale)
      CV2.F.scores <- gfBLUP::factorScores(data = d.RF[c("G", sec.RF)],
                                           loadings = L.cov,
                                           uniquenesses = PSI.cov,
                                           m = FM.fit$m,
                                           type = "genetic-thomson-repdiv",
                                           Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)
      
      CV2.d.final <- cbind(CV2.F.scores, d.RF$Y)
      names(CV2.d.final)[ncol(CV2.d.final)] <- "Y"
      names(CV2.d.final)[1] <- "G"
      
      #### 6. Selecting the relevant factors -----------------------------------------------------------------------------------------------
      selection <- gfBLUP::factorSelect(CV2.d.final, procedure = "leaps", verbose = FALSE)
      
      #### 7. Multi-trait genomic prediction -----------------------------------------------------------------------------------------------
      # CV1.temp <- gfBLUP::gfBLUP(data = CV1.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
      CV2.temp <- gfBLUP::gfBLUP(data = CV2.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
      toc(log = TRUE)
      ########################################################################
      
      # CV1.acc[run] <- cor(pred.target$pred.target, CV1.temp$preds[match(pred.target$G, names(CV1.temp$preds))])
      CV2.acc[run] <- cor(pred.target$pred.target, CV2.temp$preds[match(pred.target$G, names(CV2.temp$preds))])
      #### Runcie & Cheng 2019 correction --------------------------------------
      temp <- estimate_gcor(data = data.frame(ID = pred.target$G,
                                              obs = pred.target$pred.target,
                                              pred = CV2.temp$preds[match(pred.target$G, names(CV2.temp$preds))]),
                            Knn = K[pred.target$G, pred.target$G],
                            method = "MCMCglmm",
                            normalize = T)
      CV2.RC.acc[run] <- temp["g_cor"]
      
      penG[run] <- tempG$optPen
      penE[run] <- tempE$optPen
      subset[run] <- paste(selection, collapse = "-")
      
      extra[[run]] <- list(loadings = FM.fit$loadings,
                           uniquenesses = FM.fit$uniquenesses,
                           m = FM.fit$m,
                           m.selected = selection,
                           Sg = CV2.temp$Sg,
                           Se = CV2.temp$Se,
                           h2s = CV2.temp$h2s)
    }
    
    # Retrieve computational times:
    tictoc.logs <- tic.log(format = FALSE)
    tic.clearlog()
    comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
    
    # Collect results:
    worker.result <- list(list(result = data.frame(#CV1.acc = CV1.acc,
                                                   CV2.acc = CV2.acc,
                                                   CV2.RC.acc = CV2.RC.acc,
                                                   penG = penG,
                                                   penE = penE,
                                                   subset = subset,
                                                   comptimes = comptimes),
                               extra = extra))
    
    names(worker.result) <- sprintf("worker_%d", i)
    return(worker.result)
    
  })
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
toc()

# Restructuring parallel results for saving:
CV1.acc <- CV2.acc <- CV2.RC.acc <- numeric()
penG <- numeric()
penE <- numeric()
subset <- character()
comptimes <- numeric()
extra <- vector("list")

for (j in 1:length(work)) {
  
  # CV1.acc <- c(CV1.acc, par.results[[sprintf("worker_%d", j)]]$result$CV1.acc)
  CV2.acc <- c(CV2.acc, par.results[[sprintf("worker_%d", j)]]$result$CV2.acc)
  CV2.RC.acc <- c(CV2.RC.acc, par.results[[sprintf("worker_%d", j)]]$result$CV2.RC.acc)
  penG <- c(penG, par.results[[sprintf("worker_%d", j)]]$result$penG)
  penE <- c(penE, par.results[[sprintf("worker_%d", j)]]$result$penE)
  subset <- c(subset, par.results[[sprintf("worker_%d", j)]]$result$subset)
  comptimes <- c(comptimes, par.results[[sprintf("worker_%d", j)]]$result$comptimes)
  extra <- c(extra, par.results[[sprintf("worker_%d", j)]]$extra)
  
}

# CV1.results <- data.frame(acc = CV1.acc,
#                           penG = penG,
#                           penE = penE,
#                           subset = subset,
#                           comptimes = comptimes)

CV2.results <- data.frame(acc = CV2.acc,
                          acc.RC = CV2.RC.acc,
                          penG = penG,
                          penE = penE,
                          subset = subset,
                          comptimes = comptimes)

# Export results:
# write.csv(CV1.results, "hyper/results/3a_hyper_results_gfblup_CV1_RF.csv")
write.csv(CV2.results, sprintf("hyper/results/%s/3b_hyper_results_gfblup_CV2VEG.csv", prep))

list.save(extra, sprintf("hyper/results/%s/3_hyper_extra_results_gfblup_CV2VEG.RData", prep))




