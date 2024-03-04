# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!
# Runtime: ~ s for 250 datasets.

# Loading libraries:
library(rlist)
library(tictoc)
library(doParallel)
library(gfBLUPold)
library(gfBLUP)

# Setting seed:
set.seed(1997)

# Setting working directory:
# wd <- "C:/Users/Killian/Desktop/gfblup-methodological-paper"
wd <- "~/gfblup_methodology"
setwd(wd)

# Loading kinship:
load("K_hyper.RData")

# Both CVs:
for (CV in c("CV1", "CV2")) {
  
  # Hyperspectral data:
  tic(sprintf("gfBLUP %s", CV))
  
  n.datasets <- 250
  n.cores <- 5
  work <- split(1:n.datasets, ceiling(seq_along(1:n.datasets) / ceiling(n.datasets / n.cores)))
  doParallel::registerDoParallel(cores = n.cores)
  
  invisible(
    par.results <- foreach::foreach(i = 1:length(work), .packages = c("rlist", "tictoc", "gfBLUPold", "gfBLUP"), .combine = "c") %dopar% {
      
      par.work <- work[[i]]
      
      # Setting up result storage:
      acc <- numeric(length(par.work))
      penG <- numeric(length(par.work))
      penE <- numeric(length(par.work))
      subset <- character(length(par.work))
      extra <- vector("list", length(par.work))
      
      # Running:
      for (run in 1:length(par.work)) {
        
        # Loading hyperspectral dataset:
        datalist <- list.load(file = sprintf("analyses_datasets/hyper/hyper_dataset_%d.RData", par.work[run]))
        
        # Storing data and prediction target:
        d <- datalist$data
        pred.target <- datalist$pred.target
        
        if (CV == "CV1") {
          d[which(is.na(d$Y)), 2:ncol(d)] <- NA
        }
        
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
        
        F.scores <- gfBLUP::factorScores(data = d.RF[c("G", sec.RF)],
                                         loadings = L.cov,
                                         uniquenesses = PSI.cov,
                                         m = FM.fit$m,
                                         type = "genetic-thomson",
                                         Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)
        
        d.final <- cbind(F.scores, d.RF$Y)
        names(d.final)[ncol(d.final)] <- "Y"
        names(d.final)[1] <- "G"
        
        #### 6. Selecting the relevant factors -----------------------------------------------------------------------------------------------
        selection <- gfBLUP::factorSelect(d.final, procedure = "leaps", verbose = FALSE)
        
        #### 7. Multi-trait genomic prediction -----------------------------------------------------------------------------------------------
        temp <- gfBLUP::gfBLUP(data = d.final, selection = selection, K = K, sepExp = FALSE, verbose = FALSE)
        
        toc(log = TRUE)
        ########################################################################
        
        acc[run] <- cor(pred.target$pred.target, temp$preds[match(pred.target$G, names(temp$preds))])
        penG[run] <- tempG$optPen
        penE[run] <- tempE$optPen
        subset[run] <- paste(selection, collapse = "-")
        
        extra[[run]] <- list(loadings = FM.fit$loadings,
                             uniquenesses = FM.fit$uniquenesses,
                             m = FM.fit$m,
                             m.selected = selection,
                             Sg = temp$Sg,
                             Se = temp$Se,
                             h2s = temp$h2s)
      }
      
      # Retrieve computational times:
      tictoc.logs <- tic.log(format = FALSE)
      tic.clearlog()
      comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
      
      # Collect results:
      worker.result <- list(list(result = data.frame(acc = acc,
                                                     penG = penG,
                                                     penE = penE,
                                                     subset = subset,
                                                     comptimes = comptimes),
                                 extra = extra))
      
      names(worker.result) <- sprintf("worker_%d", i)
      return(worker.result)
      
    })
  doParallel::stopImplicitCluster()
  toc()
  
  # Restructuring parallel results for saving:
  acc <- numeric()
  penG <- numeric()
  penE <- numeric()
  subset <- character()
  comptimes <- numeric()
  extra <- vector("list")
  
  for (j in 1:length(work)) {
    
    acc <- c(acc, par.results[[sprintf("worker_%d", j)]]$result$acc)
    penG <- c(penG, par.results[[sprintf("worker_%d", j)]]$result$penG)
    penE <- c(penE, par.results[[sprintf("worker_%d", j)]]$result$penE)
    subset <- c(subset, par.results[[sprintf("worker_%d", j)]]$result$subset)
    comptimes <- c(comptimes, par.results[[sprintf("worker_%d", j)]]$result$comptimes)
    extra <- c(extra, par.results[[sprintf("worker_%d", j)]]$extra)
    
  }
  
  results <- data.frame(acc = acc,
                        penG = penG,
                        penE = penE,
                        subset = subset,
                        comptimes = comptimes)
  
  # Making correct CV label:
  if (CV == "CV1") {
    lab <- "a"
  } else if (CV == "CV2") {
    lab <- "b"
  }
  
  # Export results:
  write.csv(results, sprintf("analyses_results/hyper/3%s_hyper_results_gfblup_%s.csv", lab, CV))
  
  list.save(extra, sprintf("analyses_results/hyper/3%s_hyper_extra_results_gfblup_%s.RData", lab, CV))
  
}




