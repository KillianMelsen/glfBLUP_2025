start <- Sys.time()

# Loading libraries:
library(rlist)
library(tictoc)
library(gfBLUP)

# Setting seed:
set.seed(1997)

# Total number of secondary features:
ps <- seq(100, 1400, 100)

# Loading kinship:
load(paste0(dirname(getwd()), "/genotypes/K_sim.RData")); rm(M)

# Setting up timing result storage:
# - First dimension is for the different values of p
# - Second dimension is used to store results from 10 timing runs after a single warmup run
# - Third dimension is used to store results
runs.warmup <- 1
runs.timing <- 5
steps <- c("Redundancy filtering",
           "Regularization",
           "Factor model",
           "Factor scores",
           "Subset selection",
           "gfBLUP genomic prediction",
           "Other")

results <- expand.grid(Step = steps, Run = 1:runs.timing, p = ps)
results$Duration <- numeric(nrow(results))

result.row.start <- 1
result.row.end <- length(steps)
p <- ps[1]
for (p in ps) {
  
  # Warmup =====================================================================
  run.warmup <- 1
  for (run.warmup in 1:runs.warmup) {
    
    # Loading simulated dataset:
    datalist <- list.load(file = sprintf("datasets/timing_p%d.RData", p))
  
    # Storing data and prediction target:
    d <- datalist$data.real
    pred.target <- datalist$pred.target
    
    tic("Full model") # Full model start
    
    
    
    ## 1. Make training data and store feature/focal trait names ===============
    d.train <- droplevels(na.omit(d))
    sec <- names(d[2:(ncol(d) - 1)])
    foc <- names(d)[ncol(d)]
    
    
    
    ## 2. Mock redundancy filtering the sec. features using training data ======
    tic("Redundancy filtering") # Redundancy filtering start
    temp <- gfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 0.95, verbose = F)
    toc(log = FALSE) # Redundancy filtering stop
    
    
    
    ## 3. Regularization =======================================================
    tic("Regularization") # Regularization start
    folds <- gfBLUP::createFolds(genos = unique(as.character(d.train$G)))
    tempG <- gfBLUP::regularizedCorrelation(data = d.train[c("G", sec)],
                                            folds = folds, what = "genetic", dopar = TRUE)
    
    tempE <- gfBLUP::regularizedCorrelation(data = d.train[c("G", sec)],
                                            folds = folds, what = "residual", dopar = TRUE)
    Rg.reg <- tempG$optCor
    toc(log = FALSE) # Regularization stop
    
    
    
    ## 4. Fitting factor model =================================================
    tic("Factor model") # Factor model start
    FM.fit <- gfBLUP::factorModel(data = d.train[c("G", sec)], cormat = Rg.reg, what = "genetic")
    toc(log = FALSE) # Factor model stop
    
    
    
    ## 5. Getting factor scores (also for the test set in CV2!) ================
    tic("Factor scores") # Factor scores start
    D <- sqrt(diag(tempG$Sg)) 
    L.cov <- diag(D) %*% FM.fit$loadings
    PSI.cov <- outer(D, D) * FM.fit$uniquenesses
    
    # CV2 Factor scores:
    # First recenter/rescale the training and test secondary data together:
    d[sec] <- sapply(d[sec], scale)
    CV2.F.scores <- gfBLUP::factorScores(data = d[c("G", sec)],
                                         loadings = L.cov,
                                         uniquenesses = PSI.cov,
                                         m = FM.fit$m,
                                         type = "genetic-thomson-repdiv",
                                         Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)
    
    CV2.d.final <- cbind(CV2.F.scores, d$Y)
    names(CV2.d.final)[ncol(CV2.d.final)] <- "Y"
    names(CV2.d.final)[1] <- "G"
    toc(log = FALSE) # Factor scores stop
    
    
    
    ## 6. Selecting the relevant factors =======================================
    tic("Subset selection") # Subset selection start
    # Subset selection can only be done using training data:
    selection <- gfBLUP::factorSelect(na.omit(CV2.d.final), procedure = "leaps", verbose = F)
    toc(log = FALSE) # Subset selection stop
    
    
    
    ## 7. Multi-trait genomic prediction =======================================
    tic("gfBLUP genomic prediction") # gfBLUP genomic prediction start
    CV2.temp <- gfBLUP::gfBLUP(data = CV2.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
    toc(log = FALSE) # gfBLUP genomic prediction stop
    
    toc(log = FALSE) # Full model stop
    
    cat(sprintf("Warmup run %d done, accuracy of %.3f\n\n",
                run.warmup, cor(pred.target, CV2.temp$preds[match(names(pred.target), names(CV2.temp$preds))])))
    
  }
  
  cat(sprintf("Proceeding to timing runs for p = %s\n\n", p))
  
  # Timing =====================================================================
  run.timing <- 1
  for (run.timing in 1:runs.timing) {
    
    # Loading simulated dataset:
    datalist <- list.load(file = sprintf("datasets/timing_p%d.RData", p))
    
    # Storing data and prediction target:
    d <- datalist$data.real
    pred.target <- datalist$pred.target
    
    tic("Full model") # Full model start
    
    ## 1. Make training data and store feature/focal trait names ===============
    d.train <- droplevels(na.omit(d))
    sec <- names(d[2:(ncol(d) - 1)])
    foc <- names(d)[ncol(d)]
    
    
    
    ## 2. Mock redundancy filtering the sec. features using training data ======
    tic("Redundancy filtering") # Redundancy filtering start
    temp <- gfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 0.95, verbose = F)
    toc(log = TRUE) # Redundancy filtering stop
    
    
    
    ## 3. Regularization =======================================================
    tic("Regularization") # Regularization start
    folds <- gfBLUP::createFolds(genos = unique(as.character(d.train$G)))
    tempG <- gfBLUP::regularizedCorrelation(data = d.train[c("G", sec)],
                                            folds = folds, what = "genetic", dopar = TRUE)
    
    tempE <- gfBLUP::regularizedCorrelation(data = d.train[c("G", sec)],
                                            folds = folds, what = "residual", dopar = TRUE)
    Rg.reg <- tempG$optCor
    toc(log = TRUE) # Regularization stop
    
    
    
    ## 4. Fitting factor model =================================================
    tic("Factor model") # Factor model start
    FM.fit <- gfBLUP::factorModel(data = d.train[c("G", sec)], cormat = Rg.reg, what = "genetic")
    toc(log = TRUE) # Factor model stop
    
    
    
    ## 5. Getting factor scores (also for the test set in CV2!) ================
    tic("Factor scores") # Factor scores start
    D <- sqrt(diag(tempG$Sg)) 
    L.cov <- diag(D) %*% FM.fit$loadings
    PSI.cov <- outer(D, D) * FM.fit$uniquenesses
    
    # CV2 Factor scores:
    # First recenter/rescale the training and test secondary data together:
    d[sec] <- sapply(d[sec], scale)
    CV2.F.scores <- gfBLUP::factorScores(data = d[c("G", sec)],
                                         loadings = L.cov,
                                         uniquenesses = PSI.cov,
                                         m = FM.fit$m,
                                         type = "genetic-thomson-repdiv",
                                         Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)
    
    CV2.d.final <- cbind(CV2.F.scores, d$Y)
    names(CV2.d.final)[ncol(CV2.d.final)] <- "Y"
    names(CV2.d.final)[1] <- "G"
    toc(log = TRUE) # Factor scores stop
    
    
    
    ## 6. Selecting the relevant factors =======================================
    tic("Subset selection") # Subset selection start
    # Subset selection can only be done using training data:
    selection <- gfBLUP::factorSelect(na.omit(CV2.d.final), procedure = "leaps", verbose = F)
    toc(log = TRUE) # Subset selection stop
    
    
    
    ## 7. Multi-trait genomic prediction =======================================
    tic("gfBLUP genomic prediction") # gfBLUP genomic prediction start
    CV2.temp <- gfBLUP::gfBLUP(data = CV2.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
    toc(log = TRUE) # gfBLUP genomic prediction stop
    
    toc(log = TRUE) # Full model stop
    
    cat(sprintf("Timing run %d done, accuracy of %.3f\n\n",
                run.timing, cor(pred.target, CV2.temp$preds[match(names(pred.target), names(CV2.temp$preds))])))
    
    times <- tic.log(format = FALSE)
    tic.clearlog()
    times <- unlist(lapply(times, function(x) x$toc - x$tic))
    names(times) <- steps
    times["Other"] <- times["Other"] - sum(times[1:6])
    
    results[result.row.start:result.row.end, "Durations"] <- as.numeric(times)
    result.row.start <- result.row.end + 1
    result.row.end <- result.row.start + 6
  }
  cat(sprintf("Timing runs for p = %s done, proceeding to warmup of next p\n\n", p))
}

# Saving results:
write.csv(results, "timing.csv")

end <- Sys.time()
end - start
