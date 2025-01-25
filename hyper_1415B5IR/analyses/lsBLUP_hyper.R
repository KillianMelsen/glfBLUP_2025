# !!! IMPORTANT: SET MKL_NUMTHREADS=1 if using Intel MKL
# Runtime: ~ 85s for 36 datasets using 18 workers.
prep <- "splines"
# Loading libraries:
library(rlist)
library(tictoc)
library(doParallel)
library(gfBLUPold)
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

# Both CVs:
for (CV in c("CV1", "CV2")) {
  
  # Hyperspectral data:
  tic(sprintf("lsBLUP %s", CV))
  
  n.datasets <- 250
  n.cores <- parallel::detectCores() - 2
  work <- split(1:n.datasets, ceiling(seq_along(1:n.datasets) / ceiling(n.datasets / n.cores)))
  cl <- parallel::makeCluster(n.cores, outfile = sprintf("logs/lsBLUP_hyper_%s.txt", CV))
  doParallel::registerDoParallel(cl)
  
  invisible(
    par.results <- foreach::foreach(i = 1:length(work), .packages = c("rlist", "tictoc", "gfBLUPold"), .combine = "c") %dopar% {
      
      par.work <- work[[i]]
      set.seed(1997)
      
      # Setting up result storage:
      acc <- numeric(length(par.work))
      extra <- vector("list", length(par.work))
      
      
      # Running:
      for (run in 1:length(par.work)) {
        
        # Loading hyperspectral dataset:
        datalist <- list.load(file = sprintf("hyper_1415B5IR/datasets/%s/hyper_dataset_%d.RData", prep, par.work[run]))
        
        # Storing data and prediction target:
        d <- datalist$data
        pred.target <- datalist$pred.target
        
        # Subsetting K (only really happens for the first dataset...):
        K <- K[unique(d$G), unique(d$G)]
        
        ### Model ##############################################################
        tic(run)
        RESULT <- gfBLUPold::lsBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, sepExp = FALSE, t.RF = 0.95, verbose = FALSE)
        toc(log = TRUE)
        ########################################################################
        
        RESULT$preds <- RESULT$preds[match(pred.target$G, names(RESULT$preds))]
        
        if (CV == "CV1") {
          acc[run] <- cor(RESULT$preds, pred.target$pred.target)
        } else if (CV == "CV2") {
          #### Runcie & Cheng 2019 correction ----------------------------------
          temp <- estimate_gcor(data = data.frame(ID = pred.target$G,
                                                  obs = pred.target$pred.target,
                                                  pred = RESULT$preds),
                                Knn = K[pred.target$G, pred.target$G],
                                method = "MCMCglmm",
                                normalize = T)
          acc[run] <- temp["g_cor"]
        }
        
        if (length(RESULT) == 1) {
          extra[[run]] <- "All coefficients were 0, used univariate model..."
        } else {
          extra[[run]] <- list(coefs = RESULT$coefs,
                               Vg = RESULT$Vg,
                               Ve = RESULT$Ve,
                               H2s = RESULT$H2s,
                               LSP = RESULT$LSP,
                               AM.direct = RESULT$AM.direct,
                               AM.indirect = RESULT$AM.indirect)
        }
      }
      
      # Retrieve computational times:
      tictoc.logs <- tic.log(format = FALSE)
      tic.clearlog()
      comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
      
      # Collect results:
      worker.result <- list(list(result = data.frame(acc = acc,
                                                     comptimes = comptimes),
                                 extra = extra))
      names(worker.result) <- sprintf("worker_%d", i)
      return(worker.result)
      
    })
  doParallel::stopImplicitCluster()
  parallel::stopCluster(cl)
  toc()
  
  # Restructuring parallel results for saving:
  acc <- numeric()
  comptimes <- numeric()
  extra <- vector("list")
  
  for (j in 1:length(work)) {
    
    acc <- c(acc, par.results[[sprintf("worker_%d", j)]]$result$acc)
    comptimes <- c(comptimes, par.results[[sprintf("worker_%d", j)]]$result$comptimes)
    extra <- c(extra, par.results[[sprintf("worker_%d", j)]]$extra)
    
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
  write.csv(results, sprintf("hyper_1415B5IR/results/%s/11%s_hyper_results_lsblup_%s.csv", prep, lab, CV))
  
  list.save(extra, sprintf("hyper_1415B5IR/results/%s/11%s_hyper_extra_results_lsblup_%s.RData", prep, lab, CV))
  
}




