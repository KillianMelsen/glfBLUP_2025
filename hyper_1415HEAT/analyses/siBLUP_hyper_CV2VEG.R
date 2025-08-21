# This script produces all intermediate results for the siBLUP CV2VEG analyses on
# the hyperspectral HEAT dataset.
#
# The script is structured as follows:
# - There is an outer loop that would loop over the different CVs (CV1/CV2), but
#   this script only does CV2. So manually set `CV <- "CV2"` and `lab <- "b"`.
# - Then there is a parallel loop dividing the 250 replications over 18 parallel
#   workers.
# - The replications produced by each worker are given by the `par.work`
#   variable on line 65.
# - A seed is set inside each worker and the replicates of the worker are produced
#   inside the inner loop starting on line 74.
#
# Suppose you want to reproduce the result of replication 215. Looking at the `work`
# variable, replication 215 is handled by the 16th worker. So we can manually set
# `i <- 16` instead of running the parallel for loop. We also note that replication
# 215 is the 5th replication that is produced by the 16th worker, so in the inner
# loop we do not go from `1:length(par.work)`, but `1:5`. Then running the loop
# will produce the CV2VEG accuracy of replication 215 as the 5th value in `acc`.
# This can be compared to the 215th accuracy stored in the appropriate intermediate
# results csv file:
#
# all.equal(read.csv(sprintf("hyper_1415HEAT/results/%s/5%s_hyper_results_siblup_%sVEG.csv", prep, lab, CV))$acc[215],
#           acc[5])
#
# [1] TRUE

prep <- "nosplines"

# Loading libraries:
library(rlist)
library(tictoc)
library(doParallel)
library(gfBLUPold)
source("helper_functions.R")
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
for (CV in c("CV2")) {
  
  # Hyperspectral data:
  tic(sprintf("siBLUP %s", CV))
  
  n.datasets <- 250
  n.cores <- 18
  work <- split(1:n.datasets, ceiling(seq_along(1:n.datasets) / ceiling(n.datasets / n.cores)))
  cl <- parallel::makeCluster(n.cores, outfile = sprintf("logs/siBLUP_hyper_%sVEG.txt", CV))
  doParallel::registerDoParallel(cl)
  
  invisible(
    par.results <- foreach::foreach(i = 1:length(work), .packages = c("rlist", "tictoc", "gfBLUPold"), .combine = "c") %dopar% {
      
      par.work <- work[[i]]
      set.seed(1997)
      
      # Setting up result storage:
      acc <- numeric(length(par.work))
      pen <- numeric(length(par.work))
      extra <- vector("list", length(par.work))
      
      # Running:
      for (run in 1:length(par.work)) {
        
        # Loading hyperspectral dataset:
        datalist <- list.load(file = sprintf("hyper_1415HEAT/datasets/%s/hyper_dataset_%d.RData", prep, par.work[run]))
        
        # Storing data and prediction target:
        # 9 feb is last day of VEG, 25 feb is heading, 10 march is start of grain filling:
        dates <- c("150414")
        d <- datalist$data
        select <- which(substr(names(d), 7, 12) %in% dates)
        d <- d[c(1, select, ncol(d))]
        pred.target <- datalist$pred.target
        
        # Subsetting K (only really happens for the first dataset...):
        K <- K[unique(d$G), unique(d$G)]
        
        ### Model ##############################################################
        tic(run)
        RESULT <- gfBLUPold::siBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, t.RF = 0.95, verbose = FALSE, sepExp = FALSE)
        toc(log = TRUE)
        ########################################################################
        
        RESULT$preds <- RESULT$preds[match(pred.target$G, names(RESULT$preds))]
        
        #### Runcie & Cheng 2019 correction ----------------------------------
        temp <- estimate_gcor(data = data.frame(ID = pred.target$G,
                                                obs = pred.target$pred.target,
                                                pred = RESULT$preds),
                              Knn = K[pred.target$G, pred.target$G],
                              method = "MCMCglmm",
                              normalize = T)
        acc[run] <- temp["g_cor"]
        pen[run] <- RESULT$regPen
        
        extra[[run]] <- list(gamma = RESULT$gamma,
                             Vg = RESULT$Vg,
                             Ve = RESULT$Ve,
                             H2s = RESULT$H2s,
                             penSI = RESULT$penSI,
                             AM.direct = RESULT$AM.direct,
                             AM.indirect = RESULT$AM.indirect)
      }
      
      
      
      
      # Retrieve computational times:
      tictoc.logs <- tic.log(format = FALSE)
      tic.clearlog()
      comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
      
      # Collect results:
      worker.result <- list(list(result = data.frame(acc = acc,
                                                     pen = pen,
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
  pen <- numeric()
  comptimes <- numeric()
  extra <- vector("list")
  
  for (j in 1:length(work)) {
    
    acc <- c(acc, par.results[[sprintf("worker_%d", j)]]$result$acc)
    pen <- c(pen, par.results[[sprintf("worker_%d", j)]]$result$pen)
    comptimes <- c(comptimes, par.results[[sprintf("worker_%d", j)]]$result$comptimes)
    extra <- c(extra, par.results[[sprintf("worker_%d", j)]]$extra)
    
  }
  
  results <- data.frame(acc = acc,
                        pen = pen,
                        comptimes = comptimes)
  
  # Making correct CV label:
  if (CV == "CV1") {
    lab <- "a"
  } else if (CV == "CV2") {
    lab <- "b"
  }
  
  # Export results:
  write.csv(results, sprintf("hyper_1415HEAT/results/%s/5%s_hyper_results_siblup_%sVEG.csv", prep, lab, CV))
  
  list.save(extra, sprintf("hyper_1415HEAT/results/%s/5%s_hyper_extra_results_siblup_%sVEG.RData", prep, lab, CV))
  
}




