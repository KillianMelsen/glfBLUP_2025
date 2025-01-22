# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!
# Runtime: ~1400 s for 250 datasets.
prep <- "nosplines"
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
for (CV in c("CV2")) {
  
  # Hyperspectral data:
  tic(sprintf("siBLUP %s", CV))
  
  n.datasets <- 250
  n.cores <- parallel::detectCores() - 2
  work <- split(1:n.datasets, ceiling(seq_along(1:n.datasets) / ceiling(n.datasets / n.cores)))
  cl <- parallel::makeCluster(n.cores, outfile = sprintf("logs/siBLUP_hyper_%sVEG.txt", CV))
  doParallel::registerDoParallel(cl)
  
  invisible(
    par.results <- foreach::foreach(i = 1:length(work), .packages = c("rlist", "tictoc", "gfBLUPold"), .combine = "c") %dopar% {
      
      par.work <- work[[i]]
      set.seed(1997)
      
      # Setting up result storage:
      acc <- numeric(length(par.work))
      CV2.RC.acc <- numeric(length(par.work))
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
        CV2.RC.acc[run] <- temp["g_cor"]
        acc[run] <- cor(RESULT$preds, pred.target$pred.target)
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
                                                     CV2.RC.acc = CV2.RC.acc,
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
  CV2.RC.acc <- numeric()
  pen <- numeric()
  comptimes <- numeric()
  extra <- vector("list")
  
  for (j in 1:length(work)) {
    
    acc <- c(acc, par.results[[sprintf("worker_%d", j)]]$result$acc)
    CV2.RC.acc <- c(CV2.RC.acc, par.results[[sprintf("worker_%d", j)]]$result$CV2.RC.acc)
    pen <- c(pen, par.results[[sprintf("worker_%d", j)]]$result$pen)
    comptimes <- c(comptimes, par.results[[sprintf("worker_%d", j)]]$result$comptimes)
    extra <- c(extra, par.results[[sprintf("worker_%d", j)]]$extra)
    
  }
  
  results <- data.frame(acc = acc,
                        acc.RC = CV2.RC.acc,
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




