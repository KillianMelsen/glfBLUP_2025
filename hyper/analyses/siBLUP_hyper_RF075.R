# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!
# Runtime: ~ s for 250 datasets.

# Loading libraries:
library(rlist)
library(tictoc)
library(doParallel)
library(gfBLUPold)

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
  tic(sprintf("siBLUP RF075 %s", CV))
  
  n.datasets <- 250
  n.cores <- 16
  work <- split(1:n.datasets, ceiling(seq_along(1:n.datasets) / ceiling(n.datasets / n.cores)))
  doParallel::registerDoParallel(cores = n.cores)
  
  invisible(
    par.results <- foreach::foreach(i = 1:length(work), .packages = c("rlist", "tictoc", "gfBLUPold"), .combine = "c") %dopar% {
      
      par.work <- work[[i]]
      
      # Setting up result storage:
      acc <- numeric(length(par.work))
      pen <- numeric(length(par.work))
      extra <- vector("list", length(par.work))
      
      # Running:
      for (run in 1:length(par.work)) {
        
        # Loading hyperspectral dataset:
        datalist <- list.load(file = sprintf("analyses_datasets/hyper/hyper_dataset_%d.RData", par.work[run]))
        
        # Storing data and prediction target:
        d <- datalist$data
        pred.target <- datalist$pred.target
        
        ### Model ##############################################################
        tic(run)
        RESULT <- gfBLUPold::siBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, t.RF = 0.75, verbose = FALSE, sepExp = FALSE)
        toc(log = TRUE)
        ########################################################################
        
        RESULT$preds <- RESULT$preds[match(pred.target$G, names(RESULT$preds))]
        
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
                                                     pen = pen,
                                                     comptimes = comptimes),
                                 extra = extra))
      
      names(worker.result) <- sprintf("worker_%d", i)
      return(worker.result)
      
    })
  doParallel::stopImplicitCluster()
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
  write.csv(results, sprintf("analyses_results/hyper/5%s_hyper_results_siblup_RF075_%s.csv", lab, CV))
  
  list.save(extra, sprintf("analyses_results/hyper/5%s_hyper_extra_results_siblup_RF075_%s.RData", lab, CV))
  
}




