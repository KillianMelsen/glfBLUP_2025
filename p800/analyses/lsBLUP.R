# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!
# Should take about 1 hour and 30 minutes

# Loading libraries:
library(rlist)
library(tictoc)
library(doParallel)
library(gfBLUPold)

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

combis <- length(h2.sec) * length(CVs)
combi <- 1

par.combis <- data.frame(comm = rep(comms, each = length(h2.foc)),
                         h2y = rep(h2.foc, length(comms)))

# All simulated data for both CV1 and CV2:
tic("lsBLUP")
for (CV in CVs) {
  for (h2s in h2.sec) {
    
    tic(sprintf("lsBLUP p800 combi %d / %d, (CV = %s, h2s = %s)", combi, combis, CV, h2s))
    
    cl <- parallel::makeCluster(15, outfile = sprintf("logs/lsBLUP_sim_p800_h2s%s_%s.txt", h2s, CV))
    registerDoParallel(cl)
    
    invisible(
    foreach::foreach(i = 1:nrow(par.combis), .packages = c("rlist", "tictoc")) %dopar% {
      
      comm <- par.combis[i, "comm"]
      h2y <- par.combis[i, "h2y"]
      
      # Number of simulated datasets to load:
      n.sim <- 100
      
      # Setting up result storage:
      acc <- numeric(n.sim)
      extra <- vector("list", n.sim)
      
      # Running (SET SEED IN EACH PARALLEL WORKER!):
      set.seed(1997)
      for (sim in 1:n.sim) {
        
        cat(sprintf("Running lsBLUP on %s p800 dataset %d / %d, (h2s = %s, comm = %s, h2y = %s), combi %d / %d\n",
                    CV, sim, n.sim, h2s, comm, h2y, combi, combis))
        
        # Loading simulated dataset:
        datalist <- list.load(file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
        
        # Storing data and prediction target:
        d <- datalist$data.real
        pred.target <- datalist$pred.target
        
        ### Model ############################################################
        tic(sim)
        RESULT <- gfBLUPold::lsBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, sepExp = FALSE, t.RF = 1)
        toc(log = TRUE)
        ######################################################################
        
        RESULT$preds <- RESULT$preds[match(names(pred.target), names(RESULT$preds))]
        
        acc[sim] <- cor(RESULT$preds, pred.target)
        
        if (length(RESULT) == 1) {
          extra[[sim]] <- "All coefficients were 0, used univariate model..."
        } else {
          extra[[sim]] <- list(coefs = RESULT$coefs,
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
      results <- data.frame(acc = acc,
                            comptimes = comptimes)
      
      
      # Making correct CV label:
      if (CV == "CV1") {
        lab <- "a"
      } else if (CV == "CV2") {
        lab <- "b"
      }
      
      # Export results:
      write.csv(results, sprintf("p800/results/h2s%s/11%s_p800_results_lsblup_%s_h2y%s_comm%s_h2s%s.csv",
                                 h2s, lab, CV, h2y, comm, h2s))
      list.save(extra, file = sprintf("p800/results/h2s%s/11%s_p800_extra_results_lsblup_%s_h2y%s_comm%s_h2s%s.RData",
                                      h2s, lab, CV, h2y, comm, h2s))
      
    })
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    toc()
    combi <- combi + 1
  }
}
toc()

