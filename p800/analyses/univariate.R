# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!

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

combis <- length(h2.sec)
combi <- 1

par.combis <- data.frame(comm = rep(comms, each = length(h2.foc)),
                         h2y = rep(h2.foc, length(comms)))

# All simulated data (no CV's because univariate):
tic("Univariate")
for (h2s in h2.sec) {
  
  tic(sprintf("Univariate combi %d / %d, (h2s = %s)", combi, combis, h2s))
  
  cl <- parallel::makeCluster(15, outfile = sprintf("logs/univariate_sim_p800_h2s%s.txt", h2s))
  registerDoParallel(cl)
  
  invisible(
  foreach::foreach(i = 1:nrow(par.combis), .packages = c("rlist", "tictoc")) %dopar% {
    
    comm <- par.combis[i, "comm"]
    h2y <- par.combis[i, "h2y"]
    
    # Number of simulated datasets to load:
    n.sim <- 100
    
    # Setting up result storage:
    acc <- numeric(n.sim)
    
    # Running (SET SEED IN EACH PARALLEL WORKER!):
    set.seed(1997)
    for (sim in 1:n.sim) {
      
      cat(sprintf("Running univariate on p800 dataset %d / %d, (h2s = %s, comm = %s, h2y = %s), combi %d / %d\n",
                  sim, n.sim, h2s, comm, h2y, combi, combis))
      
      # Loading simulated dataset:
      datalist <- list.load(file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
      
      # Storing data and prediction target:
      d <- glfBLUP:::genotypeMeans(datalist$data.real)
      d[,-1] <- scale(d[,-1])
      pred.target <- datalist$pred.target
      
      ### Model ##############################################################
      tic(sim)
      RESULT <- rrBLUP::kin.blup(data = d, geno = "G", pheno = "Y", K = K)
      toc(log = TRUE)
      ########################################################################
      
      acc[sim] <- cor(RESULT$g[datalist$test.set], pred.target)
      
    }
    
    # Retrieve computational times:
    tictoc.logs <- tic.log(format = FALSE)
    tic.clearlog()
    comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
    
    # Collect results:
    results <- data.frame(acc = acc,
                          comptimes = comptimes)
    
    # Export results:
    write.csv(results, sprintf("p800/results/h2s%s/2_p800_results_univariate_h2y%s_comm%s_h2s%s.csv",
                               h2s, h2y, comm, h2s))
    
  })
  doParallel::stopImplicitCluster()
  parallel::stopCluster(cl)
  toc()
  combi <- combi + 1
}
toc()


