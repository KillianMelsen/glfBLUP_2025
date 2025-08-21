# This script produces all intermediate results for the univariate analyses on
# the simulated datasets with random residual structure (p800).
#
# The script is structured as follows:
# - There is a single outer loop for the secondary feature heritabilities ("05", "07", "09").
# - Then there is a parallel loop dividing the 3 x 5 = 15 combinations of the
#   communalities ("02", "05", "08") and focal trait heritabilities ("01", "03", "05", "07", "09")
#   over 15 parallel workers.
# - The combination of communality and focal trait heritability produced by each worker
#   are defined on lines 67 and 68.
# - A seed is then set inside each worker and the 100 replicates for the specific
#   combination of CV, secondary feature heritability, communality, and focal trait
#   heritability are produced in a final inner loop starting on line 78.
#
# Suppose you want to reproduce the result of replication 41 for `h2s = "05"`,
# `h2y = "07"`, and `comm = "05"`. The first step is to set these values
# manually instead of relying on the outer loop (h2s), and assignments on
# lines 67 (comm) and 68 (h2y). Then we set `n.sim <- 41` instead of setting it to
# 100 on line 71. Then running the inner loop will produce the correct accuracy for
# replication 41 as the 41st value in `acc`.
# This can be compared to the 41st accuracy stored in the appropriate intermediate
# results csv file:
#
# all.equal(read.csv(sprintf("p800/results/h2s%s/2_p800_results_univariate_h2y%s_comm%s_h2s%s.csv", h2s, h2y, comm, h2s))$acc[41],
#           acc[41])
#
# [1] TRUE

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


