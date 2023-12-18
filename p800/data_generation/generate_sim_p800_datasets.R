# For execution on WSL using oneMKL:
# Add `export MKL_NUM_THREADS=3` and `export MKL_DYNAMIC=FALSE` to ~/.profile (assuming a system with 20 threads)
# 5 parallel processes using foreach * 3 MKL threads = 15 out of 20 threads total.
# Runtimes for 5 x 100 datasets:
# WSL MKL1:
# WSL MKL3:

# This scripts generates 100 datasets for each of the 45 h2s - comm - h2y combinations.
# These datasets contain 500 genotypes replicated twice, with 300 training genotypes and 200 test genotypes.
# There are 4 signal factors and 4 noise factors with 100 features each, for a total of p=800 features.
# As a result the data is high-dimensional as p > n > n_g.

# Signal factors are not correlated to each other.

# Working directory and packages:
start <- Sys.time()
wd <- "~/gfblup_methodology"
# wd <- "C:/Users/killi/Desktop/gfblup-methodological-paper"
setwd(wd)
library(doParallel)
library(tictoc)

# Simulated genetic parameters:
h2.sec <- c("05", "07", "09")
comms <- c("02", "05", "08")
h2.foc <- c("01", "03", "05", "07", "09")

combis <- length(h2.sec) * length(comms)
combi <- 1

for (h2s in h2.sec) {
  for (comm in comms) {
    
    cl <- parallel::makeCluster(5, outfile = sprintf("logs/generate_sim_p800_h2s%s_comm%s.txt", h2s, comm))
    doParallel::registerDoParallel(cl)
    
    cat(sprintf("Starting parallel generation of 5 x 100 simulated datasets (h2s = %s, comm = %s, combination %d/%d)...\n",
                h2s, comm, combi, combis))
    tic(sprintf("h2s = %s, comm = %s", h2s, comm))
    
    foreach::foreach(i = 1:length(h2.foc)) %dopar% {
      
      # Getting h2F chr value:
      h2y <- h2.foc[i]
      
      # Loading libraries:
      library(rlist)
      
      # Setting seed and working directory at start of 100 simulated datasets:
      set.seed(1997)
      setwd(wd)
      
      # Loading kinship and marker data, setting number of simulations and replicates:
      load("K_sim.RData")
      n.sim <- 100
      r <- 2
      
      # Determining numeric values for the genetic parameters:
      comm.num <- as.numeric(comm) / 10
      
      sg2.s <- as.numeric(h2s) / 10
      se2.s <- 1 - sg2.s
      
      sg2.y <- as.numeric(h2y) / 10
      se2.y <- 1 - sg2.y
      
      # Running simulations:
      sim <- 1
      for (sim in 1:n.sim) {
        
        datalist <- gfBLUP::GFAsim(K = K, r = r, n.LSF = 4, n.LNF = 4,
                                   LSF.rho = 0, LNF.rho = 0, LSNF.rho = 0,
                                   S.per.LF = 100, Y.psi = 1 - comm.num,
                                   L.min = 0.3, L.max = 0.8,
                                   S.sg2 = sg2.s, S.se2 = se2.s,
                                   Y.sg2 = sg2.y, Y.se2 = se2.y)
        
        datalist$simParams <- list(r = r, n.LSF = 4, n.LNF = 4,
                                   LSF.rho = 0, LNF.rho = 0, LSNF.rho = 0,
                                   S.per.LF = 100, Y.psi = 1 - comm.num,
                                   L.min = 0.3, L.max = 0.8,
                                   S.sg2 = sg2.s, S.se2 = se2.s,
                                   Y.sg2 = sg2.y, Y.se2 = se2.y)
          
        list.save(datalist, file = sprintf("sim_p800/datasets/sim_p800_h2y%s_comm%s_h2s%s_dataset_%d.RData",
                                           h2y, comm, h2s, sim))
      }
    }
    toc()
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    combi <- combi + 1
  }
}

end <- Sys.time()
end - start


