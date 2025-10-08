# This script produces all intermediate results for the benchmark CV1 and CV2 analyses on
# analyses on the simulated datasets with low-rank residual structure (p800_lowrank).
#
# The script is structured as follows:
# - There are two outer loops: one for the CVs (CV1 and CV2), and a second one
#   for the secondary feature heritabilities ("05", "07", "09").
# - Then there is a parallel loop dividing the 3 x 5 = 15 combinations of the
#   communalities ("02", "05", "08") and focal trait heritabilities ("01", "03", "05", "07", "09")
#   over 15 parallel workers.
# - The combination of communality and focal trait heritability produced by each worker
#   are defined on lines 73 and 74.
# - A seed is then set inside each worker and the 100 replicates for the specific
#   combination of CV, secondary feature heritability, communality, and focal trait
#   heritability are produced in a final inner loop starting on line 84.
#
# Suppose you want to reproduce the result of replication 13 for `CV = "CV2"`, `lab = "b"`,
# `h2s = "07"`, `h2y = "07"`, and `comm = "08"`. The first step is to set these values
# manually instead of relying on the two outer loops (CV and h2s), and assignments on
# lines 73 (comm) and 74 (h2y). Then we set `n.sim <- 13` instead of setting it to
# 100 on line 77. Then running the inner loop will produce the correct accuracy for
# replication 13 as the 13th value in `acc`.
# This can be compared to the 13th accuracy stored in the appropriate intermediate
# results csv file:
#
# all.equal(read.csv(sprintf("p800_lowrank/results/h2s%s/1%s_p800_lowrank_results_benchmark_%s_h2y%s_comm%s_h2s%s.csv", h2s, lab, CV, h2y, comm, h2s))$acc[13],
#           acc[13])
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

# Loading helper functions:
source("helper_functions.R")

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
tic("Benchmark")
for (CV in CVs) {
  for (h2s in h2.sec) {
    
    tic(sprintf("Benchmark combi %d / %d, (CV = %s, h2s = %s)", combi, combis, CV, h2s))
    
    cl <- parallel::makeCluster(15, outfile = sprintf("logs/benchmark_sim_p800_lowrank_h2s%s_%s.txt", h2s, CV))
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
          
          cat(sprintf("Running benchmark on %s p800_lowrank dataset %d / %d, (h2s = %s, comm = %s, h2y = %s), combi %d / %d\n",
                      CV, sim, n.sim, h2s, comm, h2y, combi, combis))
          
          # Loading simulated dataset:
          datalist <- list.load(file = sprintf("p800_lowrank/datasets/p800_lowrank_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
          
          # Storing data and prediction target:
          d <- datalist$data.bm
          pred.target <- datalist$pred.target
          
          # Determining number of replicates in simulated dataset:
          r <- nrow(datalist$data.real) / length(unique(datalist$data.real$G))
          
          # Storing benchmark Vg and Ve and giving the correct names which I forgot during data generation...:
          Sg_bm <- datalist$Sg.bm
          Se_bm <- datalist$Se.bm
          colnames(Sg_bm) <- rownames(Sg_bm) <- c(paste0("F", 1:(ncol(Sg_bm) - 1)), "Y")
          colnames(Se_bm) <- rownames(Se_bm) <- c(paste0("F", 1:(ncol(Se_bm) - 1)), "Y")
          
          
          ### Model ############################################################
          if (CV == "CV1") {
            tic(sim)
            RESULT <- benchmark1(d = d, Sg = Sg_bm, Se = Se_bm, K = K, r = r)
            toc(log = TRUE)
          } else {
            tic(sim)
            RESULT <- Benchmark1_CV2(d = d, Sg = Sg_bm, Se = Se_bm, K = K, r = r)
            toc(log = TRUE)
          }
          ######################################################################
          
          acc[sim] <- cor(RESULT$preds, pred.target)
          
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
        write.csv(results, sprintf("p800_lowrank/results/h2s%s/1%s_p800_lowrank_results_benchmark_%s_h2y%s_comm%s_h2s%s.csv",
                                   h2s, lab, CV, h2y, comm, h2s))
        
      })
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    toc()
    combi <- combi + 1
  }
}
toc()

