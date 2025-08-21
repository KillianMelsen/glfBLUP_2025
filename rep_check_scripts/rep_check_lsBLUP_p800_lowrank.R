# This script is a copy of the script at `p800_lowrank/analyses/lsBLUP.R`.
# The purpose of this script is to provide an annotated example of how to perform
# a reproducibility spot check. All code that is not needed for the spot check has
# been commented out.
# Please carefully read the instructions below, and make sure you manually define
# the variables as explained below.

# This script produces all intermediate results for the lsBLUP CV1 and CV2 analyses on
# the simulated datasets with low-rank residual structure (p800_lowrank).
#
# The script is structured as follows:
# - There are two outer loops: one for the CVs (CV1 and CV2), and a second one
#   for the secondary feature heritabilities ("05", "07", "09").
# - Then there is a parallel loop dividing the 3 x 5 = 15 combinations of the
#   communalities ("02", "05", "08") and focal trait heritabilities ("01", "03", "05", "07", "09")
#   over 15 parallel workers.
# - The combination of communality and focal trait heritability produced by each worker
#   are defined on lines 71 and 72.
# - A seed is then set inside each worker and the 100 replicates for the specific
#   combination of CV, secondary feature heritability, communality, and focal trait
#   heritability are produced in a final inner loop starting on line 83.
#
# Suppose you want to reproduce the result of replication 27 for `CV = "CV1"`, `lab = "a"`,
# `h2s = "05"`, `h2y = "05"`, and `comm = "02"`. The first step is to set these values
# manually instead of relying on the two outer loops (CV and h2s), and assignments on
# lines 71 (comm) and 72 (h2y). Then we set `n.sim <- 27` instead of setting it to
# 100 on line 75. Then running the inner loop will produce the correct accuracy for
# replication 27 as the 27th value in `acc`.
# This can be compared to the 27th accuracy stored in the appropriate intermediate
# results csv file:
#
# all.equal(read.csv(sprintf("p800_lowrank/results/h2s%s/11%s_p800_lowrank_results_lsblup_%s_h2y%s_comm%s_h2s%s.csv", h2s, lab, CV, h2y, comm, h2s))$acc[27],
#           acc[27])
#
# [1] TRUE

# Define these manually:
CV = "CV1"
lab = "a"
h2s = "05"
h2y = "05"
comm = "02"
n.sim <- 27

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
# tic("lsBLUP")
# for (CV in CVs) {
  # for (h2s in h2.sec) {
    
    # tic(sprintf("lsBLUP p800_lowrank combi %d / %d, (CV = %s, h2s = %s)", combi, combis, CV, h2s))
    
    # cl <- parallel::makeCluster(15, outfile = sprintf("logs/lsBLUP_sim_p800_lowrank_h2s%s_%s.txt", h2s, CV))
    # registerDoParallel(cl)
    
    # invisible(
      # foreach::foreach(i = 1:nrow(par.combis), .packages = c("rlist", "tictoc")) %dopar% {
        
        # comm <- par.combis[i, "comm"]
        # h2y <- par.combis[i, "h2y"]
        
        # Number of simulated datasets to load:
        # n.sim <- 100
        
        # Setting up result storage:
        acc <- numeric(n.sim)
        extra <- vector("list", n.sim)
        
        # Running (SET SEED IN EACH PARALLEL WORKER!):
        set.seed(1997)
        for (sim in 1:n.sim) {
          
          cat(sprintf("Running lsBLUP on %s p800_lowrank dataset %d / %d, (h2s = %s, comm = %s, h2y = %s), combi %d / %d\n",
                      CV, sim, n.sim, h2s, comm, h2y, combi, combis))
          
          # Loading simulated dataset:
          datalist <- list.load(file = sprintf("p800_lowrank/datasets/p800_lowrank_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
          
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

cat(paste0(all.equal(read.csv(sprintf("p800_lowrank/results/h2s%s/11%s_p800_lowrank_results_lsblup_%s_h2y%s_comm%s_h2s%s.csv", h2s, lab, CV, h2y, comm, h2s))$acc[27],
                     acc[27]), "\n"))
        
# We don't need the code below for the reproducibility spot check.
#         # Retrieve computational times:
#         tictoc.logs <- tic.log(format = FALSE)
#         tic.clearlog()
#         comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
#         
#         # Collect results:
#         results <- data.frame(acc = acc,
#                               comptimes = comptimes)
#         
#         
#         # Making correct CV label:
#         if (CV == "CV1") {
#           lab <- "a"
#         } else if (CV == "CV2") {
#           lab <- "b"
#         }
#         
#         # Export results:
#         write.csv(results, sprintf("p800_lowrank/results/h2s%s/11%s_p800_lowrank_results_lsblup_%s_h2y%s_comm%s_h2s%s.csv",
#                                    h2s, lab, CV, h2y, comm, h2s))
#         list.save(extra, file = sprintf("p800_lowrank/results/h2s%s/11%s_p800_lowrank_extra_results_lsblup_%s_h2y%s_comm%s_h2s%s.RData",
#                                         h2s, lab, CV, h2y, comm, h2s))
#         
#       })
#     doParallel::stopImplicitCluster()
#     parallel::stopCluster(cl)
#     toc()
#     combi <- combi + 1
#   }
# }
# toc()

