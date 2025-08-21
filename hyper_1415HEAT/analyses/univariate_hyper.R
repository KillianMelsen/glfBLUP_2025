# This script produces all intermediate results for the univariate analyses on
# the hyperspectral HEAT dataset.
#
# The script is structured as follows:
# - There is a single parallel loop dividing the 250 replications over 18 parallel
#   workers.
# - The replications produced by each worker are given by the `par.work`
#   variable on line 54.
# - A seed is set inside each worker and the replicates of the worker are produced
#   inside the inner loop starting on line 61.
#
# Suppose you want to reproduce the result of replication 246.
# Looking at the `work` variable, replication 246 is handled by
# the 18th worker. So we can manually set `i <- 18` instead of running the parallel
# for loop. We also note that replication 246 is the 8th replication that is
# produced by the 18th worker, so in the inner loop we do not go from
# `1:length(par.work)`, but `1:8`. Then running the loop
# will produce the univariate accuracy of replication 246 as the 8th value in `acc`.
# This can be compared to the 246th accuracy stored in the appropriate intermediate
# results csv file:
#
# all.equal(read.csv("hyper_1415HEAT/results/2_hyper_results_univariate.csv")$acc[246],
#           acc[8])
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

# Loading kinship:
load("genotypes/K_hyper.RData")

# Hyperspectral data:
tic("Univariate")

n.datasets <- 250
n.cores <- 18
work <- split(1:n.datasets, ceiling(seq_along(1:n.datasets) / ceiling(n.datasets / n.cores)))
cl <- parallel::makeCluster(n.cores, outfile = "logs/univariate_hyper.txt")
doParallel::registerDoParallel(cl)

invisible(
par.results <- foreach::foreach(i = 1:length(work), .packages = c("rlist", "tictoc"), .combine = "c") %dopar% {
  
  par.work <- work[[i]]
  set.seed(1997)
  
  # Setting up result storage:
  acc <- numeric(length(par.work))
  
  # Running:
  for (run in 1:length(par.work)) {

    # Loading hyperspectral dataset:
    datalist <- list.load(file = sprintf("hyper_1415HEAT/datasets/splines/hyper_dataset_%d.RData", par.work[run]))
    
    # Storing data and prediction target:
    d <- glfBLUP:::genotypeMeans(datalist$data)
    d[,-1] <- scale(d[,-1])
    pred.target <- datalist$pred.target
    
    # Subsetting K (only really happens for the first dataset...):
    K <- K[unique(d$G), unique(d$G)]
    
    ### Model ##############################################################
    tic(run)
    RESULT <- rrBLUP::kin.blup(data = d, geno = "G", pheno = "Y", K = K)
    toc(log = TRUE)
    ########################################################################
    
    acc[run] <- cor(RESULT$g[datalist$test.set], pred.target$pred.target)
    
  }
  
  # Retrieve computational times:
  tictoc.logs <- tic.log(format = FALSE)
  tic.clearlog()
  comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
  
  # Collect results:
  result <- list(data.frame(acc = acc,
                            comptimes = comptimes))
  
  names(result) <- sprintf("worker_%d", i)
  return(result)
  
})
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
toc()

# Restructuring parallel results for saving:
acc <- numeric()
comptimes <- numeric()
for (j in 1:length(work)) {
  
  acc <- c(acc, par.results[[sprintf("worker_%d", j)]]$acc)
  comptimes <- c(comptimes, par.results[[sprintf("worker_%d", j)]]$comptimes)
  
}

results <- data.frame(acc = acc,
                      comptimes = comptimes)

# Export results:
write.csv(results, "hyper_1415HEAT/results/2_hyper_results_univariate.csv")




