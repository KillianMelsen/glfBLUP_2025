# !!! IMPORTANT: SET MKL_NUMTHREADS=1 and MKL_DYNAMIC=TRUE in ~/.profile !!!
# Runtime: ~ 330s for 250 datasets.

# Loading libraries:
library(rlist)
library(tictoc)
library(doParallel)

# Setting seed:
set.seed(1997)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading hyperspectral kinship:
load("genotypes/H_hyper.RData")

# Hyperspectral data:
tic("Phenomic")

n.datasets <- 250
n.cores <- parallel::detectCores() - 2
work <- split(1:n.datasets, ceiling(seq_along(1:n.datasets) / ceiling(n.datasets / n.cores)))
cl <- parallel::makeCluster(n.cores, outfile = "logs/phenomic_hyper.txt")
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
    datalist <- list.load(file = sprintf("hyper/datasets/hyper_dataset_%d.RData", par.work[run]))
    
    # Storing data and prediction target:
    d <- datalist$data
    pred.target <- datalist$pred.target
    
    ### Model ##############################################################
    tic(run)
    RESULT <- rrBLUP::kin.blup(data = d, geno = "G", pheno = "Y", K = H)
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
write.csv(results, "hyper/results/0_hyper_results_phenomic.csv")




