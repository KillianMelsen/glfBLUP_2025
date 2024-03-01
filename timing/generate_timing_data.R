start <- Sys.time()

# Loading libraries:
library(rlist)
library(tictoc)

# Setting seed:
set.seed(1997)

# Loading kinship:
load("genotypes/K_sim.RData"); rm(M)

# Simulated genetic parameters:
h2s <- "09"
comm <- "08"
h2y <- "05"

# Total number of secondary features:
ps <- seq(100, 1400, 100)

p <- ps[1]
for (p in ps) {
  
  # Number of latent factors (noise + signal):
  m <- 10
  n.LSF <- floor(m / 2)
  n.LNF <- m - n.LSF
  S.per.LF <- p / m
  
  # Number of replicates:
  r <- 3
  
  # Determining numeric values for the genetic parameters:
  comm.num <- as.numeric(comm) / 10
  
  sg2.s <- as.numeric(h2s) / 10
  se2.s <- 1 - sg2.s
  
  sg2.y <- as.numeric(h2y) / 10
  se2.y <- 1 - sg2.y
  
  # Running simulations:
  tic(sprintf("p = %s", p))
  datalist <- gfBLUP::GFAsim(K = K, r = r, n.LSF = n.LSF, n.LNF = n.LNF,
                             LSF.rho = 0, LNF.rho = 0, LSNF.rho = 0,
                             S.per.LF = S.per.LF, Y.psi = 1 - comm.num,
                             L.min = 0.3, L.max = 0.8,
                             S.sg2 = sg2.s, S.se2 = se2.s,
                             Y.sg2 = sg2.y, Y.se2 = se2.y)
  
  datalist$simParams <- list(r = r, n.LSF = n.LSF, n.LNF = n.LNF,
                             LSF.rho = 0, LNF.rho = 0, LSNF.rho = 0,
                             S.per.LF = S.per.LF, Y.psi = 1 - comm.num,
                             L.min = 0.3, L.max = 0.8,
                             S.sg2 = sg2.s, S.se2 = se2.s,
                             Y.sg2 = sg2.y, Y.se2 = se2.y)
  
  list.save(datalist, file = sprintf("timing/datasets/timing_p%s.RData", p))
  
  cat(sprintf("Generation of dataset for p = %s done!\n\n", p))
  toc()
  
}

end <- Sys.time()
end - start










