prep <- "nosplines"
source("helper_functions/Estimate_gcor_prediction.R")
load("genotypes/K_hyper.RData")
CV = "CV2"
run = 1
set.seed(1997)

# Loading hyperspectral dataset:
datalist <- list.load(file = sprintf("hyper_1415HEAT/datasets/%s/hyper_dataset_%d.RData", prep, run))

# Storing data and prediction target:
# 9 feb is last day of VEG, 25 feb is heading, 10 march is start of grain filling:
dates <- c("150414")
d <- datalist$data
select <- which(substr(names(d), 7, 12) %in% dates)
d <- d[c(1, select, ncol(d))]
pred.target <- datalist$pred.target

# Subsetting K (only really happens for the first dataset...):
K <- K[unique(d$G), unique(d$G)]

### Model ##############################################################
RESULT <- gfBLUPold::lsBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, sepExp = FALSE, t.RF = 0.95, verbose = FALSE)
########################################################################

RESULT$preds <- RESULT$preds[match(pred.target$G, names(RESULT$preds))]

#### Runcie & Cheng 2019 correction ----------------------------------
temp <- estimate_gcor(data = data.frame(ID = pred.target$G,
                                        obs = pred.target$pred.target,
                                        pred = RESULT$preds),
                      Knn = K[pred.target$G, pred.target$G],
                      method = "MCMCglmm",
                      normalize = T)
acc <- temp["g_cor"]

x <- read.csv(sprintf("hyper_1415HEAT/results/%s/11%s_hyper_results_lsblup_%sVEG.csv", prep, "b", CV))

cat(sprintf("\nlsBLUP accuracy reproduced correctly: %s\n\n", all.equal(as.numeric(acc), x$acc[run])))
