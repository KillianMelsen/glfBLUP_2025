prep <- "splines"
source("helper_functions/Estimate_gcor_prediction.R")
load("genotypes/K_hyper.RData")
set.seed(1997)
run = 1
datalist <- rlist::list.load(file = sprintf("hyper_1415B5IR/datasets/%s/hyper_dataset_%d.RData", prep, run))

# Storing data and prediction target:
d <- datalist$data
pred.target <- datalist$pred.target

# Subsetting K (only really happens for the first dataset...):
K <- K[unique(d$G), unique(d$G)]

### Model ##############################################################

### 1. Make training data and store feature/focal trait names ------------------------------------------------------------------------
d.train <- droplevels(na.omit(d))
sec <- names(d[2:(ncol(d) - 1)])
foc <- names(d)[ncol(d)]

### 2. Redundancy filter the secondary features using training data only -------------------------------------------------------------
temp <- glfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 0.95, verbose = FALSE)
d.train.RF <- cbind(temp$data.RF, d.train[foc])
d.RF <- d[names(d.train.RF)]
sec.RF <- names(d.RF[2:(ncol(d.RF) - 1)])

### 3. Regularization ----------------------------------------------------------------------------------------------------------------
folds <- glfBLUP::createFolds(genos = unique(as.character(d.train.RF$G)))
tempG <- glfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "genetic", dopar = FALSE, verbose = FALSE)
tempE <- glfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)], folds = folds, what = "residual", dopar = FALSE, verbose = FALSE)
Rg.RF.reg <- tempG$optCor

### 4. Fitting factor model ----------------------------------------------------------------------------------------------------------
# data is only used to determine the sample size for the MP-bound. what specifies that it's a genetic correlation matrix, so the
# number of training genotypes should be used, and not the number of training individuals (= genotypes * replicates).
FM.fit <- glfBLUP::factorModel(data = d.train.RF[c("G", sec.RF)], cormat = Rg.RF.reg, what = "genetic", verbose = FALSE)

#### 5. Getting factor scores (also for the test set in CV2!) ------------------------------------------------------------------------
# Loadings and uniquenesses were estimated on the correlation scale, but should be on the covariance scale for genetic-thomson scores:
D <- sqrt(diag(tempG$Sg)) # Getting standard deviations
L.cov <- diag(D) %*% FM.fit$loadings
PSI.cov <- outer(D, D) * FM.fit$uniquenesses

# CV1 Factor scores:
CV1.d.RF <- d.RF
CV1.d.RF[which(is.na(CV1.d.RF$Y)), 2:ncol(CV1.d.RF)] <- NA
CV1.F.scores <- glfBLUP::factorScores(data = CV1.d.RF[c("G", sec.RF)],
                                      loadings = L.cov,
                                      uniquenesses = PSI.cov,
                                      m = FM.fit$m,
                                      type = "genetic-thomson-repdiv",
                                      Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)

CV1.d.final <- cbind(CV1.F.scores, CV1.d.RF$Y)
names(CV1.d.final)[ncol(CV1.d.final)] <- "Y"
names(CV1.d.final)[1] <- "G"

# CV2 Factor scores:
# First recenter/rescale the training and test secondary data together:
d.RF[sec.RF] <- sapply(d.RF[sec.RF], scale)
CV2.F.scores <- glfBLUP::factorScores(data = d.RF[c("G", sec.RF)],
                                      loadings = L.cov,
                                      uniquenesses = PSI.cov,
                                      m = FM.fit$m,
                                      type = "genetic-thomson-repdiv",
                                      Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)

CV2.d.final <- cbind(CV2.F.scores, d.RF$Y)
names(CV2.d.final)[ncol(CV2.d.final)] <- "Y"
names(CV2.d.final)[1] <- "G"

#### 6. Selecting the relevant factors -----------------------------------------------------------------------------------------------
selection <- glfBLUP::factorSelect(CV1.d.final, procedure = "leaps", verbose = FALSE)

#### 7. Multi-trait genomic prediction -----------------------------------------------------------------------------------------------
CV1.temp <- glfBLUP::glfBLUP(data = CV1.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
CV2.temp <- glfBLUP::glfBLUP(data = CV2.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
########################################################################

CV1.acc <- cor(pred.target$pred.target, CV1.temp$preds[match(pred.target$G, names(CV1.temp$preds))])

#### Runcie & Cheng 2019 correction --------------------------------------
temp <- estimate_gcor(data = data.frame(ID = pred.target$G,
                                        obs = pred.target$pred.target,
                                        pred = CV2.temp$preds[match(pred.target$G, names(CV2.temp$preds))]),
                      Knn = K[pred.target$G, pred.target$G],
                      method = "MCMCglmm",
                      normalize = T)
CV2.acc <- temp["g_cor"]

x1 <- read.csv(sprintf("hyper_1415B5IR/results/%s/3a_hyper_results_glfblup_CV1.csv", prep))
x2 <- read.csv(sprintf("hyper_1415B5IR/results/%s/3b_hyper_results_glfblup_CV2.csv", prep))

cat(sprintf("\nglfBLUP accuracy reproduced correctly for CV1: %s\n", all.equal(CV1.acc, x1$acc[run])))
cat(sprintf("glfBLUP accuracy reproduced correctly for CV2: %s\n\n", all.equal(as.numeric(CV2.acc), x2$acc[run])))


