h2y = "05"
comm = "05"
h2s = "05"
sim = 1

set.seed(1997)

load("genotypes/K_sim.RData")

# Loading simulated dataset:
datalist <- rlist::list.load(file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))

# Storing data and prediction target:
d <- datalist$data.real
pred.target <- datalist$pred.target

d.train <- droplevels(na.omit(d))
sec <- names(d[2:(ncol(d) - 1)])
foc <- names(d)[ncol(d)]

temp <- glfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 1, verbose = F)
d.train.RF <- cbind(temp$data.RF, d.train[foc])
d.RF <- d[names(d.train.RF)]
sec.RF <- names(d.RF[2:(ncol(d.RF) - 1)])

folds <- glfBLUP::createFolds(genos = unique(as.character(d.train.RF$G)))
tempG <- glfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)],
                                         folds = folds, what = "genetic",
                                         dopar = FALSE, verbose = F)
tempE <- glfBLUP::regularizedCorrelation(data = d.train.RF[c("G", sec.RF)],
                                         folds = folds, what = "residual",
                                         dopar = FALSE, verbose = F)
Rg.RF.reg <- tempG$optCor

FM.fit <- glfBLUP::factorModel(data = d.train.RF[c("G", sec.RF)], cormat = Rg.RF.reg, what = "genetic", verbose = F)

D <- sqrt(diag(tempG$Sg))
L.cov <- diag(D) %*% FM.fit$loadings
PSI.cov <- outer(D, D) * FM.fit$uniquenesses

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

selection <- glfBLUP::factorSelect(CV1.d.final, procedure = "leaps", verbose = F)

CV1.temp <- glfBLUP::glfBLUP(data = CV1.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)
CV2.temp <- glfBLUP::glfBLUP(data = CV2.d.final, selection = selection, K = K, sepExp = FALSE, verbose = F)

CV1.acc <- cor(pred.target, CV1.temp$preds[match(names(pred.target), names(CV1.temp$preds))])
CV2.acc <- cor(pred.target, CV2.temp$preds[match(names(pred.target), names(CV2.temp$preds))])

x1 <- read.csv(sprintf("p800/results/h2s%s/3a_p800_results_glfblup_CV1_h2y%s_comm%s_h2s%s.csv", h2s, h2y, comm, h2s))
x2 <- read.csv(sprintf("p800/results/h2s%s/3b_p800_results_glfblup_CV2_h2y%s_comm%s_h2s%s.csv", h2s, h2y, comm, h2s))

cat(sprintf("\nglfBLUP accuracy reproduced correctly for CV1: %s\n", all.equal(CV1.acc, x1$acc[sim])))
cat(sprintf("glfBLUP accuracy reproduced correctly for CV2: %s\n\n", all.equal(CV2.acc, x2$acc[sim])))


