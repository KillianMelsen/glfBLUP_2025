CV = "CV2"
h2y = "05"
comm = "02"
h2s = "09"
sim = 1
lab = "b"

set.seed(1997)

load("genotypes/K_sim.RData")
datalist <- rlist::list.load(file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
d <- datalist$data.real
pred.target <- datalist$pred.target

RESULT <- gfBLUPold::lsBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, sepExp = FALSE, t.RF = 1)

RESULT$preds <- RESULT$preds[match(names(pred.target), names(RESULT$preds))]
acc <- cor(RESULT$preds, pred.target)

x <- read.csv(sprintf("p800/results/h2s%s/11%s_p800_results_lsblup_%s_h2y%s_comm%s_h2s%s.csv",
                      h2s, lab, CV, h2y, comm, h2s))

cat(sprintf("lsBLUP accuracy reproduced correctly: %s", all.equal(acc, x$acc[sim])))
