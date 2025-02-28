CV = "CV1"
h2y = "03"
comm = "05"
h2s = "07"
sim = 1
lab = "a"

set.seed(1997)

load("genotypes/K_sim.RData")
datalist <- rlist::list.load(file = sprintf("p800/datasets/p800_h2y%s_comm%s_h2s%s_dataset_%d.RData", h2y, comm, h2s, sim))
d <- datalist$data.real
pred.target <- datalist$pred.target

RESULT <- gfBLUPold::siBLUP(d = d, K = K, CV = CV, do.parallel = FALSE, t.RF = 1, verbose = FALSE)

RESULT$preds <- RESULT$preds[match(names(pred.target), names(RESULT$preds))]

acc <- cor(RESULT$preds, pred.target)

x <- read.csv(sprintf("p800/results/h2s%s/5%s_p800_results_siblup_%s_h2y%s_comm%s_h2s%s_1to50.csv",
                      h2s, lab, CV, h2y, comm, h2s))

cat(sprintf("siBLUP accuracy reproduced correctly: %s", all.equal(acc, x$acc[sim])))
