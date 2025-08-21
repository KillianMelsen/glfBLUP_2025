# This script uses the pre-processed data to generate the datafiles for the 250
# replications for the hyperspectral B5IR data.
#
# To check reproducibility of any of the generated datafiles, set the `n.datasets`
# variable to the replication you would like to check.
#
# Suppose we want to check the datafile for replication 11. We first manually set `n.datasets <- 11`
# instead of setting it to 250 on line 28. We then temporarily comment out the code
# that saves the datafiles to disk on lines 90 and 101 to avoid overwriting. We can then
# simply run the loop.
# We can then check the data stored in the variables `datalist.splines` and `datalist.VEGsplines`
# against the datafiles already stored on disk:
#
# all.equal(rlist::list.load(sprintf("hyper_1415B5IR/datasets/splines/hyper_dataset_%d.RData", run)),
#           datalist.splines)
# [1] TRUE
#
# all.equal(rlist::list.load(sprintf("hyper_1415B5IR/datasets/VEGsplines/hyper_dataset_%d.RData", run)),
#           datalist.VEGsplines)
# [1] TRUE

library(tictoc)

# Setting seed:
set.seed(1997)

# Number of datasets:
n.datasets <- 250

# Loading data:
pseudoCRD.splines <- readRDS("hyper_1415B5IR/data_generation/pseudoCRD_splines.rds")
pseudoCRD.VEGsplines <- readRDS("hyper_1415B5IR/data_generation/pseudoCRD_VEGsplines.rds")

# 38 trials:
trials  <- unique(pseudoCRD.splines$trial)

tic()
for (run in 1:n.datasets) {
  # 12 trials go to test set, 26 to train set:
  test <- sample(trials, 12)
  train <- trials[!trials %in% test]
  
  # Reloading data:
  pseudoCRD.splines <- readRDS("hyper_1415B5IR/data_generation/pseudoCRD_splines.rds")
  pseudoCRD.VEGsplines <- readRDS("hyper_1415B5IR/data_generation/pseudoCRD_VEGsplines.rds")
  
  # Fixing colnames:
  colnames(pseudoCRD.splines)[4:(ncol(pseudoCRD.splines) - 2)] <-
    substr(colnames(pseudoCRD.splines)[4:(ncol(pseudoCRD.splines) - 2)], 2, 13)
  
  colnames(pseudoCRD.VEGsplines)[4:(ncol(pseudoCRD.VEGsplines) - 2)] <-
    substr(colnames(pseudoCRD.VEGsplines)[4:(ncol(pseudoCRD.VEGsplines) - 2)], 2, 13)
  
  pseudoCRD.splines.test <- pseudoCRD.splines[which(pseudoCRD.splines$trial %in% test),]
  pseudoCRD.splines.train <- pseudoCRD.splines[which(pseudoCRD.splines$trial %in% train),]
  
  pseudoCRD.VEGsplines.test <- pseudoCRD.VEGsplines[which(pseudoCRD.VEGsplines$trial %in% test),]
  pseudoCRD.VEGsplines.train <- pseudoCRD.VEGsplines[which(pseudoCRD.VEGsplines$trial %in% train),]
  
  # Scaling secondary data separately for training and test set:
  pseudoCRD.splines.test[, 4:(ncol(pseudoCRD.splines.test) - 1)] <-
    sapply(pseudoCRD.splines.test[, 4:(ncol(pseudoCRD.splines.test) - 1)], scale)
  pseudoCRD.splines.train[, 4:(ncol(pseudoCRD.splines.train) - 1)] <-
    sapply(pseudoCRD.splines.train[, 4:(ncol(pseudoCRD.splines.train) - 1)], scale)
  
  pseudoCRD.VEGsplines.test[, 4:(ncol(pseudoCRD.VEGsplines.test) - 1)] <-
    sapply(pseudoCRD.VEGsplines.test[, 4:(ncol(pseudoCRD.VEGsplines.test) - 1)], scale)
  pseudoCRD.VEGsplines.train[, 4:(ncol(pseudoCRD.VEGsplines.train) - 1)] <-
    sapply(pseudoCRD.VEGsplines.train[, 4:(ncol(pseudoCRD.VEGsplines.train) - 1)], scale)
  
  # Recreating full datasets after scaling:
  pseudoCRD.splines <- rbind(pseudoCRD.splines.test, pseudoCRD.splines.train)
  pseudoCRD.VEGsplines <- rbind(pseudoCRD.VEGsplines.test, pseudoCRD.VEGsplines.train)
  
  # Setting gy_adjusted (reps from pseudo-CRD) to NA for test genotypes:
  pseudoCRD.splines[which(pseudoCRD.splines$gid %in% unique(pseudoCRD.splines.test$gid)),
                      "gy_adjusted"] <- NA
  pseudoCRD.VEGsplines[which(pseudoCRD.VEGsplines$gid %in% unique(pseudoCRD.VEGsplines.test$gid)),
                    "gy_adjusted"] <- NA
  
  # Saving the data (time splines):
  datalist.splines <- list(data = as.data.frame(pseudoCRD.splines[, 3:624]),
                           pred.target = unique(as.data.frame(pseudoCRD.splines.test[, c(3, 625)])),
                           test.set = unique(pseudoCRD.splines$gid[which(is.na(pseudoCRD.splines$gy_adjusted))]),
                           train.set = unique(pseudoCRD.splines$gid[which(!is.na(pseudoCRD.splines$gy_adjusted))]))
  
  names(datalist.splines$pred.target) <- c("G", "pred.target")
  names(datalist.splines$data)[1] <- "G"
  names(datalist.splines$data)[622] <- "Y"
  rlist::list.save(datalist.splines, file = sprintf("hyper_1415B5IR/datasets/splines/hyper_dataset_%d.RData", run))
  
  # Saving the data (VEG time splines):
  datalist.VEGsplines <- list(data = as.data.frame(pseudoCRD.VEGsplines[, 3:252]),
                              pred.target = unique(as.data.frame(pseudoCRD.VEGsplines.test[, c(3, 253)])),
                              test.set = unique(pseudoCRD.VEGsplines$gid[which(is.na(pseudoCRD.VEGsplines$gy_adjusted))]),
                              train.set = unique(pseudoCRD.VEGsplines$gid[which(!is.na(pseudoCRD.VEGsplines$gy_adjusted))]))
  
  names(datalist.VEGsplines$pred.target) <- c("G", "pred.target")
  names(datalist.VEGsplines$data)[1] <- "G"
  names(datalist.VEGsplines$data)[250] <- "Y"
  rlist::list.save(datalist.VEGsplines, file = sprintf("hyper_1415B5IR/datasets/VEGsplines/hyper_dataset_%d.RData", run))
}
toc()

