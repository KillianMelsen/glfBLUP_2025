library(tictoc)

# Setting seed:
set.seed(1997)

# Number of datasets:
n.datasets <- 250

# Loading data:
pseudoCRD.nosplines <- readRDS("hyper_1415HEAT/data_generation/pseudoCRD_nosplines.rds")
pseudoCRD.splines <- readRDS("hyper_1415HEAT/data_generation/pseudoCRD_splines.rds")
# pseudoCRD.VEGsplines <- readRDS("hyper_1415HEAT/data_generation/pseudoCRD_VEGsplines.rds")

# 37 trials:
trials  <- unique(pseudoCRD.splines$trial)

tic()
for (run in 1:n.datasets) {
  # 12 trials go to test set, 25 to train set:
  test <- sample(trials, 12)
  train <- trials[!trials %in% test]
  
  # Reloading data:
  pseudoCRD.nosplines <- readRDS("hyper_1415HEAT/data_generation/pseudoCRD_nosplines.rds")
  pseudoCRD.splines <- readRDS("hyper_1415HEAT/data_generation/pseudoCRD_splines.rds")
  # pseudoCRD.VEGsplines <- readRDS("hyper_1415HEAT/data_generation/pseudoCRD_VEGsplines.rds")
  
  pseudoCRD.nosplines.test <- pseudoCRD.nosplines[which(pseudoCRD.nosplines$trial %in% test),]
  pseudoCRD.nosplines.train <- pseudoCRD.nosplines[which(pseudoCRD.nosplines$trial %in% train),]
  
  pseudoCRD.splines.test <- pseudoCRD.splines[which(pseudoCRD.splines$trial %in% test),]
  pseudoCRD.splines.train <- pseudoCRD.splines[which(pseudoCRD.splines$trial %in% train),]
  
  # pseudoCRD.VEGsplines.test <- pseudoCRD.VEGsplines[which(pseudoCRD.VEGsplines$trial %in% test),]
  # pseudoCRD.VEGsplines.train <- pseudoCRD.VEGsplines[which(pseudoCRD.VEGsplines$trial %in% train),]
  
  # Scaling secondary data separately for training and test set:
  pseudoCRD.nosplines.test[, 4:(ncol(pseudoCRD.nosplines.test) - 1)] <-
    sapply(pseudoCRD.nosplines.test[, 4:(ncol(pseudoCRD.nosplines.test) - 1)], scale)
  pseudoCRD.nosplines.train[, 4:(ncol(pseudoCRD.nosplines.train) - 1)] <-
    sapply(pseudoCRD.nosplines.train[, 4:(ncol(pseudoCRD.nosplines.train) - 1)], scale)
  
  pseudoCRD.splines.test[, 4:(ncol(pseudoCRD.splines.test) - 1)] <-
    sapply(pseudoCRD.splines.test[, 4:(ncol(pseudoCRD.splines.test) - 1)], scale)
  pseudoCRD.splines.train[, 4:(ncol(pseudoCRD.splines.train) - 1)] <-
    sapply(pseudoCRD.splines.train[, 4:(ncol(pseudoCRD.splines.train) - 1)], scale)
  
  # pseudoCRD.VEGsplines.test[, 4:(ncol(pseudoCRD.VEGsplines.test) - 1)] <-
  #   sapply(pseudoCRD.VEGsplines.test[, 4:(ncol(pseudoCRD.VEGsplines.test) - 1)], scale)
  # pseudoCRD.VEGsplines.train[, 4:(ncol(pseudoCRD.VEGsplines.train) - 1)] <-
  #   sapply(pseudoCRD.VEGsplines.train[, 4:(ncol(pseudoCRD.VEGsplines.train) - 1)], scale)
  
  # Recreating full datasets after scaling:
  pseudoCRD.nosplines <- rbind(pseudoCRD.nosplines.test, pseudoCRD.nosplines.train)
  pseudoCRD.splines <- rbind(pseudoCRD.splines.test, pseudoCRD.splines.train)
  # pseudoCRD.VEGsplines <- rbind(pseudoCRD.VEGsplines.test, pseudoCRD.VEGsplines.train)
  
  # Setting gy_adjusted (reps from pseudo-CRD) to NA for test genotypes:
  pseudoCRD.nosplines[which(pseudoCRD.nosplines$gid %in% unique(pseudoCRD.nosplines.test$gid)),
                      "gy_adjusted"] <- NA
  pseudoCRD.splines[which(pseudoCRD.splines$gid %in% unique(pseudoCRD.splines.test$gid)),
                      "gy_adjusted"] <- NA
  # pseudoCRD.VEGsplines[which(pseudoCRD.VEGsplines$gid %in% unique(pseudoCRD.VEGsplines.test$gid)),
  #                   "gy_adjusted"] <- NA
  
  # Saving the data (no time splines):
  datalist.nosplines <- list(data = as.data.frame(pseudoCRD.nosplines[, 3:(ncol(pseudoCRD.nosplines) - 1)]),
                             pred.target = unique(as.data.frame(pseudoCRD.nosplines.test[, c(3, ncol(pseudoCRD.nosplines.test))])),
                             test.set = unique(pseudoCRD.nosplines$gid[which(is.na(pseudoCRD.nosplines$gy_adjusted))]),
                             train.set = unique(pseudoCRD.nosplines$gid[which(!is.na(pseudoCRD.nosplines$gy_adjusted))]))
  
  names(datalist.nosplines$pred.target) <- c("G", "pred.target")
  names(datalist.nosplines$data)[1] <- "G"
  names(datalist.nosplines$data)[ncol(datalist.nosplines$data)] <- "Y"
  rlist::list.save(datalist.nosplines, file = sprintf("hyper_1415HEAT/datasets/nosplines/hyper_dataset_%d.RData", run))
  
  # Saving the data (time splines):
  datalist.splines <- list(data = as.data.frame(pseudoCRD.splines[, 3:(ncol(pseudoCRD.splines) - 1)]),
                           pred.target = unique(as.data.frame(pseudoCRD.splines.test[, c(3, ncol(pseudoCRD.splines.test))])),
                           test.set = unique(pseudoCRD.splines$gid[which(is.na(pseudoCRD.splines$gy_adjusted))]),
                           train.set = unique(pseudoCRD.splines$gid[which(!is.na(pseudoCRD.splines$gy_adjusted))]))
  
  names(datalist.splines$pred.target) <- c("G", "pred.target")
  names(datalist.splines$data)[1] <- "G"
  names(datalist.splines$data)[ncol(datalist.splines$data)] <- "Y"
  rlist::list.save(datalist.splines, file = sprintf("hyper_1415HEAT/datasets/splines/hyper_dataset_%d.RData", run))
  
  # Saving the data (VEG time splines):
  # datalist.VEGsplines <- list(data = as.data.frame(pseudoCRD.VEGsplines[, 3:(ncol(pseudoCRD.VEGsplines) - 1)]),
  #                             pred.target = unique(as.data.frame(pseudoCRD.VEGsplines.test[, c(3, ncol(pseudoCRD.VEGsplines.test))])),
  #                             test.set = unique(pseudoCRD.VEGsplines$gid[which(is.na(pseudoCRD.VEGsplines$gy_adjusted))]),
  #                             train.set = unique(pseudoCRD.VEGsplines$gid[which(!is.na(pseudoCRD.VEGsplines$gy_adjusted))]))
  # 
  # names(datalist.VEGsplines$pred.target) <- c("G", "pred.target")
  # names(datalist.VEGsplines$data)[1] <- "G"
  # names(datalist.VEGsplines$data)[ncol(datalist.VEGsplines$data)] <- "Y"
  # rlist::list.save(datalist.VEGsplines, file = sprintf("hyper_1415DRIP/datasets/VEGsplines/hyper_dataset_%d.RData", run))
}
toc()

