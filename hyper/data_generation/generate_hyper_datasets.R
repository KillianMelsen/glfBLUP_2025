################################################################################
# Script that generates the 250 random train-test divisions of the
# hyperspectral data.
#
# - 26 training trials, 13 test trials

# Loading libraries:
library(gfBLUPold)
library(rlist)
library(vroom)
library(tidyverse)
library(stringr)
library(statgenGWAS)
library(lme4)
library(doParallel)
library(tictoc)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading functions:
for (file.name in list.files("helper_functions")) {
  source(paste0("helper_functions/", file.name))
}

# Setting seed for the main instance:
set.seed(1997)

# Number of datasets:
n.datasets <- 25

n.cores <- parallel::detectCores() - 2
work <- split(1:n.datasets, ceiling(seq_along(1:n.datasets) / ceiling(n.datasets / n.cores)))
seeds <- sample(1:1000, n.cores)
cl <- parallel::makeCluster(n.cores, outfile = "logs/hyper_data_generation.txt")
doParallel::registerDoParallel(cl)

tic("Dataset generation")
invisible(
foreach::foreach(i = 1:length(work), .packages = c("magrittr")) %dopar% {
  
  # Getting the dataset numbers that this worker client will generate:
  par.work <- work[[i]]
  
  # Setting the seed:
  set.seed(seeds[i])
  
  # Loading data:
  yield_1415 <- readxl::read_xlsx("hyper/data_generation/EYT_all_data_Y13_14_to_Y15_16.xlsx", sheet = 2)
  genotypes <- vroom::vroom("hyper/data_generation/Krause_et_al_2018_Genotypes.csv")
  pseudoCRD.unscaled <- vroom::vroom("hyper/data_generation/hyper_pseudoCRD_unscaled.csv")
  
  genotypes$GID <- as.character(genotypes$GID)
  pseudoCRD.unscaled$gid <- as.character(pseudoCRD.unscaled$gid)
  yield_1415$GID <- as.character(yield_1415$GID)
  
  # Store yield data from treatment B5IR:
  yield_1415_B5IR <- dplyr::filter(yield_1415, env == "B5IR")
  
  # Creating kinship matrix if required:
  if (!("K_hyper.RData" %in% list.files("genotypes"))) {
    # 1033 genotypes overlap between yield_1415_B5IR and genotypes file.
    K <- as.data.frame(genotypes[which(genotypes$GID %in% yield_1415_B5IR$GID), ])
    rownames(K) <- K$GID
    K <- K[colnames(K) != "GID"]
    K_gdata <- createGData(geno = K)
    K_gdata <- codeMarkers(K_gdata)
    
    # Note that the final kinship matrix has 1033 genotypes:
    K <- kinship(K_gdata$markers, "vanRaden")
    
    # Check for genotypes with missing yield (they will be removed):
    missing_yield <- as.character(yield_1415_B5IR[which(is.na(yield_1415_B5IR$gy)), "GID"])
    
    K <- K[setdiff(rownames(K), missing_yield), setdiff(colnames(K), missing_yield)]
    M_hyper <- K_gdata$markers[setdiff(rownames(K_gdata$markers), missing_yield),]
    
    save(K, "genotypes/K_hyper.RData")
    save(M_hyper, "genotypes/M_hyper.RData")
  } else {
    load("genotypes/K_hyper.RData")
    load("genotypes/M_hyper.RData")
  }
  
  for (run in par.work) {
    
    # Genotypes with missing yield:
    missing_yield <- yield_1415_B5IR[which(is.na(yield_1415_B5IR$gy)), 'GID']

    # Create training and test splits:
    # 39 trials.
    trials  <- unique(pseudoCRD.unscaled$trial)
    
    # 1/3 goes to test set, 2/3 to train set.
    test <- sample(trials, 13)
    train <- trials[!trials %in% test]
    
    # Subsetting all data:
    # Subsetting test data on secondary traits (1/3):
    test_data_sec <- pseudoCRD.unscaled[which(pseudoCRD.unscaled$trial %in% test),]
    
    # Subsetting test data on focal trait (1/3, contains 366 unique genotypes):
    # 1170 - (2 * 3 * 13) = 1092. 1092 / 13 / 3 = 28 other genotypes per trial.
    # So 28 * 13 = 364. Add two checks and we get 364 + 2 = 366 genotypes in total.
    test_data_foc <- yield_1415_B5IR[which(yield_1415_B5IR$trial.name %in% test),]
    
    # Subsetting training data on secondary traits (2/3):
    train_data_sec <- pseudoCRD.unscaled[which(pseudoCRD.unscaled$trial %in% train),]
    
    # Subsetting training data on focal trait (2/3, contains 730 unique genotypes):
    # Again, 2340 - (2 * 3 * 26) = 2184. 2148 / 26 / 3 = 28 other genotypes per trial.
    # So 28 * 26 = 728. Add two checks and we get 728 + 2 = 730 genotypes in total.
    train_data_foc <- yield_1415_B5IR[which(yield_1415_B5IR$trial.name %in% train),]
    
    # Generating adjusted yield BLUEs for the test set (removing design factors):
    lmm.yield <- lme4::lmer(gy ~ GID + (1|trial.no) + (1|trial.no:rep) + (1|trial.no:rep:subblock), 
                            control = lme4::lmerControl(calc.derivs = FALSE),
                            data = test_data_foc)
    
    # Storing yield BLUEs for the test set:
    BLUES_yield <- lme4::fixef(lmm.yield)
    
    # Setting correct names:
    names(BLUES_yield) <- levels(lmm.yield@frame[["GID"]])
    
    # Adding intercept to all but the first genotype:
    BLUES_yield <- c(BLUES_yield[1],(BLUES_yield[2:length(BLUES_yield)] + BLUES_yield[1]))
    
    # Creating BLUEs vector of the correct length (considering replicates):
    test_BLUES_yield <- BLUES_yield[match(lmm.yield@frame[["GID"]], names(BLUES_yield))]
    
    # Adding (possibly alpha-scaled) residuals to BLUEs for the test set:
    adjusted_BLUES <- BLUE_adjust(lmm.yield, "GID", alpha_scaling = FALSE)
    
    # Generate right sized dataframe:
    # Remove any rows with missing yield:
    gp_ready_foc_test <- dplyr::filter(test_data_foc, !is.na(gy))
    
    # Adjusted_BLUES and test_BLUES_yield are named. Is that necessary?
    gp_ready_foc_test$yield_adjusted_BLUES <- as.numeric(adjusted_BLUES)
    gp_ready_foc_test$BLUES_yield <- as.numeric(test_BLUES_yield)
    
    # Remove test genotypes with missing yield data in both yield and secondary dataframes:
    gp_ready_foc_test <- dplyr::filter(gp_ready_foc_test, !GID %in% missing_yield$GID)
    gp_ready_sec_test <- dplyr::filter(test_data_sec, !gid %in% missing_yield$GID)

    # Make sure the final dataframe has right size to merge with secondary traits:
    # Take the gp_read_foc_test dataframe, select the columns trial.name, GID, yield_adjusted_BLUES,
    # and BLUES_yield, then only keep the distinct rows. The result is written over gp_ready_foc_test:
    gp_ready_foc_test <- gp_ready_foc_test %>% 
      dplyr::select(trial.name, GID, yield_adjusted_BLUES:BLUES_yield) %>% 
      dplyr::distinct()
    
    # Finally bind the test dataframes containing (adjusted) BLUEs and secondary traits:
    #####
    # NOTE: Does bind_cols match dataframes based on order. As in, the rows (plants) need to
    # be in the same order in both dataframes? Otherwise secondary phenotypes and yields don't match.
    # That may not be an issue in CV1, but in CV2 it would be.
    #####
    # ANOTHER NOTE: length(unique(gp_ready_sec_test$gid)) = 366. Why not 27 genotypes * 13 trials + 2 checks
    # = 27*13+2 = 353? 
    # View(table(gp_ready_sec_test$gid)$Freq) shows 366 genotypes with two having a freq. of 39 (checks),
    # and 364 having a freq of 3. View(table(gp_ready_sec_test$trial)) shows each trial has a freq. of
    # 90. So 30 - 2 = 28 unique genotypes per trial...
    gp_ready_test <- gp_ready_sec_test %>% 
      dplyr::select(-trial) %>% 
      dplyr::bind_cols(., gp_ready_foc_test[, 3:4])
    
    # Doing similar things for the training data:
    lmm.yield <- lme4::lmer(gy ~ GID + (1|trial.no) + (1|trial.no:rep) + (1|trial.no:rep:subblock), 
                            control = lme4::lmerControl(calc.derivs = FALSE),
                            data = train_data_foc)
    
    BLUES_yield <- lme4::fixef(lmm.yield)
    names(BLUES_yield) <- levels(lmm.yield@frame[["GID"]])
    BLUES_yield <- c(BLUES_yield[1],(BLUES_yield[2:length(BLUES_yield)] + BLUES_yield[1]))
    train_BLUES_yield <- BLUES_yield[match(lmm.yield@frame[["GID"]], names(BLUES_yield))]
    
    # Adding (possibly alpha-scaled) residuals:
    ######
    # NOTE: why 2339 and not 728 * 3 + 2 * 26 * 3 = 2340?
    adjusted_BLUES <- BLUE_adjust(lmm.yield, "GID", alpha_scaling = FALSE)
    
    # Generate right sized dataframe:
    # Remove any rows with missing yield:
    gp_ready_foc_train <- dplyr::filter(train_data_foc, !is.na(gy))
    
    # Adjusted_BLUES and test_BLUES_yield are named. Is that necessary?
    gp_ready_foc_train$yield_adjusted_BLUES <- as.numeric(adjusted_BLUES)
    gp_ready_foc_train$BLUES_yield <- as.numeric(train_BLUES_yield)
    
    # remove genotypes with missing yield data
    # Finally end up with 2337 rows; 28 * 26 * 3 = 2184. Add 2 checks; 2184 + 2 * 26 * 3 = 2340.
    # Remove one genotype with one plot missing yield; 2340 - 1 * 3 = 2337.
    gp_ready_foc_train <- dplyr::filter(gp_ready_foc_train, !GID %in% missing_yield$GID)
    gp_ready_sec_train<- dplyr::filter(train_data_sec, !gid %in% missing_yield$GID)

    # make sure the final dataframe has right size to merge with secondary traits
    gp_ready_foc_train <- gp_ready_foc_train %>%
      dplyr::select(trial.name, GID, yield_adjusted_BLUES:BLUES_yield) %>%
      dplyr::distinct()
    
    gp_ready_train <- gp_ready_sec_train %>% 
      dplyr::select(-trial) %>% 
      dplyr::bind_cols(., gp_ready_foc_train[3:4])
    
    # Filter out any genotypes for which we do not have marker data:
    # We lost the 2 checks (in both sets) + 6 in the test set + 53 in the train set.
    # Total number of rows lost is 2 * 39 * 3 + 6 * 3 + 53 * 3 = 411.
    # 1170 + 2337 = 3507. 1074 + 2022 = 3096. 3507 - 411 is indeed 3096.
    gp_ready_train <- dplyr::filter(gp_ready_train, gid %in% genotypes$GID)
    gp_ready_test <- dplyr::filter(gp_ready_test, gid %in% genotypes$GID)
    
    # Scaling secondary trait data per column (time x wavelength combination):
    # Not scaling genotype BLUES:
    gp_ready_train[, 3:(length(gp_ready_train)-1)] <- sapply(gp_ready_train[, 3:(length(gp_ready_train)-1)], scale)
    gp_ready_test[, 3:(length(gp_ready_test)-1)] <- sapply(gp_ready_test[, 3:(length(gp_ready_test)-1)], scale)
    
    # create full dataset
    gp_full <- dplyr::bind_rows(gp_ready_test, gp_ready_train)
    
    # remove scaling attributes
    gp_full[,3:length(gp_full)] <- sapply(gp_full[,3:length(gp_full)], as.numeric)
    
    # Setting yield_adjusted_BLUES (reps from pseudo-CRD) to NA for test genotypes.
    gp_full[gp_full$gid %in% unique(gp_ready_test$gid), 623:length(gp_full)] <- NA
    
    # Saving the data:
    datalist <- list(data = as.data.frame(gp_full[,2:623]),
                     pred.target = unique(as.data.frame(gp_ready_test[, c(2, 624)])),
                     test.set = unique(gp_full$gid[which(is.na(gp_full$yield_adjusted_BLUES))]),
                     train.set = unique(gp_full$gid[which(!is.na(gp_full$yield_adjusted_BLUES))]))
    
    names(datalist$pred.target) <- c("G", "pred.target")
    names(datalist$data)[1] <- "G"
    names(datalist$data)[622] <- "Y"
    
    rlist::list.save(datalist, file = sprintf("hyper/datasets/hyper_dataset_%d.RData", run))
    
  }
})
doParallel::stopImplicitCluster()
parallel::stopCluster(cl)
toc()


