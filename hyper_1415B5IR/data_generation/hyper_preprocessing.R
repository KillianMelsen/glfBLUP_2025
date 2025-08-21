# This script does the pre-processing of the raw hyperspectral B5IR data.
# It produces a number of .rds files stored in `hyper_1415B5IR/data_generation/`
# that are subsequently used to generate the datafiles for the 250 replications.

# Libraries:
library(tidyverse)
library(LMMsolver)
library(tictoc)
library(statgenHTP)
set.seed(1997)

# Subsetting to Y1415 and B5IR irrigation treatment:
df <- read_csv("hyper_datafiles/Y1415_Hyperspectral_Raw_24Mar2022.csv", na = "*")
df <- df %>% dplyr::filter(year == "Y1415" & conditions == "Bed5IR")
df$gid <- as.factor(df$gid)

# Change the names of the wavelengths so we have a letter first:
colnames(df)[11:72] <- gsub("(...)(nm)", "\\2\\1", colnames(df)[11:72])

# Making trial names consistent with yield data:
df$trial <- gsub(" ", "", df$trial)

# Setting a fixed plot order so we can make sure all corrected data is combined
# in the right way:
plot.order <- unique(df$plot)

# Loop through dates:
dates <- unique(df$date)
corr.list <- vector("list", length(dates))
names(corr.list) <- dates
d <- dates[1]
i <- 1
tic()
for (d in dates) {
  df.sel <- df %>% dplyr::filter(date == d)
  
  # Making a copy to store the corrected data:
  df.corr <- df.sel
  df.corr[, 11:ncol(df.corr)] <- NA
  
  df.sel$R <- as.factor(df.sel$row)
  df.sel$C <- as.factor(df.sel$col)
  
  # Just to be sure: match subsetted dataframes to plot order:
  df.sel <- df.sel[match(plot.order, df.sel$plot),]
  df.corr <- df.corr[match(plot.order, df.corr$plot),]
  
  # Loop through the wavelengths within each date:
  features <- colnames(df.sel)[11:(ncol(df.sel) - 2)]
  f <- features[1]
  j <- 1
  for (f in features) {
    
    cat(sprintf("Working on feature %d / 62 for date %d / 10\n\n", j, i))
    
    # Fixed part includes intercept and genotypes:
    fixed <- as.formula(paste0(f, " ~ 1 + gid"))
    fit <- LMMsolve(fixed = fixed,
                    random = ~ trial + R + C,
                    spline = ~ spl2D(row, col, nseg = c(25, 70)),
                    data = df.sel,
                    maxit = 1000)
    
    # Getting estimates:
    estimates <- coef(fit)
    
    # Intercept:
    int <- as.numeric(estimates$`(Intercept)`)
    
    # Genotype BLUEs:
    BLUEs <- as.numeric(estimates$gid)
    names(BLUEs) <- names(estimates$gid)
    names(BLUEs) <- gsub("gid_(.*)", "\\1", names(BLUEs))
    
    # Residuals:
    resids <- as.numeric(fit$residuals)
    
    # Getting the order of genotypes in the data vector and expanding the BLUEs:
    order <- as.character(df.sel$gid)
    BLUEs <- BLUEs[match(order, names(BLUEs))]
    
    # Computing the corrected data:
    corr <- as.numeric(rep(int, nrow(df.sel)) + BLUEs + resids)
    
    # Putting the corrected data into df.corr:
    df.corr[, f] <- corr
    j <- j + 1
  }
  corr.list[[as.character(d)]] <- df.corr
  
  # Saving a list of plot numbers and corresponding row and column numbers cause
  # we need those later:
  if (i == 1) {
    coords <- df.sel[, c("plot", "row", "col")]
    saveRDS(coords, "hyper_1415B5IR/data_generation/coords.rds")
  }
  i <- i + 1
}
toc()
saveRDS(corr.list, "hyper_1415B5IR/data_generation/corr_list.rds")

rm(list = setdiff(ls(), c("dates")))
corr.list <- readRDS("hyper_1415B5IR/data_generation/corr_list.rds")

# Putting all dataframes below each other:
data <- corr.list[[as.character(dates[1])]]
for (i in 2:length(dates)) {
  data <- rbind(data, corr.list[[as.character(dates[i])]])
}

# Fitting splines through 10 timepoints and taking fitted values -----
data.splines <- data
data.splines[, 11:72] <- NA

# Converting the data for use with statgenHTP:
data$date <- as.Date(as.character(data$date), format = "%y%m%d")
data <- as.data.frame(data)
data$plot <- as.factor(data$plot)
names(data)[which(names(data) == "date")] <- "timePoint"
data$timePoint <- as.POSIXct(data$timePoint)
data$pop <- "pop"; data$pop <- as.factor(data$pop)
names(data)[which(names(data) == "row")] <- "rowId"
names(data)[which(names(data) == "col")] <- "colId"

features <- colnames(data)[11:72]
f <- features[1]
i <- 1
for (f in features) {
  cat(sprintf("Working on feature %d / 62\n\n", i))
  # Fit P-Splines Hierarchical Curve Data Model for all genotypes:
  fit.psHDM  <- fitSplineHDM(inDat = data,
                             trait = f,
                             pop = "pop",
                             genotype = "gid",
                             plotId = "plot",
                             difVar = list(geno = FALSE, plot = FALSE),
                             smoothPop = list(nseg = 5, bdeg = 2, pord = 2),
                             smoothGeno = list(nseg = 5, bdeg = 3, pord = 2),
                             smoothPlot = list(nseg = 7, bdeg = 3, pord = 2),
                             trace = FALSE)

  # Predict the P-Splines Hierarchical Curve Data Model on a dense grid:
  pred.psHDM <- predict(object = fit.psHDM,
                        pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                        se = list(pop = FALSE, geno = FALSE, plot = FALSE),
                        trace = FALSE)

  # Making sure it's all in the same order:
  order <- paste(as.character(data$timePoint), as.character(data$plot), sep = "_")

  output <- pred.psHDM$plotLevel[, c("timePoint", "plotId", "genotype", "fPlot", "obsPlot")]
  output$order <- paste(as.character(output$timePoint), as.character(output$plotId), sep = "_")
  output <- output[match(order, output$order),]
  data.splines[, f] <- as.numeric(output$fPlot)
  i <- i + 1
}

data.splines.wide <- as.data.frame(tidyr::pivot_wider(data.splines,
                                                      names_from = date,
                                                      names_glue = "{.value}_{date}",
                                                      values_from = 11:72))
saveRDS(data.splines.wide, "hyper_1415B5IR/data_generation/hyper_pseudoCRD_splines.rds")

# Fitting splines through the vegetative timepoints only and taking fitted values -----
# Putting all dataframes below each other:
rm(list = setdiff(ls(), c("dates")))
corr.list <- readRDS("hyper_1415B5IR/data_generation/corr_list.rds")

# Putting all dataframes below each other:
data <- corr.list[[as.character(dates[1])]]
for (i in 2:length(dates)) {
  data <- rbind(data, corr.list[[as.character(dates[i])]])
}
data <- droplevels(data[which(data$date %in% c("150110", "150119", "150204", "150209")),])

# Making an empty copy for storing the data:
data.splines <- data
data.splines[, 11:72] <- NA

# Converting the data for use with statgenHTP:
data$date <- as.Date(as.character(data$date), format = "%y%m%d")
data <- as.data.frame(data)
data$plot <- as.factor(data$plot)
names(data)[which(names(data) == "date")] <- "timePoint"
data$timePoint <- as.POSIXct(data$timePoint)
data$pop <- "pop"; data$pop <- as.factor(data$pop)
names(data)[which(names(data) == "row")] <- "rowId"
names(data)[which(names(data) == "col")] <- "colId"

features <- colnames(data)[11:72]
f <- features[1]
i <- 1
for (f in features) {
  cat(sprintf("Working on feature %d / 62\n\n", i))
  # Fit P-Splines Hierarchical Curve Data Model for all genotypes:
  fit.psHDM  <- fitSplineHDM(inDat = data,
                             trait = f,
                             pop = "pop",
                             genotype = "gid",
                             plotId = "plot",
                             difVar = list(geno = FALSE, plot = FALSE),
                             smoothPop = list(nseg = 5, bdeg = 2, pord = 2),
                             smoothGeno = list(nseg = 5, bdeg = 3, pord = 2),
                             smoothPlot = list(nseg = 7, bdeg = 3, pord = 2),
                             trace = FALSE)
  
  # Predict the P-Splines Hierarchical Curve Data Model on a dense grid:
  pred.psHDM <- predict(object = fit.psHDM,
                        pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                        se = list(pop = FALSE, geno = FALSE, plot = FALSE),
                        trace = FALSE)
  
  # Making sure it's all in the same order:
  order <- paste(as.character(data$timePoint), as.character(data$plot), sep = "_")
  
  output <- pred.psHDM$plotLevel[, c("timePoint", "plotId", "genotype", "fPlot", "obsPlot")]
  output$order <- paste(as.character(output$timePoint), as.character(output$plotId), sep = "_")
  output <- output[match(order, output$order),]
  data.splines[, f] <- as.numeric(output$fPlot)
  i <- i + 1
}

data.splines.wide <- as.data.frame(tidyr::pivot_wider(data.splines,
                                                      names_from = date,
                                                      names_glue = "{.value}_{date}",
                                                      values_from = 11:72))
saveRDS(data.splines.wide, "hyper_1415B5IR/data_generation/hyper_pseudoCRD_VEGsplines.rds")

# Pre-processing yield data: ===================================================
rm(list = ls())
library(statgenGWAS)
set.seed(1997)

# Loading data:
yield_1415 <- readxl::read_xlsx("hyper_datafiles/EYT_all_data_Y13_14_to_Y15_16.xlsx", sheet = 2)
genotypes <- vroom::vroom("hyper_datafiles/Krause_et_al_2018_Genotypes.csv")

pseudoCRD.splines <- readRDS("hyper_1415B5IR/data_generation/hyper_pseudoCRD_splines.rds")
pseudoCRD.VEGsplines <- readRDS("hyper_1415B5IR/data_generation/hyper_pseudoCRD_VEGsplines.rds")
coords <- readRDS("hyper_1415B5IR/data_generation/coords.rds")

genotypes$GID <- as.character(genotypes$GID)
pseudoCRD.splines$gid <- as.character(pseudoCRD.splines$gid)
pseudoCRD.VEGsplines$gid <- as.character(pseudoCRD.VEGsplines$gid)
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
  M_hyper <- K_gdata$markers
  
  save(K, file = "genotypes/K_hyper.RData")
  save(M_hyper, file = "genotypes/M_hyper.RData")
} else {
  load("genotypes/K_hyper.RData")
  load("genotypes/M_hyper.RData")
}

# Remove genotypes with missing yield (only a single variety):
missing_yield <- yield_1415_B5IR[which(is.na(yield_1415_B5IR$gy)), 'GID']
yield_1415_B5IR <- droplevels(yield_1415_B5IR[which(!(yield_1415_B5IR$GID %in% missing_yield$GID)),])

# Focal trait dataframes have row and col numbers within trials, so we replace
# them with overall row/col numbers using the plot numbers and hyper data:
yield_1415_B5IR[, c("row", "column")] <- coords[match(yield_1415_B5IR$plot, coords$plot), c("row", "col")]

# Making sure the design factors are actually factors and not numerical values:
yield_1415_B5IR$trial.name <- as.factor(yield_1415_B5IR$trial.name)
yield_1415_B5IR$rep <- as.factor(yield_1415_B5IR$rep)
yield_1415_B5IR$subblock <- as.factor(yield_1415_B5IR$subblock)
yield_1415_B5IR$R <- as.factor(yield_1415_B5IR$row)
yield_1415_B5IR$C <- as.factor(yield_1415_B5IR$column)

lmm.yield <- LMMsolve(fixed = gy ~ GID,
                      random = ~ trial.name + R + C,
                      spline = ~ spl2D(row, column, nseg = c(25, 70)),
                      data = yield_1415_B5IR,
                      maxit = 1000)

# Getting estimates:
estimates <- coef(lmm.yield)

# Intercept:
int <- as.numeric(estimates$`(Intercept)`)

# Genotype BLUEs:
BLUEs <- as.numeric(estimates$GID)
names(BLUEs) <- names(estimates$GID)
names(BLUEs) <- gsub("GID_(.*)", "\\1", names(BLUEs))

# Residuals:
resids <- as.numeric(lmm.yield$residuals)

# Getting the order of genotypes in the data vector and expanding the BLUEs:
order <- as.character(yield_1415_B5IR$GID)
BLUEs <- BLUEs[match(order, names(BLUEs))]

# Computing the corrected data:
corr <- rep(int, nrow(yield_1415_B5IR)) + BLUEs + resids

# Generate right sized dataframe:
# Remove any rows with missing yield (left it here, but should just copy as
# we already removed the one genotype with missing yield):
gp_ready_foc <- dplyr::filter(yield_1415_B5IR, !is.na(gy))
gp_ready_foc$gy_adjusted <- as.numeric(corr)
gp_ready_foc$target <- as.numeric(BLUEs)

# Remove genotypes with missing yield data in both yield and secondary dataframes:
# This should also not do anything...
gp_ready_foc <- dplyr::filter(gp_ready_foc, !GID %in% missing_yield$GID)
# This one will do something (remove three rows):
pseudoCRD.splines <- dplyr::filter(pseudoCRD.splines, !gid %in% missing_yield$GID)
pseudoCRD.VEGsplines <- dplyr::filter(pseudoCRD.VEGsplines, !gid %in% missing_yield$GID)

# Select relevant columns:
gp_ready_foc <- gp_ready_foc %>% 
  dplyr::select(trial.name, plot, GID, gy_adjusted, target) %>% 
  dplyr::distinct()

# Finally bind the dataframes containing adjusted grain yield and secondary traits:
# Making sure that the secondary and focal trait data of each plot gets matched correctly:
pseudoCRD.splines <- pseudoCRD.splines[match(gp_ready_foc$plot, pseudoCRD.splines$plot),]
pseudoCRD.VEGsplines <- pseudoCRD.VEGsplines[match(gp_ready_foc$plot, pseudoCRD.VEGsplines$plot),]

# Final dataframe will have 9 design columns, 620 sec columns, and focal trait + target columns:
gp_ready.splines <- pseudoCRD.splines %>% 
  dplyr::bind_cols(., gp_ready_foc[, 4:5])

gp_ready.VEGsplines <- pseudoCRD.VEGsplines %>% 
  dplyr::bind_cols(., gp_ready_foc[, 4:5])

# Filter out any genotypes for which we do not have marker data:
# We lost the 2 checks (in both sets) + 6 in the test set + 53 in the train set.
# Total number of rows lost is 2 * 39 * 3 + 6 * 3 + 53 * 3 = 411.
# 1170 + 2337 = 3507. 1074 + 2022 = 3096. 3507 - 411 is indeed 3096.
gp_ready.splines <- dplyr::filter(gp_ready.splines, gid %in% genotypes$GID)
gp_ready.VEGsplines <- dplyr::filter(gp_ready.VEGsplines, gid %in% genotypes$GID)

gp_ready.splines <- gp_ready.splines[, c(4:5, 8, 10:631)]
gp_ready.VEGsplines <- gp_ready.VEGsplines[, c(4:5, 8, 10:259)]

saveRDS(gp_ready.splines, "hyper_1415B5IR/data_generation/pseudoCRD_splines.rds")
saveRDS(gp_ready.VEGsplines, "hyper_1415B5IR/data_generation/pseudoCRD_VEGsplines.rds")






