# Fitting the models and getting the corrected data for a single feature (wl-date combination)
# takes ~ 9 seconds, so should be ~ 90 mins for all features.
# Libraries:
library(tidyverse)
library(LMMsolver)
library(tictoc)
library(statgenHTP)
set.seed(1997)

# Subsetting to Y1415 and Bed5IR irrigation treatment:
df <- read_csv("hyper_datafiles/Y1415_Hyperspectral_Raw_24Mar2022.csv", na = "*")
df <- df %>% dplyr::filter(year == "Y1415" & conditions == "Heat")
df <- df %>% dplyr::filter(date %in% c(150414, 150424, 150428, 150506)) # The last data in May is not mentioned in the paper about this data...
df$gid <- as.factor(df$gid)

# Discarding varieties with missing measurements:
df <- droplevels(na.omit(df))
temp <- as.data.frame(table(df$gid))
discard <- as.character(temp$Var1[which(temp$Freq < 12)])
df <- droplevels(df[which(!(as.character(df$gid) %in% discard)),])

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
    
    cat(sprintf("Working on feature %d / 62 for date %d / %d\n\n", j, i, length(dates)))
    
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
    saveRDS(coords, "hyper_1415HEAT/data_generation/coords.rds")
  }
  i <- i + 1
}
toc()
saveRDS(corr.list, "hyper_1415HEAT/data_generation/corr_list.rds")



rm(list = setdiff(ls(), c("dates")))
corr.list <- readRDS("hyper_1415HEAT/data_generation/corr_list.rds")

# Putting all dataframes below each other:
data <- corr.list[[as.character(dates[1])]]
for (i in 2:length(dates)) {
  data <- rbind(data, corr.list[[as.character(dates[i])]])
}

# Saving data without fitting splines through the 4 timepoints for CV2VEG:
data.wide <- as.data.frame(tidyr::pivot_wider(data,
                                              names_from = date,
                                              names_glue = "{.value}_{date}",
                                              values_from = 11:72))
saveRDS(data.wide, "hyper_1415HEAT/data_generation/hyper_pseudoCRD_nosplines.rds")

# Fitting splines through 4 timepoints and taking fitted values -----
# Making an empty copy for storing the data:
rm(list = setdiff(ls(), c("dates")))
corr.list <- readRDS("hyper_1415HEAT/data_generation/corr_list.rds")

# Putting all dataframes below each other:
data <- corr.list[[as.character(dates[1])]]
for (i in 2:length(dates)) {
  data <- rbind(data, corr.list[[as.character(dates[i])]])
}
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
                        # newtimes = seq(min(fit.psHDM$time[["timeNumber"]]),
                        #                max(fit.psHDM$time[["timeNumber"]]),
                        #                length.out = 100),
                        pred = list(pop = TRUE, geno = TRUE, plot = TRUE),
                        se = list(pop = FALSE, geno = FALSE, plot = FALSE),
                        trace = FALSE)

  # ## Population-specific trajectories.
  # plot(pred.psHDM, plotType = "popTra", themeSizeHDM = 10)
  #
  # ## Population and genotype-specific trajectories.
  # plot(pred.psHDM, plotType = "popGenoTra", themeSizeHDM = 10)
  #
  # ## As an example we used ten randomly selected genotypes
  # set.seed(1)
  # plot.genos  <- sample(pred.psHDM$genoLevel$genotype,10, replace = FALSE)
  # names.genos <- as.character(plot.genos)
  #
  # ## Genotype- and plot-specific trajectories.
  # plot(pred.psHDM,
  #      plotType = "genoPlotTra",
  #      genotypes = plot.genos, genotypeNames = names.genos,
  #      themeSizeHDM = 10)

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
saveRDS(data.splines.wide, "hyper_1415HEAT/data_generation/hyper_pseudoCRD_splines.rds")

# Pre-processing yield data: ===================================================
rm(list = ls())
set.seed(1997)

# Loading data:
yield_1415 <- readxl::read_xlsx("hyper_datafiles/EYT_all_data_Y13_14_to_Y15_16.xlsx", sheet = 2)
genotypes <- vroom::vroom("hyper_datafiles/Krause_et_al_2018_Genotypes.csv")
pseudoCRD.nosplines <- readRDS("hyper_1415HEAT/data_generation/hyper_pseudoCRD_nosplines.rds")
pseudoCRD.splines <- readRDS("hyper_1415HEAT/data_generation/hyper_pseudoCRD_splines.rds")
coords <- readRDS("hyper_1415HEAT/data_generation/coords.rds")

genotypes$GID <- as.character(genotypes$GID)
pseudoCRD.nosplines$gid <- as.character(pseudoCRD.nosplines$gid)
pseudoCRD.splines$gid <- as.character(pseudoCRD.splines$gid)
yield_1415$GID <- as.character(yield_1415$GID)

# Store yield data from treatment HEAT:
yield_1415_B5IR <- dplyr::filter(yield_1415, env == "BLHT")

# Assuming they have already been made in the B5IR data generation script:
load("genotypes/K_hyper.RData")
load("genotypes/M_hyper.RData")

# Remove plots with missing yield, not whole sets of genotypes cause we would discard
# all checks:
missing_yield <- yield_1415_B5IR[which(is.na(yield_1415_B5IR$gy)), c('plot', "GID")]
yield_1415_B5IR <- droplevels(yield_1415_B5IR[which(!(yield_1415_B5IR$plot %in% missing_yield$plot)),])

# Removing genotypes from yield df that we did not save coords for because they had missing
# hyperspectral data:
keep <- unique(as.character(pseudoCRD.nosplines$gid))
yield_1415_B5IR <- droplevels(yield_1415_B5IR[which(yield_1415_B5IR$GID %in% keep),])

# Removing genotypes that we do not have hyperspectral data or yield for:
# keep <- intersect(unique(as.character(pseudoCRD.nosplines$gid)), unique(as.character(yield_1415_B5IR$GID)))
# yield_1415_B5IR <- droplevels(yield_1415_B5IR[which(yield_1415_B5IR$GID %in% keep),])

# Focal trait dataframes have row and col numbers within trials, so we replace
# them with overall row/col numbers using the plot numbers and hyper data:
yield_1415_B5IR[, c("row", "column")] <- coords[match(yield_1415_B5IR$plot, coords$plot), c("row", "col")]

# Making sure the design factors are actually factors and not numerical values:
yield_1415_B5IR$trial.name <- as.factor(yield_1415_B5IR$trial.name)
yield_1415_B5IR$rep <- as.factor(yield_1415_B5IR$rep)
yield_1415_B5IR$subblock <- as.factor(yield_1415_B5IR$subblock)
yield_1415_B5IR$R <- as.factor(yield_1415_B5IR$row)
yield_1415_B5IR$C <- as.factor(yield_1415_B5IR$column)

# Generating adjusted yield BLUEs for the test set (removing design factors):
# ggplot(yield_1415_B5IR, aes(x = column, y = row, fill = !! rlang::sym("gy"))) +
#   geom_tile(show.legend = TRUE) +
#   scale_fill_gradientn(colours = topo.colors(100))+
#   labs(title = "Raw yield") +
#   coord_fixed() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())

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
pseudoCRD.nosplines <- dplyr::filter(pseudoCRD.nosplines, !gid %in% missing_yield$GID)
pseudoCRD.splines <- dplyr::filter(pseudoCRD.splines, !gid %in% missing_yield$GID)

# Select relevant columns:
gp_ready_foc <- gp_ready_foc %>% 
  dplyr::select(trial.name, plot, GID, gy_adjusted, target) %>% 
  dplyr::distinct()

# Finally bind the dataframes containing adjusted grain yield and secondary traits:
# Making sure that the secondary and focal trait data of each plot gets matched correctly:
pseudoCRD.nosplines <- pseudoCRD.nosplines[match(gp_ready_foc$plot, pseudoCRD.nosplines$plot),]
pseudoCRD.splines <- pseudoCRD.splines[match(gp_ready_foc$plot, pseudoCRD.splines$plot),]

# Final dataframe will have 9 design columns, 62*t sec columns, and focal trait + target columns:
gp_ready.nosplines <- pseudoCRD.nosplines %>%
  dplyr::bind_cols(., gp_ready_foc[, 4:5])

gp_ready.splines <- pseudoCRD.splines %>% 
  dplyr::bind_cols(., gp_ready_foc[, 4:5])

# Filter out any genotypes for which we do not have marker data:
# We lost the 2 checks (in both sets) + 6 in the test set + 53 in the train set.
# Total number of rows lost is 2 * 39 * 3 + 6 * 3 + 53 * 3 = 411.
# 1170 + 2337 = 3507. 1074 + 2022 = 3096. 3507 - 411 is indeed 3096.
gp_ready.nosplines <- dplyr::filter(gp_ready.nosplines, gid %in% genotypes$GID)
gp_ready.splines <- dplyr::filter(gp_ready.splines, gid %in% genotypes$GID)

gp_ready.nosplines <- dplyr::filter(gp_ready.nosplines, gid %in% rownames(K))
gp_ready.splines <- dplyr::filter(gp_ready.splines, gid %in% rownames(K))

gp_ready.nosplines <- gp_ready.nosplines[, c(4:5, 8, 10:ncol(gp_ready.nosplines))]
gp_ready.splines <- gp_ready.splines[, c(4:5, 8, 10:ncol(gp_ready.splines))]

saveRDS(gp_ready.nosplines, "hyper_1415HEAT/data_generation/pseudoCRD_nosplines.rds")
saveRDS(gp_ready.splines, "hyper_1415HEAT/data_generation/pseudoCRD_splines.rds")






