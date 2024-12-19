# Fitting the models and getting the corrected data for a single feature (wl-date combination)
# takes ~ 9 seconds, so should be ~ 90 mins for all features.
# Libraries:
library(tidyverse)
library(LMMsolver)
library(tictoc)
library(statgenHTP)
set.seed(1997)

# Subsetting to Y1415 and Bed5IR irrigation treatment:
df <- read_csv("hyper/data_generation/Y1415_Y1516_Y1617_Hyperspectral_Raw_24Mar2022.csv", na = "*")
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
    
    # newdata <- data.frame(gid = unique(df.sel$gid),
    #                       row = 0,
    #                       col = 0)
    # preds <- predict(fit, newdata = newdata)
    
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
    saveRDS(coords, "hyper/data_generation/coords.rds")
  }
  i <- i + 1
}
toc()
saveRDS(corr.list, "hyper/data_generation/corr_list.rds")
corr.list <- readRDS("hyper/data_generation/corr_list.rds")

# Putting all dataframes below each other:
data <- corr.list[[as.character(dates[1])]]
for (i in 2:length(dates)) {
  data <- rbind(data, corr.list[[as.character(dates[i])]])
}

# Option 1: not fitting splines through the 10 timepoints ----------------------
data.wide <- as.data.frame(tidyr::pivot_wider(data,
                                              names_from = date,
                                              names_glue = "{.value}_{date}",
                                              values_from = 11:72))
saveRDS(data.wide, "hyper/data_generation/hyper_pseudoCRD_nosplines.rds")

# Option 2: fitting splines through 10 timepoints and taking fitted values -----
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
                        newtimes = seq(min(fit.psHDM$time[["timeNumber"]]),
                                       max(fit.psHDM$time[["timeNumber"]]),
                                       length.out = 100),
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
saveRDS(data.splines.wide, "hyper/data_generation/hyper_pseudoCRD_splines.rds")










