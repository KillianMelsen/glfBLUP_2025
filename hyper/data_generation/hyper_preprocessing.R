# Fitting the models and getting the corrected data for a single feature (wl-date combination)
# takes ~ 9 seconds, so should be ~ 90 mins for all features.
# Libraries:
library(tidyverse)
library(LMMsolver)
library(tictoc)

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





