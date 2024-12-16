library(LMMsolver)

# Reading raw data:
data <- read.csv("hyper/data_generation/Y1415_Y1516_Y1617_Hyperspectral_Raw_24Mar2022.csv")

# Subsetting to 14-15, removing year column, and subsetting to Bed5IR and removing
# conditions column:
data <- droplevels(data[which(data$year == "Y1415"),])[, -1]
data <- droplevels(data[which(data$conditions == "Bed5IR"),])[, -9]

# Converting hyperspectral measurements to numeric data:
data[, 9:70] <- sapply(data[, 9:70], as.numeric)

# Changing column names for the wavelengths:
colnames(data)[9:70] <- gsub("X(...)(nm)", "\\2\\1", colnames(data)[9:70])

# Widening the dataframe:
data <- as.data.frame(tidyr::pivot_wider(data,
                                         names_from = date,
                                         names_glue = "{.value}_{date}",
                                         values_from = 9:70))

# making trial names consistent with yield data:
data$trial <- gsub(" ", "", data$trial)






