library(rlist)
set.seed(1997)

# Load the first hyperspectral dataset:
data <- list.load("hyper/datasets/hyper_dataset_1.RData")$data

# Remove yield:
data <- data[, -ncol(data)]

# Convert to means (= BLUEs):
means <- gfBLUP:::genotypeMeans(data)

# Remove the G column:
rownames(means) <- means$G
means <- means[, -1]

# Mean-centering and scaling:
means.scaled <- scale(means)

# Calculating H matrix:
H <- (means.scaled %*% t(means.scaled)) / 620

# Saving:
save(H, file = "genotypes/H_hyper.RData")


# Now doing the same while only using the hyperspectral measurements from the
# vegetative phase:
# 9 feb is last day of VEG, 25 feb is heading, 10 march is start of grain filling:
dates <- c("150110", "150119", "150204", "150209")
select <- which(substr(colnames(means.scaled), 7, 12) %in% dates)
means.scaled <- means.scaled[, select]

# Calculating H matrix:
H <- (means.scaled %*% t(means.scaled)) / 248

# Saving:
save(H, file = "genotypes/H_hyper_VEG.RData")
