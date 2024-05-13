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
