# K = 10 =======================================================================
rm(list = ls())
K <- 10

setwd("C:/Users/killi/Desktop/LINUX")

# Loading posterior loadings and predictions samples:
load(sprintf("megalmm_states/hyper/CV2_hyper_dataset_1_K%s/LAMBDA.RData", K))
load(sprintf("megalmm_states/hyper/CV2_hyper_dataset_1_K%s/LAMBDA_Y.RData", K))
load(sprintf("megalmm_states/hyper/CV2_hyper_dataset_1_K%s/PREDS.RData", K))

# Which factors are important for predicting yield?
(useful <- names(loadings.y.rounded)[which(abs(loadings.y.rounded) > 0.05)])

# 20000 x 10 x 63 array (samples x factors x features)
dim(Lambda_samples)
dimnames(Lambda_samples)
dimnames(Lambda_samples)[[2]] <- paste0("F", 1:dim(Lambda_samples)[2])

# 20000 x 1032 x 1 array (samples x genotypes x yield)
dim(pred_samples)
pred_samples <- pred_samples[, , 1] # Collapse into simple matrix
dimnames(pred_samples)[[2]] <- substr(dimnames(pred_samples)[[2]], 1, 7) # Remove weird part of G names

# Some quick examples of the posterior loading samples of different features on the first useful factor (F1):
plot(Lambda_samples[, "F1", "Y"])
plot(Lambda_samples[, "F1", "nm685_150310"])
plot(Lambda_samples[, "F1", "nm707_150310"])
plot(Lambda_samples[, "F1", "nm751_150310"])

# Some quick examples of the posterior predictions of different genotypes:
colnames(pred_samples)[1:3]
plot(pred_samples[, 1])
plot(pred_samples[, 2])
plot(pred_samples[, 3])

library(ggplot2)

# Subsetting to the factors useful for predicting yield, and yield itself:
## Yield =======================================================================
subset <- as.data.frame(Lambda_samples[, useful, "Y"])
subset$Iteration <- 1:nrow(subset)
subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Factor", values_to = "Loading")

# Plotting the significant Y loadings:
ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Factor)) +
  geom_point() +
  theme_dark() +
  ggtitle("Yield")
ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_YIELD.jpg", height = 20, width = 40, units = "cm")

# The important wavelengths around the "switching" point in gfBLUP:
(WL <- dimnames(Lambda_samples)[[3]][41:46])

# Plotting these wavelengths for one of the factors relevant to Y at a time:
## WLs on F1 ===================================================================
subset <- as.data.frame(Lambda_samples[, useful[1], WL])
subset$Iteration <- 1:nrow(subset)
subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")

ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
  geom_point() +
  theme_dark() +
  ggtitle("F1")
ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_WL_F1.jpg", height = 20, width = 40, units = "cm")

## WLs on F2 ===================================================================
subset <- as.data.frame(Lambda_samples[, useful[2], WL])
subset$Iteration <- 1:nrow(subset)
subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")

ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
  geom_point() +
  theme_dark() +
  ggtitle("F2")
ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_WL_F2.jpg", height = 20, width = 40, units = "cm")

## WLs on F3 ===================================================================
subset <- as.data.frame(Lambda_samples[, useful[3], WL])
subset$Iteration <- 1:nrow(subset)
subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")

ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
  geom_point() +
  theme_dark() +
  ggtitle("F3")
ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_WL_F3.jpg", height = 20, width = 40, units = "cm")

## WLs on F6 ===================================================================
subset <- as.data.frame(Lambda_samples[, useful[4], WL])
subset$Iteration <- 1:nrow(subset)
subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")

ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
  geom_point() +
  theme_dark() +
  ggtitle("F6")
ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_WL_F6.jpg", height = 20, width = 40, units = "cm")

## WLs on F9 ===================================================================
subset <- as.data.frame(Lambda_samples[, useful[5], WL])
subset$Iteration <- 1:nrow(subset)
subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")

ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
  geom_point() +
  theme_dark() +
  ggtitle("F9")
ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_WL_F9.jpg", height = 20, width = 40, units = "cm")

# # K = 50 =======================================================================
# rm(list = ls())
# K <- 50
# 
# setwd("C:/Users/killi/Desktop/LINUX")
# load(sprintf("megalmm_states/hyper/CV2_hyper_dataset_1_K%s/LAMBDA.RData", K))
# load(sprintf("megalmm_states/hyper/CV2_hyper_dataset_1_K%s/LAMBDA_Y.RData", K))
# load(sprintf("megalmm_states/hyper/CV2_hyper_dataset_1_K%s/PREDS.RData", K))
# 
# # Which factors are important for predicting yield?
# (useful <- names(loadings.y.rounded)[which(abs(loadings.y.rounded) > 0.05)])
# 
# # 20000 x 50 x 63 array (samples x factors x features)
# dim(Lambda_samples)
# dimnames(Lambda_samples)
# dimnames(Lambda_samples)[[2]] <- paste0("F", 1:dim(Lambda_samples)[2])
# 
# # 20000 x 1032 x 1 array (samples x genotypes x yield)
# dim(pred_samples)
# pred_samples <- pred_samples[, , 1] # Collapse into simple matrix
# dimnames(pred_samples)[[2]] <- substr(dimnames(pred_samples)[[2]], 1, 7) # Remove weird part of G names
# 
# # Some quick examples of the posterior loading samples of different features on the first useful factor (F1):
# plot(Lambda_samples[, "F1", "Y"])
# plot(Lambda_samples[, "F1", "nm685_150310"])
# plot(Lambda_samples[, "F1", "nm707_150310"])
# plot(Lambda_samples[, "F1", "nm751_150310"])
# 
# # Some quick examples of the posterior predictions of different genotypes:
# colnames(pred_samples)[1:3]
# 
# line <- lm(pred_samples[, 1] ~ matrix(1:nrow(pred_samples)))
# plot(pred_samples[, 1]); abline(line, col = "red")
# 
# line <- lm(pred_samples[, 2] ~ matrix(1:nrow(pred_samples)))
# plot(pred_samples[, 2]); abline(line, col = "red")
# 
# line <- lm(pred_samples[, 3] ~ matrix(1:nrow(pred_samples)))
# plot(pred_samples[, 3]); abline(line, col = "red")
# 
# library(ggplot2)
# 
# # Subsetting to the factors useful for predicting yield, and yield itself:
# ## Yield =======================================================================
# subset <- as.data.frame(Lambda_samples[, useful, "Y"])
# subset$Iteration <- 1:nrow(subset)
# subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Factor", values_to = "Loading")
# 
# # Plotting the significant Y loadings:
# ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Factor)) +
#   geom_point() +
#   theme_dark() +
#   ggtitle("Yield")
# ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K50_YIELD.jpg", height = 20, width = 40, units = "cm")
# 
# # The important wavelengths around the "switching" point in gfBLUP:
# (WL <- dimnames(Lambda_samples)[[3]][41:46])
# 
# # Plotting these wavelengths for one of the factors relevant to Y at a time:
# ## WLs on F1 ===================================================================
# subset <- as.data.frame(Lambda_samples[, useful[1], WL])
# subset$Iteration <- 1:nrow(subset)
# subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")
# 
# ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
#   geom_point() +
#   theme_dark() +
#   ggtitle("F1")
# ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K50_WL_F1.jpg", height = 20, width = 40, units = "cm")
# 
# ## WLs on F2 ===================================================================
# subset <- as.data.frame(Lambda_samples[, useful[2], WL])
# subset$Iteration <- 1:nrow(subset)
# subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")
# 
# ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
#   geom_point() +
#   theme_dark() +
#   ggtitle("F2")
# ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_WL_F2.jpg", height = 20, width = 40, units = "cm")
# 
# ## WLs on F3 ===================================================================
# subset <- as.data.frame(Lambda_samples[, useful[3], WL])
# subset$Iteration <- 1:nrow(subset)
# subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")
# 
# ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
#   geom_point() +
#   theme_dark() +
#   ggtitle("F3")
# ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_WL_F3.jpg", height = 20, width = 40, units = "cm")
# 
# ## WLs on F6 ===================================================================
# subset <- as.data.frame(Lambda_samples[, useful[4], WL])
# subset$Iteration <- 1:nrow(subset)
# subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")
# 
# ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
#   geom_point() +
#   theme_dark() +
#   ggtitle("F6")
# ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_WL_F6.jpg", height = 20, width = 40, units = "cm")
# 
# ## WLs on F9 ===================================================================
# subset <- as.data.frame(Lambda_samples[, useful[5], WL])
# subset$Iteration <- 1:nrow(subset)
# subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Wavelength", values_to = "Loading")
# 
# ggplot(data = subset, mapping = aes(x = Iteration, y = Loading, color = Wavelength)) +
#   geom_point() +
#   theme_dark() +
#   ggtitle("F9")
# ggsave(filename = "MegaLMM_ILLUSTRATION/MegaLMM_K10_WL_F9.jpg", height = 20, width = 40, units = "cm")




