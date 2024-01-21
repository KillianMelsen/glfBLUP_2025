# Numbers of MegaLMM factors:
# Ms <- c(5, 10)
Ms <- c(50)

for (M in Ms) {
  
  # Loading posterior samples for loadings and predictions, and the mean posterior yield loadings.
  load(sprintf("posterior_arrays/LAMBDA_M%s.RData", M))
  load(sprintf("mean_posterior_loadings_Y/LAMBDA_Y_M%s.RData", M))
  load(sprintf("posterior_arrays/PREDS_M%s.RData", M))
  
  # Which factors are important for predicting yield?
  (useful <- names(loadings.y)[which(abs(loadings.y) > 0.05)])
  
  # 25000 x 50 x 63 array (samples x factors x features)
  dim(Lambda_samples)
  dimnames(Lambda_samples)
  dimnames(Lambda_samples)[[2]] <- paste0("F", 1:dim(Lambda_samples)[2])
  
  # 25000 x 1032 x 1 array (samples x genotypes x yield)
  dim(pred_samples)
  pred_samples <- pred_samples[, , 1] # Collapse into simple matrix
  dimnames(pred_samples)[[2]] <- substr(dimnames(pred_samples)[[2]], 1, 7) # Remove weird part of G names
  
  # Some quick examples of the posterior loading samples of different features on the first useful factor (F1):
  plot(Lambda_samples[, useful[1], "Y"])
  plot(Lambda_samples[, useful[1], "nm685_150310"])
  plot(Lambda_samples[, useful[1], "nm707_150310"])
  plot(Lambda_samples[, useful[1], "nm751_150310"])
  
  # Some quick examples of the posterior predictions of different genotypes:
  colnames(pred_samples)[1:3]
  
  line <- lm(pred_samples[, 1] ~ matrix(1:nrow(pred_samples)))
  plot(pred_samples[, 1]); abline(line, col = "red")
  
  line <- lm(pred_samples[, 2] ~ matrix(1:nrow(pred_samples)))
  plot(pred_samples[, 2]); abline(line, col = "red")
  
  line <- lm(pred_samples[, 3] ~ matrix(1:nrow(pred_samples)))
  plot(pred_samples[, 3]); abline(line, col = "red")
  
  library(ggplot2)
  
  # Subsetting to the factors useful for predicting yield, and yield itself:
  ## Yield =======================================================================
  subset <- as.data.frame(Lambda_samples[, useful, "Y"])
  subset$Iteration <- 1:nrow(subset)
  subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1), names_to = "Factor", values_to = "Loading")
  subset$Factor <- factor(subset$Factor, levels = useful)
  
  pal <- RColorBrewer::brewer.pal(12, "Paired")
  
  # Shuffling colors:
  set.seed(1997)
  cols <- sample(pal[1:length(useful)], length(pal[1:length(useful)]))
  
  # Plotting the significant Y loadings:
  ggplot(data = subset, mapping = aes(x = Iteration, y = Loading)) +
    geom_point(mapping = aes(fill = Factor), color = "black", pch = 21, stroke = 0.1) +
    scale_fill_manual(values = pal[1:length(useful)]) +
    theme_classic() +
    ggtitle("Yield")
  
  ggsave(filename = sprintf("mean_posterior_loadings_Y/M%s_Y.jpg", M), height = 20, width = 40, units = "cm")
  
  # The important wavelengths around the "switching" point in gfBLUP:
  (WL <- dimnames(Lambda_samples)[[3]][41:45])
  
  # Shuffling colors:
  set.seed(1997)
  cols <- sample(pal[1:length(WL)], length(pal[1:length(WL)]))
  
  for (f in useful) {
    
    subset <- as.data.frame(Lambda_samples[, f, WL])
    subset$Iteration <- 1:nrow(subset)
    subset <- tidyr::pivot_longer(subset, cols = 1:(ncol(subset) - 1),
                                  names_to = "Wavelength", values_to = "Loading")
    
    ggplot(data = subset, mapping = aes(x = Iteration, y = Loading)) +
      geom_point(mapping = aes(fill = Wavelength), color = "black", pch = 21, stroke = 0.1) +
      scale_fill_manual(values = cols) +
      theme_classic() +
      ggtitle(f)
    
    ggsave(filename = sprintf("custom_traceplots/M%s_WL_%s.jpg", M, f), height = 20, width = 40, units = "cm")
    
  }
}


