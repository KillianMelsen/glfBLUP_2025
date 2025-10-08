# This script produces several figures of the manuscript supplementary material
# in the following order:
# - Figure S2B
# - Figure S2C
# - Figure S2A
# - Figure S7B
# - Figure S4
# - Figure S7C
# - Figure S5
# - Figure S6
# - Figure S7A
# - Figure S3
# To reproduce the figures, make sure that the datafiles loaded on line 36, 37,
# 101, 102, 175, 176, 233, 234, 316, 317, 423, and 424 are available.
# These datafiles can be produced using the script at
# `hyper_1415B5IR/analyses_single_date/MegaLMM_hyper_single_date.R`.
# Note that these figures may not look exactly like the ones in the manuscript
# due to the random nature of the MCMC procedures used by MegaLMM. This behaviour
# of MegaLMM is unfortunately outside of our control.

# Loading libraries:
library(rlist)
library(glfBLUP)
library(ggplot2)
library(grid)

# Setting seed:
set.seed(1997)

# Setting working directory:
wd <- getwd()
setwd(wd)

# M = 5 ====
# Calculating posterior means:
load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/LAMBDA_M5.RData")
load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/PREDS_M5.RData")
mean_Lambda <- MegaLMM::get_posterior_mean(Lambda_samples)
mean_pred <- MegaLMM::get_posterior_mean(pred_samples)

loadings <- as.data.frame(t(mean_Lambda))
names(loadings) <- paste0("F", 1:nrow(mean_Lambda))
loadings.wl <- loadings[-1,]
loadings.y <- loadings[1,]

loadings.wl$Wavelength <- as.numeric(rownames(loadings.wl))

loadings.wl.long <- tidyr::pivot_longer(loadings.wl, names_to = "Factor",
                                        cols = 1:(ncol(loadings.wl) - 1),
                                        values_to = "Loading")

loadings.wl.long$Factor <- factor(loadings.wl.long$Factor, levels = paste0("F", 1:nrow(mean_Lambda)))

# Some data for the spectral background:
conesdata <- read.csv("http://www.cvrl.org/database/data/cones/linss10e_5.csv")
names(conesdata) <- c("Wavelength", "Red", "Green", "Blue")
conesdata[is.na(conesdata)] <- 0
conesdata$colour <- rgb(conesdata$Red, conesdata$Green, conesdata$Blue, alpha = 0.6)
gradient <- t(conesdata$colour[conesdata$Wavelength >= 400 & conesdata$Wavelength <= 800])
g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Plotting:
ggplot(data = loadings.wl.long, mapping = aes(x = Wavelength, y = Loading, color = Factor)) +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = NatParksPalettes::natparks.pals("Acadia")[c(3, 4, 5, 7, 9)]) +
  ylim(c(-0.5, 1.1)) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_blank(),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  annotate("text", x = 757, y = 1.05, label = paste(("lambda(F1*', '* Y) * ' = ' *"), round(loadings.y["F1"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.75, label = paste(("lambda(F2*', '* Y) * ' = ' *"), round(loadings.y["F2"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.40, label = paste(("lambda(F3*', '* Y) * ' = ' *"), round(loadings.y["F3"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.225, label = paste(("lambda(F4*', '* Y) * ' = ' *"), round(loadings.y["F4"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.575, label = paste(("lambda(F5*', '* Y) * ' = ' *"), round(loadings.y["F5"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  xlab(NULL) +
  annotate("text", x = 400, y = 1, label = "B", color = "white", parse = TRUE, size = 4, hjust = 0)

ggsave("plots/SM_figure_S2B.png", width = 24, height = 7.5, units = "cm")

# M = 10 =======================================================================
# Setting seed:
set.seed(1997)
rm(Lambda_samples, pred_samples)

# Calculating posterior means:
load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/LAMBDA_M10.RData")
load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/PREDS_M10.RData")
mean_Lambda <- MegaLMM::get_posterior_mean(Lambda_samples)
mean_pred <- MegaLMM::get_posterior_mean(pred_samples)

loadings <- as.data.frame(t(mean_Lambda))
names(loadings) <- paste0("F", 1:nrow(mean_Lambda))
loadings.wl <- loadings[-1,]
loadings.y <- loadings[1,]

loadings.wl$Wavelength <- as.numeric(rownames(loadings.wl))

loadings.wl.long <- tidyr::pivot_longer(loadings.wl, names_to = "Factor",
                                        cols = 1:(ncol(loadings.wl) - 1),
                                        values_to = "Loading")

loadings.wl.long$Factor <- factor(loadings.wl.long$Factor, levels = paste0("F", 1:nrow(mean_Lambda)))

# Some data for the spectral background:
conesdata <- read.csv("http://www.cvrl.org/database/data/cones/linss10e_5.csv")
names(conesdata) <- c("Wavelength", "Red", "Green", "Blue")
conesdata[is.na(conesdata)] <- 0
conesdata$colour <- rgb(conesdata$Red, conesdata$Green, conesdata$Blue, alpha = 0.6)
gradient <- t(conesdata$colour[conesdata$Wavelength >= 400 & conesdata$Wavelength <= 800])
g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Plotting:
ggplot(data = loadings.wl.long, mapping = aes(x = Wavelength, y = Loading, color = Factor)) +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = NatParksPalettes::natparks.pals("Acadia", 10)) +
  ylim(c(-1.1, 1.1)) +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  annotate("text", x = 757, y = 0.75, label = paste(("lambda(F1*', '* Y) * ' = ' *"), round(loadings.y["F1"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  annotate("text", x = 757, y = -0.6, label = paste(("lambda(F2 *', '* Y) * ' = ' *"), round(loadings.y["F2"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  annotate("text", x = 757, y = 1, label = paste(("lambda(F3*', '* Y) * ' = ' *"), round(loadings.y["F3"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  annotate("text", x = 757, y = -0.3, label = paste(("lambda(F4*', '* Y) * ' = ' *"), round(loadings.y["F4"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  annotate("text", x = 757, y = -0.75, label = paste(("lambda(F5*', '* Y) * ' = ' *"), round(loadings.y["F5"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  annotate("text", x = 757, y = -0.9, label = paste(("lambda(F6*', '* Y) * ' = ' *"), round(loadings.y["F6"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  annotate("text", x = 757, y = 0.12, label = paste(("lambda(F7*', '* Y) * ' = ' *"), round(loadings.y["F7"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  annotate("text", x = 757, y = 0.6, label = paste(("lambda(F8*', '* Y) * ' = ' *"), round(loadings.y["F8"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  annotate("text", x = 757, y = 0.45, label = paste(("lambda(F9*', '* Y) * ' = ' *"), round(loadings.y["F9"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  annotate("text", x = 757, y = -0.45, label = paste(("lambda(F10*', '* Y) * ' = ' *"), round(loadings.y["F10"], 2)),
           color = "white", parse = TRUE, size = 3, hjust = 0) +
  xlab("Wavelength (nm)") +
  annotate("text", x = 400, y = 1, label = "C", color = "white", parse = TRUE, size = 4, hjust = 0)

ggsave("plots/SM_figure_S2C.png", width = 24, height = 8.3, units = "cm")

# M = 3 ========================================================================
# Setting seed:
set.seed(1997)
rm(Lambda_samples, pred_samples)

# Calculating posterior means:
load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/LAMBDA_M3.RData")
load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/PREDS_M3.RData")
mean_Lambda <- MegaLMM::get_posterior_mean(Lambda_samples)
mean_pred <- MegaLMM::get_posterior_mean(pred_samples)

loadings <- as.data.frame(t(mean_Lambda))
names(loadings) <- paste0("F", 1:nrow(mean_Lambda))
loadings.wl <- loadings[-1,]
loadings.y <- loadings[1,]

loadings.wl$Wavelength <- as.numeric(rownames(loadings.wl))

loadings.wl.long <- tidyr::pivot_longer(loadings.wl, names_to = "Factor",
                                        cols = 1:(ncol(loadings.wl) - 1),
                                        values_to = "Loading")

loadings.wl.long$Factor <- factor(loadings.wl.long$Factor, levels = paste0("F", 1:nrow(mean_Lambda)))

# Some data for the spectral background:
conesdata <- read.csv("http://www.cvrl.org/database/data/cones/linss10e_5.csv")
names(conesdata) <- c("Wavelength", "Red", "Green", "Blue")
conesdata[is.na(conesdata)] <- 0
conesdata$colour <- rgb(conesdata$Red, conesdata$Green, conesdata$Blue, alpha = 0.6)
gradient <- t(conesdata$colour[conesdata$Wavelength >= 400 & conesdata$Wavelength <= 800])
g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Plotting:
ggplot(data = loadings.wl.long, mapping = aes(x = Wavelength, y = Loading, color = Factor)) +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = NatParksPalettes::natparks.pals("Acadia", 10)[c(2, 8, 10)]) +
  ylim(c(-1.1, 1.1)) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_blank(),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  annotate("text", x = 757, y = -0.75, label = paste(("lambda(F1*', '* Y) * ' = ' *"), round(loadings.y["F1"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = -0.10, label = paste(("lambda(F2 *', '* Y) * ' = ' *"), round(loadings.y["F2"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.5, label = paste(("lambda(F3*', '* Y) * ' = ' *"), round(loadings.y["F3"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  xlab(NULL) +
  annotate("text", x = 400, y = 1, label = "A", color = "white", parse = TRUE, size = 4, hjust = 0)

ggsave("plots/SM_figure_S2A.png", width = 24, height = 8.3, units = "cm")

# M = 5 traceplotting ==========================================================
rm(list = ls())
if (!("Lambda_samples" %in% ls() & "pred_samples" %in% ls())) {
  load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/LAMBDA_M5.RData")
  load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/PREDS_M5.RData")
}

dimnames(Lambda_samples)[[2]] <- paste0("F", 1:dim(Lambda_samples)[2])
pred_samples <- pred_samples[, , 1]
dimnames(pred_samples)[[2]] <- substr(dimnames(pred_samples)[[2]], 1, 7)

loadings.y <- as.data.frame(Lambda_samples[, , "Y"])
loadings.y$Iteration <- 1:nrow(loadings.y)
loadings.y <- tidyr::pivot_longer(loadings.y, cols = 1:(ncol(loadings.y) - 1), names_to = "Factor", values_to = "Loading")
loadings.y$Factor <- factor(loadings.y$Factor, levels = paste0("F", 1:dim(Lambda_samples)[2]))

# Plotting the significant Y loadings:
ggplot(loadings.y, aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Factor), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia")[3:7]) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_blank(),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  xlab(NULL) +
  annotate("text", x = 1000, y = 0.4, label = "B", color = "black", parse = TRUE, size = 4, hjust = 0)

ggsave("plots/SM_figure_S7B.png", width = 24, height = 8, units = "cm")

# The important wavelengths around the "switching" point in glfBLUP:
(WL <- dimnames(Lambda_samples)[[3]][41:50])

factors <- paste0("F", 1:dim(Lambda_samples)[2])
niter <- dim(Lambda_samples)[1]
data <- as.data.frame(Lambda_samples[, factors[1], WL])
data$Factor <- rep(factors[1], niter)
data$Iteration <- 1:niter

for (f in factors[-1]) {
  subset <- as.data.frame(Lambda_samples[, f, WL])
  subset$Factor <- rep(f, niter)
  subset$Iteration <- 1:niter
  data <- rbind(data, subset)
}

data <- tidyr::pivot_longer(data, cols = 1:length(WL),
                            names_to = "Wavelength", values_to = "Loading")
data$Wavelength <- factor(data$Wavelength, levels = WL, labels = paste0(WL, " nm"))
data$Factor <- factor(data$Factor, levels = factors)

scaleFUN <- function(x) sprintf("%.2f", x)

p <- ggplot(data, aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Wavelength), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia", n = 10)) +
  theme_classic(base_size = 11) +
  facet_grid(rows = vars(Factor), scales = "free_y") +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  ylab("Loading") + xlab("Posterior sample (thin = 2)") +
  scale_y_continuous(labels=scaleFUN)

ggsave(plot = p, filename = "plots/SM_figure_S4.png",
       width = 24, height = 31.3, units = "cm")

# M = 10 traceplotting =========================================================
rm(list = ls())
if (!("Lambda_samples" %in% ls() & "pred_samples" %in% ls())) {
  load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/LAMBDA_M10.RData")
  load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/PREDS_M10.RData")
}

dimnames(Lambda_samples)[[2]] <- paste0("F", 1:dim(Lambda_samples)[2])
pred_samples <- pred_samples[, , 1]
dimnames(pred_samples)[[2]] <- substr(dimnames(pred_samples)[[2]], 1, 7)

loadings.y <- as.data.frame(Lambda_samples[, , "Y"])
loadings.y$Iteration <- 1:nrow(loadings.y)
loadings.y <- tidyr::pivot_longer(loadings.y, cols = 1:(ncol(loadings.y) - 1), names_to = "Factor", values_to = "Loading")
loadings.y$Factor <- factor(loadings.y$Factor, levels = paste0("F", 1:dim(Lambda_samples)[2]))

scaleFUN <- function(x) sprintf("%.2f", x)

# Plotting the significant Y loadings:
ggplot(loadings.y, aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Factor), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia", 10)) +
  theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  xlab("Posterior sample (thin = 2)") +
  scale_y_continuous(labels = scaleFUN) +
  annotate("text", x = 1000, y = 0.55, label = "C", color = "black", parse = TRUE, size = 4, hjust = 0)

ggsave("plots/SM_figure_S7C.png", width = 24, height = 8, units = "cm")

# The important wavelengths around the "switching" point in glfBLUP:
(WL <- dimnames(Lambda_samples)[[3]][41:50])

factors <- paste0("F", 1:dim(Lambda_samples)[2])
niter <- dim(Lambda_samples)[1]
data <- as.data.frame(Lambda_samples[, factors[1], WL])
data$Factor <- rep(factors[1], niter)
data$Iteration <- 1:niter

for (f in factors[-1]) {
  subset <- as.data.frame(Lambda_samples[, f, WL])
  subset$Factor <- rep(f, niter)
  subset$Iteration <- 1:niter
  data <- rbind(data, subset)
}

data <- tidyr::pivot_longer(data, cols = 1:length(WL),
                            names_to = "Wavelength", values_to = "Loading")
data$Wavelength <- factor(data$Wavelength, levels = WL, labels = paste0(WL, " nm"))
data$Factor <- factor(data$Factor, levels = factors)

scaleFUN <- function(x) sprintf("%.2f", x)

p1 <- ggplot(data[which(data$Factor %in% factors[1:5]),], aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Wavelength), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia", 10)) +
  theme_classic(base_size = 11) +
  facet_grid(rows = vars(Factor), scales = "free_y") +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  ylab("Loading") + xlab("Posterior sample (thin = 2)") +
  scale_y_continuous(labels=scaleFUN)

p2 <- ggplot(data[which(data$Factor %in% factors[6:10]),], aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Wavelength), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia", 10)) +
  theme_classic(base_size = 11) +
  facet_grid(rows = vars(Factor), scales = "free_y") +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  ylab("Loading") + xlab("Posterior sample (thin = 2)") +
  scale_y_continuous(labels=scaleFUN)

ggsave(plot = p1, filename = "plots/SM_figure_S5.png",
       width = 24, height = 31.3, units = "cm")

ggsave(plot = p2, filename = "plots/SM_figure_S6.png",
       width = 24, height = 31.3, units = "cm")

# M = 3 traceplotting ==========================================================
rm(list = ls())
if (!("Lambda_samples" %in% ls() & "pred_samples" %in% ls())) {
  load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/LAMBDA_M3.RData")
  load("hyper_1415B5IR/megalmm_hyper_single_date_arrays/PREDS_M3.RData")
}

dimnames(Lambda_samples)[[2]] <- paste0("F", 1:dim(Lambda_samples)[2])
pred_samples <- pred_samples[, , 1]
dimnames(pred_samples)[[2]] <- substr(dimnames(pred_samples)[[2]], 1, 7)

loadings.y <- as.data.frame(Lambda_samples[, , "Y"])
loadings.y$Iteration <- 1:nrow(loadings.y)
loadings.y <- tidyr::pivot_longer(loadings.y, cols = 1:(ncol(loadings.y) - 1), names_to = "Factor", values_to = "Loading")
loadings.y$Factor <- factor(loadings.y$Factor, levels = paste0("F", 1:dim(Lambda_samples)[2]))

scaleFUN <- function(x) sprintf("%.2f", x)

# Plotting the significant Y loadings:
ggplot(loadings.y, aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Factor), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia", 10)[c(2, 8, 10)]) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_blank(),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  xlab(NULL) +
  scale_y_continuous(labels=scaleFUN) +
  annotate("text", x = 1000, y = 0.62, label = "A", color = "black", parse = TRUE, size = 4, hjust = 0)

ggsave("plots/SM_figure_S7A.png", width = 24, height = 8, units = "cm")

# The important wavelengths around the "switching" point in glfBLUP:
(WL <- dimnames(Lambda_samples)[[3]][41:50])

factors <- paste0("F", 1:dim(Lambda_samples)[2])
niter <- dim(Lambda_samples)[1]
data <- as.data.frame(Lambda_samples[, factors[1], WL])
data$Factor <- rep(factors[1], niter)
data$Iteration <- 1:niter

for (f in factors[-1]) {
  subset <- as.data.frame(Lambda_samples[, f, WL])
  subset$Factor <- rep(f, niter)
  subset$Iteration <- 1:niter
  data <- rbind(data, subset)
}

data <- tidyr::pivot_longer(data, cols = 1:length(WL),
                            names_to = "Wavelength", values_to = "Loading")
data$Wavelength <- factor(data$Wavelength, levels = WL, labels = paste0(WL, " nm"))
data$Factor <- factor(data$Factor, levels = factors)

scaleFUN <- function(x) sprintf("%.2f", x)

p <- ggplot(data, aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Wavelength), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia", 10)) +
  theme_classic(base_size = 11) +
  facet_grid(rows = vars(Factor), scales = "free_y") +
  theme(axis.text = element_text(color = "black", size = 11),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  ylab("Loading") + xlab("Posterior sample (thin = 2)") +
  scale_y_continuous(labels=scaleFUN)

ggsave(plot = p, filename = "plots/SM_figure_S3.png",
       width = 24, height = 31.3, units = "cm")



