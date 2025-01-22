# Loading libraries:
library(rlist)
library(gfBLUP)
library(ggplot2)
library(grid)

# Setting seed:
set.seed(1997)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading kinship:
load("genotypes/K_hyper.RData")

# Which dataset should we look at:
dataset <- 1
      
# Loading hyperspectral dataset:
datalist <- list.load(file = sprintf("hyper/datasets/hyper_dataset_%d.RData", dataset))

# Storing data:
d <- datalist$data

# Subsetting to 10-03-2015:
d <- d[c(1, grep(".*_150310", names(d)), ncol(d))]

# Make training data and store feature/focal trait names:
d.train <- droplevels(na.omit(d))
sec <- names(d[2:(ncol(d) - 1)])
foc <- names(d)[ncol(d)]

# Regularization:
folds <- gfBLUP::createFolds(genos = unique(as.character(d$G)))
tempG <- gfBLUP::regularizedCorrelation(data = d[c("G", sec)], folds = folds, what = "genetic", dopar = TRUE, verbose = FALSE)
tempE <- gfBLUP::regularizedCorrelation(data = d[c("G", sec)], folds = folds, what = "residual", dopar = TRUE, verbose = FALSE)
Rg.reg <- tempG$optCor

# Fitting factor model:
FM.fit <- gfBLUP::factorModel(data = d[c("G", sec)], cormat = Rg.reg, what = "genetic", verbose = FALSE)

# Getting factor scores:
# scaling back loadings and uniquenesses:
D <- sqrt(diag(tempG$Sg))
L.cov <- diag(D) %*% FM.fit$loadings
PSI.cov <- outer(D, D) * FM.fit$uniquenesses

# Determining factor scores:
d[sec] <- sapply(d[sec], scale)
scores <- gfBLUP::factorScores(data = d[c("G", sec)],
                               loadings = L.cov,
                               uniquenesses = PSI.cov,
                               m = FM.fit$m,
                               type = "genetic-thomson-repdiv",
                               Se = outer(sqrt(diag(tempE$Se)), sqrt(diag(tempE$Se))) * tempE$optCor)

d <- cbind(scores, d$Y)
names(d)[ncol(d)] <- "Y"
names(d)[1] <- "G"

# Getting the genetic correlations and loadings:
gencors <- cov2cor(gfBLUP::covSS(na.omit(d))$Sg)
loadings <- as.data.frame(FM.fit$loadings[,])
rownames(loadings) <- substr(rownames(loadings), 3, 5)
loadings$Wavelength <- rownames(loadings)
loadings <- tidyr::pivot_longer(loadings, 1:3, names_to = "Factor", values_to = "Loading")
loadings$Wavelength <- as.numeric(loadings$Wavelength)
loadings$Factor <- factor(loadings$Factor, levels = c("F1", "F2", "F3"))

# Some data for the spectral background:
conesdata <- read.csv("http://www.cvrl.org/database/data/cones/linss10e_5.csv")
names(conesdata) <- c("Wavelength", "Red", "Green", "Blue")
conesdata[is.na(conesdata)] <- 0
conesdata$colour <- rgb(conesdata$Red, conesdata$Green, conesdata$Blue, alpha = 0.6)
gradient <- t(conesdata$colour[conesdata$Wavelength >= 400 & conesdata$Wavelength <= 800])
g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# gradient <- t(photobiology::w_length2rgb(400:800))
# g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Plotting:
ggplot(data = loadings, mapping = aes(x = Wavelength, y = Loading, color = Factor)) +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = c("F1" = "#63a7ff", "F2" = "#ffb667", "F3" = "yellowgreen")) +
  ylim(c(-0.5, 1)) +
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
  annotate("text", x = 757, y = 0.8, label = paste(("rho[(F2*', '* Y)]^g * ' = ' *"), round(gencors["F2", "Y"], 2)),
           color = "white", parse = TRUE, size = 6, hjust = 0) +
  annotate("text", x = 757, y = -0.3, label = paste(("rho[(F1*', '* Y)]^g * ' = ' *"), round(gencors["F1", "Y"], 2)),
           color = "white", parse = TRUE, size = 6, hjust = 0) +
  annotate("text", x = 757, y = 0.2, label = paste(("rho[(F3*', '* Y)]^g * ' = ' *"), round(gencors["F3", "Y"], 2)),
           color = "white", parse = TRUE, size = 6, hjust = 0)

ggsave("plots/gfBLUP_hyper_single_date.png", width = 24, height = 8, units = "cm")


