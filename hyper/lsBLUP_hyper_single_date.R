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
sec <- names(d[2:(ncol(d) - 1)])
foc <- names(d)[ncol(d)]

# Run siBLUP:
RESULT <- gfBLUPold::lsBLUP(d = d, K = K, CV = "CV2", do.parallel = FALSE, sepExp = FALSE, t.RF = 1, verbose = TRUE)

# Getting the SI-Y genetic correlation and SI weights (gamma):
coefs <- data.frame(Coefficient = as.numeric(RESULT$coefs))
rownames(coefs) <- names(RESULT$coefs)
gencors <- cov2cor(RESULT$Vg)
rownames(coefs) <- substr(rownames(coefs), 3, 5)
coefs$Wavelength <- rownames(coefs)
coefs$Wavelength <- as.numeric(coefs$Wavelength)

# Some data for the spectral background:
# conesdata <- read.csv("http://www.cvrl.org/database/data/cones/linss10e_5.csv")
# names(conesdata) <- c("Wavelength", "Red", "Green", "Blue")
# conesdata[is.na(conesdata)] <- 0
# conesdata$colour <- rgb(conesdata$Red, conesdata$Green, conesdata$Blue, alpha = 0.8)   
# gradient <- t(conesdata$colour[conesdata$Wavelength >= 400 & conesdata$Wavelength <= 800])
# g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

gradient <- t(photobiology::w_length2rgb(400:800))
g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Plotting:
ggplot(data = coefs, mapping = aes(x = Wavelength, y = Coefficient)) +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_line(linewidth = 1.5, color = "gray") +
  ylim(c(-1.5, 1.5)) +
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
  annotate("text", x = 757, y = -0.8, label = paste(("rho[(LSP*', '* Y)]^g * ' = ' *"), round(gencors["LSP", "Y"], 2)),
           color = "white", parse = TRUE, size = 6, hjust = 0) +
  ylab("LSP Coefficient")

ggsave("plots/lsBLUP_hyper_single_date.png", width = 24, height = 8, units = "cm")


