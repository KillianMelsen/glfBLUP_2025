set.seed(1997)
results <-
  rlist::list.load("hyper_1415HEAT/results/splines/11a_hyper_extra_results_lsblup_CV1.RData")

i <- sample(1:250, 1)
results <- results[[i]]
results <- names(results$coefs)

temp <- colnames(rlist::list.load("hyper_1415HEAT/datasets/splines/hyper_dataset_1.RData")$data)
temp <- temp[2:(length(temp) - 1)]
wl <- as.numeric(unique(substr(temp, 3, 5)))
dates <- unique(substr(temp, 7, 12))
rm(temp, i)

data <- expand.grid(Date = dates,
                    Wavelength = wl,
                    Retained = FALSE)

for (i in results) {
  wavelength <- as.numeric(substr(i, 3, 5))
  date <- substr(i, 7, 12)
  data[which(data$Date == date & data$Wavelength == wavelength), "Retained"] <- TRUE
}

library(ggplot2)
library(grid)
data$Date <- factor(data$Date, levels = levels(data$Date), labels = c("2015-04-14", "2015-04-24",
                                                                      "2015-04-28", "2015-05-06"))

# Some data for the spectral background:
conesdata <- read.csv("http://www.cvrl.org/database/data/cones/linss10e_5.csv")
names(conesdata) <- c("Wavelength", "Red", "Green", "Blue")
conesdata[is.na(conesdata)] <- 0
conesdata$colour <- rgb(conesdata$Red, conesdata$Green, conesdata$Blue, alpha = 0.6)
gradient <- t(conesdata$colour[conesdata$Wavelength >= 400 & conesdata$Wavelength <= 800])
g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

ggplot(data, aes(x = Wavelength, y = 0, color = Retained)) +
  facet_wrap(vars(Date), ncol = 2) +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_point(size = 1) +
  scale_color_manual(values = c("FALSE" = "red", "TRUE" = "green")) +
  theme_classic() +
  xlab("Wavelength (nm)") + ylab(NULL) +
  theme(legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

ggsave(filename = "plots/hyper_RF_1415HEAT.png", dpi = 640, width = 15, height = 6, units = "cm")



