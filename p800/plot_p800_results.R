# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading libraries:
library(ggplot2)

models <- data.frame(name = c("benchmark", "univariate", "gfblup", "MegaLMM", "lsblup", "siblup", "multiMLP"),
                     label = c("1", "2", "3", "12", "11", "5", "8"))

h2.foc <- c("01", "03", "05", "07", "09")
h2.sec <- c("05", "07", "09")
comms <- c("02", "05", "08")

CVs <- c("CV1", "CV2")

n.row <- ((nrow(models) - 1) * length(h2.foc) * length(h2.sec) * length(comms) * length(CVs)) + 1 * (length(h2.foc) * length(h2.sec) * length(comms))

medians <- data.frame(Model = character(n.row),
                      h2y = numeric(n.row),
                      h2s = character(n.row),
                      comm = character(n.row),
                      CV = character(n.row),
                      Accuracy = numeric(n.row))

### Loading all results:
i <- 1; h2y <- h2.foc[1]; h2s <- h2.sec[1]; comm <- comms[1]; CV <- CVs[1]
for (model in models$name) {
  label <- models[which(models$name == model), "label"]
  for (h2y in h2.foc) {
    for (h2s in h2.sec) {
      for (comm in comms) {
        if (model != "univariate") {
          for (CV in CVs) {
            
            if (CV == "CV1") {
              label.letter <- "a"
            } else {
              label.letter <- "b"
            }
              
            if (model == "MegaLMM") {
              acc <- 0
            } else {
              accs <- read.csv(sprintf("p800/results/h2s%s/%s%s_p800_results_%s_%s_h2y%s_comm%s_h2s%s.csv",
                                       h2s, label, label.letter, model, CV, h2y, comm, h2s))$acc
              stopifnot(length(accs) == 20 | length(accs) == 100)
              acc <- median(accs)
            }
            
            medians[i, "Model"] <- model
            medians[i, "h2y"] <- as.numeric(h2y) / 10
            medians[i, "h2s"] <- paste0("h2s", h2s)
            medians[i, "comm"] <- paste0("comm", comm)
            medians[i, "CV"] <- CV
            medians[i, "Accuracy"] <- acc
            
            i <- i + 1
          }
        } else if (model == "univariate") {
          
          accs <- read.csv(sprintf("p800/results/h2s%s/2_p800_results_univariate_h2y%s_comm%s_h2s%s.csv", h2s, h2y, comm, h2s))$acc
          stopifnot(length(accs) == 20 | length(accs) == 100)
          acc <- median(read.csv(sprintf("p800/results/h2s%s/2_p800_results_univariate_h2y%s_comm%s_h2s%s.csv", h2s, h2y, comm, h2s))$acc)
          
          medians[i, "Model"] <- "Univariate"
          medians[i, "h2y"] <- as.numeric(h2y) / 10
          medians[i, "h2s"] <- paste0("h2s", h2s)
          medians[i, "comm"] <- paste0("comm", comm)
          medians[i, "CV"] <- "Univariate"
          medians[i, "Accuracy"] <- acc
          
          i <- i + 1
          
        }
      }
    }
  }
}

medians$Model <- factor(medians$Model, levels = c("Univariate", "benchmark", "gfblup", "MegaLMM", "lsblup", "siblup", "multiMLP"),
                        labels = c("Univariate", "Benchmark", "gfBLUP", "MegaLMM", "lsBLUP", "siBLUP", "multiMLP"))

h = 20
w = 25

labs <- list("comm02" = expression(bold("\u03A8")[~y]~"= 0.8"),
             "comm05" = expression(bold("\u03A8")[~y]~"= 0.5"),
             "comm08" = expression(bold("\u03A8")[~y]~"= 0.2"),
             "h2s05" = expression(h^2*(s)~"= 0.5"),
             "h2s07" = expression(h^2*(s)~"= 0.7"),
             "h2s09" = expression(h^2*(s)~"= 0.9"))

labs_labeller <- function(variable, value) {
  return(labs[value])
}

colors <- RColorBrewer::brewer.pal(7, "Dark2")
colors <- NatParksPalettes::natparks.pals("Charmonix")


ggplot(data = medians, mapping = aes(x = h2y, y = Accuracy, color = Model)) +
  theme_classic() +
  geom_point(size = 1) +
  geom_line(mapping = aes(x = h2y, y = Accuracy, color = Model, linetype = CV), linewidth = 0.8) +
  ylim(-0.1, 1) +
  xlab("Focal trait heritability") +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_linetype_manual(breaks = c("CV1", "CV2", "Univariate"),
                        values = c("CV1" = 1,
                                   "CV2" = 2,
                                   "Univariate" = 3)) +
  scale_color_manual(values = c("Univariate" = "black",
                                "Benchmark" = "red",
                                "gfBLUP" = colors[2],
                                "MegaLMM" = colors[3],
                                "siBLUP" = "purple4",
                                "lsBLUP" = colors[5],
                                "multiMLP" = colors[6])) +
  facet_grid(cols = vars(comm), rows = vars(h2s), labeller = labs_labeller) +
  theme(strip.background = element_blank()) +
  geom_hline(yintercept = seq(-0.1, 1, 0.1), color = "gray", linetype = 1, linewidth = 0.1) +
  theme(legend.position = "right",
        legend.spacing.x = unit(0.1, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11),
        panel.background = element_blank(),
        plot.title = element_text(size = 11),
        legend.key.size = unit(2,"line")) +
  guides(Model = guide_legend(title.position = "top", title.hjust = 0.5)) +
  labs(tag = "Secondary feature heritability") +
  theme(plot.tag.position = c(0.855, 0.48),
        plot.tag = element_text(angle = 270, face = "bold", size = 11),
        plot.margin = margin(l = 0.5, r = 1, t = 0.5, b = 0.5, unit = "cm")) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5,
                              nrow = 8, byrow = TRUE),
         linetype = guide_legend(title.position = "top", title.hjust = 0.5,
                                 nrow = 3, byrow = TRUE)) +
  ggtitle("Focal trait unique variance")

ggsave(filename = "plots/p800.png", dpi = 640, width = 25, height = 25, units = "cm")

# Without multiMLP:
ggplot(data = medians[which(medians$Model != "multiMLP"),], mapping = aes(x = h2y, y = Accuracy, color = Model)) +
  theme_classic() +
  geom_point(size = 1) +
  geom_line(mapping = aes(x = h2y, y = Accuracy, color = Model, linetype = CV), linewidth = 0.8) +
  ylim(0.3, 1) +
  xlab("Focal trait heritability") +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_linetype_manual(breaks = c("CV1", "CV2", "Univariate"),
                        values = c("CV1" = 1,
                                   "CV2" = 2,
                                   "Univariate" = 3)) +
  scale_color_manual(values = c("Univariate" = colors[1],
                                "Benchmark" = colors[2],
                                "gfBLUP" = colors[3],
                                "MegaLMM" = colors[4],
                                "siBLUP" = colors[5],
                                "lsBLUP" = colors[6])) +
  facet_grid(cols = vars(comm), rows = vars(h2s), labeller = labs_labeller) +
  theme(strip.background = element_blank()) +
  geom_hline(yintercept = seq(0.3, 1, 0.1), color = "gray", linetype = 1, linewidth = 0.1) +
  theme(legend.position = "right",
        legend.spacing.x = unit(0.1, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11),
        panel.background = element_blank(),
        plot.title = element_text(size = 11)) +
  guides(Model = guide_legend(title.position = "top", title.hjust = 0.5)) +
  labs(tag = "Secondary feature heritability") +
  theme(plot.tag.position = c(0.862, 0.48),
        plot.tag = element_text(angle = 270, face = "bold", size = 11),
        plot.margin = margin(l = 0.5, r = 1, t = 0.5, b = 0.5, unit = "cm")) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5,
                              nrow = 8, byrow = TRUE),
         linetype = guide_legend(title.position = "top", title.hjust = 0.5,
                                 nrow = 3, byrow = TRUE)) +
  ggtitle("Focal trait communality")

ggsave(filename = "plots/p800_2.png", dpi = 640, width = 25, height = 25, units = "cm")


