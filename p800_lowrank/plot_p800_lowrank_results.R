# This script produces figure 6 of the manuscript main text.
# Make sure that all intermediate result files loaded by the code block on lines 30 to
# 81 are available.

# Setting working directory:
wd <- getwd()
setwd(wd)

# Loading libraries:
library(ggplot2)

models <- data.frame(name = c("benchmark", "univariate", "glfblup", "MegaLMM", "lsblup", "siblup", "multiMLP"),
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
              
            accs <- read.csv(sprintf("p800_lowrank/results/h2s%s/%s%s_p800_lowrank_results_%s_%s_h2y%s_comm%s_h2s%s.csv",
                                     h2s, label, label.letter, model, CV, h2y, comm, h2s))$acc
            stopifnot(length(accs) == 20 | length(accs) == 100)
            stopifnot(length(unique(accs)) == 20 | length(unique(accs)) == 100)
            acc <- median(accs)
            
            medians[i, "Model"] <- model
            medians[i, "h2y"] <- as.numeric(h2y) / 10
            medians[i, "h2s"] <- paste0("h2s", h2s)
            medians[i, "comm"] <- paste0("comm", comm)
            medians[i, "CV"] <- CV
            medians[i, "Accuracy"] <- acc
            
            i <- i + 1
          }
        } else if (model == "univariate") {
          
          accs <- read.csv(sprintf("p800_lowrank/results/h2s%s/2_p800_lowrank_results_univariate_h2y%s_comm%s_h2s%s.csv", h2s, h2y, comm, h2s))$acc
          stopifnot(length(accs) == 20 | length(accs) == 100)
          stopifnot(length(unique(accs)) == 20 | length(unique(accs)) == 100)
          acc <- median(read.csv(sprintf("p800_lowrank/results/h2s%s/2_p800_lowrank_results_univariate_h2y%s_comm%s_h2s%s.csv", h2s, h2y, comm, h2s))$acc)
          
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

medians$Model <- factor(medians$Model, levels = c("Univariate", "benchmark", "glfblup", "MegaLMM", "lsblup", "siblup", "multiMLP"),
                        labels = c("Univariate", "Benchmark", "glfBLUP", "MegaLMM", "lsBLUP", "siBLUP", "multiMLP"))

h = 20
w = 25

labs <- list("comm02" = expression(bold("\u03A8")[~y]~"= 0.8"),
             "comm05" = expression(bold("\u03A8")[~y]~"= 0.5"),
             "comm08" = expression(bold("\u03A8")[~y]~"= 0.2"),
             "h2s05" = expression(h^2*(s)~"= 0.5"),
             "h2s07" = expression(h^2*(s)~"= 0.7"),
             "h2s09" = expression(h^2*(s)~"= 0.9"),
             "CV1" = "CV1", "CV2" = "CV2", "Univariate" = "Univariate")

labs_labeller <- function(variable, value) {
  return(labs[value])
}

ggplot(data = medians, mapping = aes(x = h2y, y = Accuracy, color = Model)) +
  theme_classic() +
  geom_hline(yintercept = seq(-0.05, 1, 0.05), color = "gray30", linetype = 1, linewidth = 0.1) +
  geom_line(mapping = aes(x = h2y, y = Accuracy, color = Model, linetype = CV), linewidth = 0.8) +
  geom_point(size = 3, mapping = aes(shape = CV)) +
  facet_grid(cols = vars(comm), rows = vars(h2s), labeller = labs_labeller) +
  xlab("Focal trait heritability") +
  labs(tag = "Secondary feature heritability") +
  ggtitle("Focal trait unique variance") +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold", size = 11, hjust = 0.5),
        legend.position = "right",
        legend.spacing.x = unit(0.1, "cm"),
        legend.text = element_text(size = 9),
        legend.key.size = unit(2,"line"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5),
        plot.tag.position = c(0.855, 0.48),
        plot.tag = element_text(angle = 270, face = "bold", size = 11),
        plot.margin = margin(l = 0.5, r = 1, t = 0.5, b = 0.5, unit = "cm"),
        strip.background = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("Univariate" = "#000000", "Benchmark" = "#E69F00",
                                "glfBLUP" = "#56B4E9", "MegaLMM" = "#009E73",
                                "siBLUP" = "#CC79A7", "lsBLUP" = "#D55E00",
                                "multiMLP" = "#0072B2")) +
  guides(Model = guide_legend(title.position = "top"),
         color = guide_legend(title.position = "top", nrow = 8, byrow = TRUE,
                              override.aes = list(shape = NA, linewidth = 2)),
         linetype = guide_legend(title.position = "top", nrow = 3, byrow = TRUE))

ggsave(filename = "plots/main_text_figure_6.png", dpi = 640, width = 25, height = 27.5, units = "cm")
