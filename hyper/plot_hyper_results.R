# Plotting the hyperspectral results:
library(ggplot2)
library(ggforce)
wd <- getwd()
setwd(wd)

# Settings:
# models <- c("Univariate", "Phenomic", "gfBLUP", "MegaLMM", "lsBLUP", "siBLUP", "MultiMLP")
models <- c("Univariate", "gfBLUP", "MegaLMM", "lsBLUP", "siBLUP", "MultiMLP")
scenarios <- c("CV1", "CV2", "CV2VEG")

# Loading results:
results <- expand.grid(Model = models[2:length(models)],
                       Scenario = scenarios,
                       Run = 1:250,
                       Accuracy = 0)

results <- rbind(results,
                 expand.grid(Model = models[1],
                             Scenario = "N/A",
                             Run = 1:250,
                             Accuracy = 0))

# results <- rbind(results,
#                  expand.grid(Model = models[2],
#                              Scenario = scenarios[2:3],
#                              Run = 1:250,
#                              Accuracy = 0))

results[results$Model == "Univariate" & results$Scenario == "N/A", 4] <- read.csv("hyper/results/2_hyper_results_univariate.csv")$acc
# results[results$Model == "Phenomic" & results$Scenario == "CV2", 4] <- read.csv("hyper/results/0_hyper_results_phenomic.csv")$acc
# results[results$Model == "Phenomic" & results$Scenario == "CV2VEG", 4] <- read.csv("hyper/results/0_hyper_results_phenomic_CV2VEG.csv")$acc

results[results$Model == "gfBLUP" & results$Scenario == "CV1", 4] <- read.csv("hyper/results/3a_hyper_results_gfblup_CV1_RF.csv")$acc
results[results$Model == "gfBLUP" & results$Scenario == "CV2", 4] <- read.csv("hyper/results/3b_hyper_results_gfblup_CV2_RF.csv")$acc
results[results$Model == "gfBLUP" & results$Scenario == "CV2VEG", 4] <- read.csv("hyper/results/3b_hyper_results_gfblup_CV2VEG.csv")$acc

results[results$Model == "MegaLMM" & results$Scenario == "CV1", 4] <- read.csv("hyper/results/12a_hyper_results_megalmm_CV1_RF.csv")$acc
results[results$Model == "MegaLMM" & results$Scenario == "CV2", 4] <- read.csv("hyper/results/12b_hyper_results_megalmm_CV2_RF.csv")$acc
results[results$Model == "MegaLMM" & results$Scenario == "CV2VEG", 4] <- read.csv("hyper/results/12b_hyper_results_megalmm_CV2VEG_RF.csv")$acc

results[results$Model == "lsBLUP" & results$Scenario == "CV1", 4] <- read.csv("hyper/results/11a_hyper_results_lsblup_CV1.csv")$acc
results[results$Model == "lsBLUP" & results$Scenario == "CV2", 4] <- read.csv("hyper/results/11b_hyper_results_lsblup_CV2.csv")$acc
results[results$Model == "lsBLUP" & results$Scenario == "CV2VEG", 4] <- read.csv("hyper/results/11b_hyper_results_lsblup_CV2VEG.csv")$acc

results[results$Model == "siBLUP" & results$Scenario == "CV1", 4] <- read.csv("hyper/results/5a_hyper_results_siblup_CV1.csv")$acc
results[results$Model == "siBLUP" & results$Scenario == "CV2", 4] <- read.csv("hyper/results/5b_hyper_results_siblup_CV2.csv")$acc
results[results$Model == "siBLUP" & results$Scenario == "CV2VEG", 4] <- read.csv("hyper/results/5b_hyper_results_siblup_CV2VEG.csv")$acc

results[results$Model == "MultiMLP" & results$Scenario == "CV1", 4] <- read.csv("hyper/results/8a_hyper_results_multiMLP_CV1.csv")$acc
results[results$Model == "MultiMLP" & results$Scenario == "CV2", 4] <- read.csv("hyper/results/8b_hyper_results_multiMLP_CV2.csv")$acc
results[results$Model == "MultiMLP" & results$Scenario == "CV2VEG", 4] <- read.csv("hyper/results/8b_hyper_results_multiMLP_CV2VEG.csv")$acc

# results$Model <- factor(results$Model, levels = c("Univariate", "Phenomic", "gfBLUP", "MegaLMM", "siBLUP", "lsBLUP", "MultiMLP"))
results$Model <- factor(results$Model, levels = c("Univariate", "gfBLUP", "MegaLMM", "siBLUP", "lsBLUP", "MultiMLP"))
results$Scenario <- factor(results$Scenario, levels = c("CV1", "CV2", "CV2VEG", "N/A"))

medians <- aggregate(Accuracy ~ Scenario + Model, data = results, FUN = median)
medians$max <- aggregate(Accuracy ~ Scenario + Model, data = results, FUN = max)$Accuracy

# medians$Model <- factor(medians$Model, levels = c("Univariate", "Phenomic", "gfBLUP", "MegaLMM", "siBLUP", "lsBLUP", "MultiMLP"))
medians$Model <- factor(medians$Model, levels = c("Univariate", "gfBLUP", "MegaLMM", "siBLUP", "lsBLUP", "MultiMLP"))
medians$Scenario <- factor(medians$Scenario, levels = c("CV1", "CV2", "CV2VEG", "N/A"))

# Plotting:
ggplot(data = results, mapping = aes(x = Model, y = Accuracy, fill = Scenario)) +
  stat_boxplot(geom = "errorbar", linewidth = 0.25) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25) +
  geom_text(data = medians, aes(label = round(Accuracy, 2), y = max + 0.05),
            position = position_dodge(width = 0.75), size = 2, color = "gray45") +
  theme_classic(base_size = 11) +
  scale_fill_manual(breaks = c("CV1", "CV2", "CV2VEG"),
                    values = c("CV1" = "#63a7ff", "CV2" = "#ffb667",
                               "NA" = "gray", "CV2VEG" = "yellowgreen")) +
  ylim(min(min(results$Accuracy), 0), 1.0) +
  theme(axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  geom_sina(data = results, mapping = aes(x = Model, y = Accuracy, fill = Scenario),
            size = 0.3, position = position_dodge(width = 0.75), shape = 21, stroke = 0.05, maxwidth = 0.5, alpha = 0.25)

ggsave(filename = "plots/hyper_CV2VEG.png", dpi = 640, width = 15, height = 6, units = "cm")



# Plotting:
ggplot(data = droplevels(results[which(results$Scenario != "CV2VEG"),]), mapping = aes(x = Model, y = Accuracy, fill = Scenario)) +
  stat_boxplot(geom = "errorbar", linewidth = 0.25) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25) +
  geom_text(data = droplevels(medians[which(medians$Scenario != "CV2VEG"),]), aes(label = round(Accuracy, 2), y = max + 0.05),
            position = position_dodge(width = 0.75), size = 2, color = "gray45") +
  theme_classic(base_size = 11) +
  scale_fill_manual(breaks = c("CV1", "CV2"),
                    values = c("CV1" = "#63a7ff", "CV2" = "#ffb667",
                               "NA" = "gray")) +
  ylim(min(min(results$Accuracy), 0), 1.0) +
  theme(axis.text = element_text(color = "black", size = 6),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  geom_sina(data = droplevels(results[which(results$Scenario != "CV2VEG"),]), mapping = aes(x = Model, y = Accuracy, fill = Scenario),
            size = 0.3, position = position_dodge(width = 0.75), shape = 21, stroke = 0.05, maxwidth = 0.5, alpha = 0.25)

ggsave(filename = "plots/hyper.png", dpi = 640, width = 15, height = 6, units = "cm")





