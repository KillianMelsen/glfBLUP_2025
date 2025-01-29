# Plotting the hyperspectral results:
library(ggplot2)
library(ggforce)
wd <- getwd()
setwd(wd)
set.seed(1997)

# Settings:
models <- c("Univariate", "gfBLUP", "MegaLMM", "lsBLUP", "siBLUP", "MultiMLP")
# models <- c("Univariate", "gfBLUP", "MegaLMM", "lsBLUP", "siBLUP")
scenarios <- c("CV1", "CV2", "CV2VEG")
irrigation <- c("HEAT", "B5IR")

# Loading results:
results <- expand.grid(Model = models[2:length(models)],
                       Scenario = scenarios,
                       Irrigation = irrigation,
                       Run = 1:250,
                       Accuracy = 0)

results <- rbind(results,
                 expand.grid(Model = models[1],
                             Scenario = "N/A",
                             Irrigation = "HEAT",
                             Run = 1:250,
                             Accuracy = 0))

results <- rbind(results,
                 expand.grid(Model = models[1],
                             Scenario = "N/A",
                             Irrigation = "B5IR",
                             Run = 1:250,
                             Accuracy = 0))

# Results for HEAT:
results[results$Model == "Univariate" & results$Scenario == "N/A" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/2_hyper_results_univariate.csv")$acc

results[results$Model == "gfBLUP" & results$Scenario == "CV1" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/3a_hyper_results_gfblup_CV1.csv")$acc
results[results$Model == "gfBLUP" & results$Scenario == "CV2" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/3b_hyper_results_gfblup_CV2.csv")$acc
results[results$Model == "gfBLUP" & results$Scenario == "CV2VEG" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/nosplines/3b_hyper_results_gfblup_CV2VEG.csv")$acc

results[results$Model == "MegaLMM" & results$Scenario == "CV1" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/12a_hyper_results_megalmm_CV1.csv")$acc
results[results$Model == "MegaLMM" & results$Scenario == "CV2" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/12b_hyper_results_megalmm_CV2.csv")$acc
results[results$Model == "MegaLMM" & results$Scenario == "CV2VEG" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/nosplines/12b_hyper_results_megalmm_CV2VEG.csv")$acc

results[results$Model == "lsBLUP" & results$Scenario == "CV1" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/11a_hyper_results_lsblup_CV1.csv")$acc
results[results$Model == "lsBLUP" & results$Scenario == "CV2" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/11b_hyper_results_lsblup_CV2.csv")$acc
results[results$Model == "lsBLUP" & results$Scenario == "CV2VEG" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/nosplines/11b_hyper_results_lsblup_CV2VEG.csv")$acc

results[results$Model == "siBLUP" & results$Scenario == "CV1" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/5a_hyper_results_siblup_CV1.csv")$acc
results[results$Model == "siBLUP" & results$Scenario == "CV2" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/5b_hyper_results_siblup_CV2.csv")$acc
results[results$Model == "siBLUP" & results$Scenario == "CV2VEG" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/nosplines/5b_hyper_results_siblup_CV2VEG.csv")$acc

results[results$Model == "MultiMLP" & results$Scenario == "CV1" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/8a_hyper_results_multiMLP_CV1.csv")$acc
results[results$Model == "MultiMLP" & results$Scenario == "CV2" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/splines/8b_hyper_results_multiMLP_CV2.csv")$acc
results[results$Model == "MultiMLP" & results$Scenario == "CV2VEG" & results$Irrigation == "HEAT", "Accuracy"] <-
  read.csv("hyper_1415HEAT/results/nosplines/8b_hyper_results_multiMLP_CV2VEG.csv")$acc

# Results for B5IR:
results[results$Model == "Univariate" & results$Scenario == "N/A" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/2_hyper_results_univariate.csv")$acc

results[results$Model == "gfBLUP" & results$Scenario == "CV1" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/3a_hyper_results_gfblup_CV1.csv")$acc
results[results$Model == "gfBLUP" & results$Scenario == "CV2" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/3b_hyper_results_gfblup_CV2.csv")$acc
results[results$Model == "gfBLUP" & results$Scenario == "CV2VEG" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/VEGsplines/3b_hyper_results_gfblup_CV2VEG.csv")$acc

results[results$Model == "MegaLMM" & results$Scenario == "CV1" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/12a_hyper_results_megalmm_CV1.csv")$acc
results[results$Model == "MegaLMM" & results$Scenario == "CV2" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/12b_hyper_results_megalmm_CV2.csv")$acc
results[results$Model == "MegaLMM" & results$Scenario == "CV2VEG" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/VEGsplines/12b_hyper_results_megalmm_CV2VEG.csv")$acc

results[results$Model == "lsBLUP" & results$Scenario == "CV1" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/11a_hyper_results_lsblup_CV1.csv")$acc
results[results$Model == "lsBLUP" & results$Scenario == "CV2" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/11b_hyper_results_lsblup_CV2.csv")$acc
results[results$Model == "lsBLUP" & results$Scenario == "CV2VEG" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/VEGsplines/11b_hyper_results_lsblup_CV2VEG.csv")$acc

results[results$Model == "siBLUP" & results$Scenario == "CV1" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/5a_hyper_results_siblup_CV1.csv")$acc
results[results$Model == "siBLUP" & results$Scenario == "CV2" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/5b_hyper_results_siblup_CV2.csv")$acc
results[results$Model == "siBLUP" & results$Scenario == "CV2VEG" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/VEGsplines/5b_hyper_results_siblup_CV2VEG.csv")$acc

results[results$Model == "MultiMLP" & results$Scenario == "CV1" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/8a_hyper_results_multiMLP_CV1.csv")$acc
results[results$Model == "MultiMLP" & results$Scenario == "CV2" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/splines/8b_hyper_results_multiMLP_CV2.csv")$acc
results[results$Model == "MultiMLP" & results$Scenario == "CV2VEG" & results$Irrigation == "B5IR", "Accuracy"] <-
  read.csv("hyper_1415B5IR/results/VEGsplines/8b_hyper_results_multiMLP_CV2VEG.csv")$acc

results$Model <- factor(results$Model, levels = c("Univariate", "gfBLUP", "MegaLMM", "siBLUP", "lsBLUP", "MultiMLP"))
# results$Model <- factor(results$Model, levels = c("Univariate", "gfBLUP", "MegaLMM", "siBLUP", "lsBLUP"))
results$Scenario <- factor(results$Scenario, levels = c("CV1", "CV2", "CV2VEG", "N/A"))
results$Irrigation <- factor(results$Irrigation, levels = c("HEAT", "B5IR"))

medians <- aggregate(Accuracy ~ Scenario + Model + Irrigation, data = results, FUN = median)
medians$max <- aggregate(Accuracy ~ Scenario + Model + Irrigation, data = results, FUN = max)$Accuracy

medians$Model <- factor(medians$Model, levels = c("Univariate", "gfBLUP", "MegaLMM", "siBLUP", "lsBLUP", "MultiMLP"))
# medians$Model <- factor(medians$Model, levels = c("Univariate", "gfBLUP", "MegaLMM", "siBLUP", "lsBLUP"))
medians$Scenario <- factor(medians$Scenario, levels = c("CV1", "CV2", "CV2VEG", "N/A"))
medians$Irrigation <- factor(medians$Irrigation, levels = c("HEAT", "B5IR"))

dummy <- data.frame(Accuracy = c(0.0, 0.8, 0.0, 1.0),
                    Irrigation = c("HEAT", "HEAT", "B5IR", "B5IR"),
                    Model = "Univariate",
                    Scenario = "N/A")

dummy$Irrigation <- factor(dummy$Irrigation, levels = c("HEAT", "B5IR"))

# Plotting with CV2VEG results:
ggplot(data = results, mapping = aes(x = Model, y = Accuracy, fill = Scenario)) +
  facet_wrap(vars(Irrigation), ncol = 1, scales = "free_y") +
  stat_boxplot(geom = "errorbar", linewidth = 0.25) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25) +
  geom_text(data = medians, aes(label = round(Accuracy, 2), y = max + 0.05),
            position = position_dodge(width = 0.75), size = 2, color = "gray45") +
  theme_classic(base_size = 11) +
  scale_fill_manual(breaks = c("CV1", "CV2", "CV2VEG"),
                    values = c("CV1" = "#63a7ff", "CV2" = "#ffb667",
                               "NA" = "gray", "CV2VEG" = "yellowgreen")) +
  # ylim(min(min(results$Accuracy), 0), 1) +
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
            size = 0.3, position = position_dodge(width = 0.75), shape = 21, stroke = 0.05, maxwidth = 0.5, alpha = 0.25) +
  geom_blank(data = dummy)


ggsave(filename = "plots/hyper.png", dpi = 640, width = 15, height = 12, units = "cm")



