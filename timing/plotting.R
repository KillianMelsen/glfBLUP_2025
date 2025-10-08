# This script uses the timing measurents produced and saved by the script at
# `timing/timing.R` and uses them to produce figure S8 of the manuscript
# supplementary material.

# Plotting the timing results:
library(ggplot2)
results <- read.csv("timing/timing.csv")
runs <- length(unique(results$Run))
results <- results[, -1]

totals <- expand.grid(Step = "Total",
                      Run = 1:runs,
                      p = seq(100, max(results$p), 100))

totals$Duration <- numeric(nrow(totals))

for (R in 1:runs) {
  for (P in seq(100, max(results$p), 100)) {
    temp <- results[which(results$Run == R & results$p == P), ]
    totals[which(totals$Run == R & totals$p == P), "Duration"] <- sum(temp$Duration)
  }
}; rm(temp, P, R)

results <- rbind(results, totals)
results.avg <- aggregate(Duration ~ Step + p, data = results, mean)

results.avg$Step <- factor(results.avg$Step,
                           levels = c("Total", "Redundancy filtering", "Genetic regularization",
                                      "Residual regularization", "Factor model",
                                      "Factor scores", "Subset selection",
                                      "gfBLUP genomic prediction", "Other"),
                           labels = c("Total", "Redundancy filtering", "Genetic regularization",
                                      "Residual regularization", "Factor model",
                                      "Factor scores", "Subset selection",
                                      "glfBLUP genomic prediction", "Other"))

results.avg$Percentage <- numeric(nrow(results.avg))

for (p in unique(results.avg$p)) {
  for (step in unique(results.avg$Step)[-9]) {
    results.avg[results.avg$p == p & results.avg$Step == step, "Percentage"] <-
      results.avg[results.avg$p == p & results.avg$Step == step, "Duration"] / results.avg[results.avg$p == p & results.avg$Step == "Total", "Duration"]
  }
}

SF <- 1 / (max(results.avg[which(results.avg$Step == "Total"), "Duration"]) * 1.1)

ggplot(mapping = aes(x = p)) +
  geom_area(aes(y = Percentage, fill = Step), droplevels(results.avg[which(!(results.avg$Step == "Total")),])) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia")[2:9]) +
  geom_line(aes(y = Duration * SF), droplevels(results.avg[which(results.avg$Step == "Total"),]), col = "red3", linewidth = 1) +
  geom_point(aes(y = Duration * SF), droplevels(results.avg[which(results.avg$Step == "Total"),]), col = "red3", size = 1.5) +
  theme_classic(base_size = 11) +
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
  scale_y_continuous(labels = scales::percent, sec.axis = sec_axis(~./SF, name = "Total runtime (s)"))

ggsave("plots/SM_figure_S8.png", width = 24, height = 15, units = "cm")


