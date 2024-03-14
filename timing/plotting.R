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
                                      "gfBLUP genomic prediction", "Other"))

pal <- RColorBrewer::brewer.pal(9, "Paired")

ggplot(data = results.avg,
       mapping = aes(x = p, y = Duration, color = Step)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  theme_classic(base_size = 11) +
  ylab("Duration in seconds") +
  scale_color_manual(values = pal) +
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
        legend.key.width = unit(0.5, "cm"))

ggsave("timing/timing.jpg", width = 30, height = 10, units = "cm")


results.avg$percentage <- numeric(nrow(results.avg))

for (p in unique(results.avg$p)) {
  for (step in unique(results.avg$Step)[-9]) {
    results.avg[results.avg$p == p & results.avg$Step == step, "percentage"] <-
      results.avg[results.avg$p == p & results.avg$Step == step, "Duration"] / results.avg[results.avg$p == p & results.avg$Step == "Total", "Duration"]
  }
}

results.avg <- droplevels(results.avg[which(!(results.avg$Step == "Total")),])

ggplot(data = results.avg, mapping = aes(x = p, y = percentage, fill = Step)) +
  geom_area() +
  theme_classic(base_size = 11) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia")[2:9]) +
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
  scale_y_continuous(labels = scales::percent)

ggsave("timing/timing_percentages.jpg", width = 30, height = 10, units = "cm")





