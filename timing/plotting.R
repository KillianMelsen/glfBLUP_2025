# Plotting the timing results:
library(ggplot2)
results <- read.csv("timing/timing.csv")
runs <- length(unique(results$Run))
results <- results[, -c(1, 5)]

totals <- expand.grid(Step = "Total",
                      Run = 1:runs,
                      p = seq(100, max(results$p), 100))

totals$Durations <- numeric(nrow(totals))

for (R in 1:runs) {
  for (P in seq(100, max(results$p), 100)) {
    temp <- results[which(results$Run == R & results$p == P), ]
    totals[which(totals$Run == R & totals$p == P), "Durations"] <- sum(temp$Durations)
  }
}; rm(temp, P, R)

results <- rbind(results, totals)
results.avg <- aggregate(Durations ~ Step + p, data = results, mean)

results.avg$Step <- factor(results.avg$Step,
                           levels = c("Total", "Redundancy filtering", "Regularization",
                                      "Factor model", "Factor scores", "Subset selection",
                                      "gfBLUP genomic prediction", "Other"))

pal <- RColorBrewer::brewer.pal(8, "Paired")

ggplot(data = results.avg,
       mapping = aes(x = p, y = Durations, color = Step)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  theme_dark() +
  ylab("Duration in seconds") +
  scale_color_manual(values = pal)

ggsave("timing/timing.jpg", width = 30, height = 10, units = "cm")


results.avg$percentage <- numeric(nrow(results.avg))

for (p in unique(results.avg$p)) {
  for (step in unique(results.avg$Step)[-8]) {
    results.avg[results.avg$p == p & results.avg$Step == step, "percentage"] <-
      results.avg[results.avg$p == p & results.avg$Step == step, "Durations"] / results.avg[results.avg$p == p & results.avg$Step == "Total", "Durations"] * 100
  }
}

results.avg <- droplevels(results.avg[which(!(results.avg$Step == "Total")),])

ggplot(data = results.avg, mapping = aes(x = p, y = percentage, fill = Step)) +
  geom_area() +
  theme_classic() +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia")[2:8])





