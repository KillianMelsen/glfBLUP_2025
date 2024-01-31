library(ggplot2)
results <- read.csv("timing.csv")
results <- results[which(results$p %in% seq(100, 1400, 100)), -c(1, 5)]

totals <- expand.grid(Step = "Total",
                      Run = 1:3,
                      p = seq(100, 1400, 100))

totals$Durations <- numeric(nrow(totals))

for (R in 1:3) {
  for (P in seq(100, 1400, 100)) {
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


ggplot(data = results[which(results$Step == "Total"), ],
       mapping = aes(x = p, y = Durations)) +
  geom_point() +
  theme_classic() +
  ylab("Duration in seconds")

pal <- RColorBrewer::brewer.pal(8, "Paired")

ggplot(data = results.avg,
       mapping = aes(x = p, y = Durations, color = Step)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  theme_dark() +
  ylab("Duration in seconds") +
  scale_color_manual(values = pal)

ggsave("timing.jpg", width = 30, height = 10, units = "cm")
