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

# M = 5 ========================================================================
# Some MegaLMM parameters:
M <- 5
burnin <- 10000
posterior <- 100000
thin <- 2

# Data manipulation:
d.train <- droplevels(d[which(!is.na(d$Y)), ])
d.test <- droplevels(d[which(is.na(d$Y)), ])

# Calculating genotypic means (BLUEs):
d.train <- gfBLUP:::genotypeMeans(d.train)
d.test <- gfBLUP:::genotypeMeans(d.test)

# Rescaling now we only have means:
d.train[, 2:ncol(d.train)] <- sapply(d.train[, 2:ncol(d.train)], scale)
d.test[, 2:ncol(d.test)] <- sapply(d.test[, 2:ncol(d.test)], scale)

d <- rbind(d.test, d.train)

# Reordering the columns (G, Y, SEC):
d <- cbind(d["G"], d["Y"], d[sec])
names(d)[3:ncol(d)] <- substr(names(d)[3:ncol(d)], 3, 5)

# MegaLMM config:
run_parameters <- MegaLMM::MegaLMM_control(
  drop0_tol = 1e-10,
  scale_Y = FALSE,
  h2_divisions = 20,
  h2_step_size = NULL,
  burn = 0,
  K = M,
  save_current_state = TRUE,
  thin = thin
)

priors = MegaLMM::MegaLMM_priors(
  tot_Y_var = list(V = 0.5, nu = 10),
  tot_F_var = list(V = 18/20, nu = 20),
  Lambda_prior = list(
    sampler = MegaLMM::sample_Lambda_prec_horseshoe,
    prop_0 = 0.1,
    delta = list(shape = 3, scale = 1),
    delta_iterations_factor = 100
  ),
  h2_priors_resids_fun = function(h2s, n) 1,
  h2_priors_factors_fun = function(h2s, n) 1
)

# Creating run ID:
run_ID <- "hyper/megalmm_states/MegaLMM_hyper_single_date_M5"

# Initializing MegaLMM:
MegaLMM_state = MegaLMM::setup_model_MegaLMM(d[, 2:ncol(d)],
                                             ~ 1 + (1|G),
                                             data = d,
                                             relmat = list(G = K),
                                             run_parameters = run_parameters,
                                             run_ID = run_ID)

maps = MegaLMM::make_Missing_data_map(MegaLMM_state, verbose = T)
MegaLMM_state <- MegaLMM::set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)

MegaLMM_state <- MegaLMM::set_priors_MegaLMM(MegaLMM_state, priors)
MegaLMM_state <- MegaLMM::initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state <- MegaLMM::initialize_MegaLMM(MegaLMM_state, verbose = T)

MegaLMM_state$Posterior$posteriorSample_params <- c("Lambda")
MegaLMM_state$Posterior$posteriorFunctions <- list(pred = "U_R[,1] + U_F %*% Lambda[,1]")

# Clearing posterior samples if they already exist for some reason:
if (file.exists(run_ID)) {
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
}

# Burn-in, collect samples, reorder factors, clear posterior, repeat 10 times
n_iter <- burnin / 10
n_burn_in <- 10
for (i in 1:n_burn_in) {
  MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state, drop_cor_threshold = 0.6)
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
}

# Clearing the burn-in samples:
MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)

# Collecting the posterior samples in a single go:
n_iter <- posterior
n_sampling <- 1
for (i in 1:n_sampling) {
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
  MegaLMM_state <- MegaLMM::save_posterior_chunk(MegaLMM_state)
}

# Reloading the saved posterior samples:
Lambda_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "Lambda")
pred_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "pred")

# Saving the arrays containing all posterior predictions and loadings
save(Lambda_samples, file = "hyper/megalmm_hyper_single_date_arrays/LAMBDA_M5.RData")
save(pred_samples, file = "hyper/megalmm_hyper_single_date_arrays/PREDS_M5.RData")

# Deleting MegaLMM state files:
unlink(run_ID, recursive = TRUE)

# Calculating posterior means:
mean_Lambda <- MegaLMM::get_posterior_mean(Lambda_samples)
mean_pred <- MegaLMM::get_posterior_mean(pred_samples)

loadings <- as.data.frame(t(mean_Lambda))
names(loadings) <- paste0("F", 1:nrow(mean_Lambda))
loadings.wl <- loadings[-1,]
loadings.y <- loadings[1,]

loadings.wl$Wavelength <- as.numeric(rownames(loadings.wl))

loadings.wl.long <- tidyr::pivot_longer(loadings.wl, names_to = "Factor",
                                        cols = 1:(ncol(loadings.wl) - 1),
                                        values_to = "Loading")

loadings.wl.long$Factor <- factor(loadings.wl.long$Factor, levels = paste0("F", 1:nrow(mean_Lambda)))

# Some data for the spectral background:
conesdata <- read.csv("http://www.cvrl.org/database/data/cones/linss10e_5.csv")
names(conesdata) <- c("Wavelength", "Red", "Green", "Blue")
conesdata[is.na(conesdata)] <- 0
conesdata$colour <- rgb(conesdata$Red, conesdata$Green, conesdata$Blue, alpha = 0.6)
gradient <- t(conesdata$colour[conesdata$Wavelength >= 400 & conesdata$Wavelength <= 800])
g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# gradient <- t(photobiology::w_length2rgb(400:800))
# g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Plotting:
ggplot(data = loadings.wl.long, mapping = aes(x = Wavelength, y = Loading, color = Factor)) +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = NatParksPalettes::natparks.pals("Acadia")[c(3, 4, 5, 7, 9)]) +
  ylim(c(-0.5, 1.1)) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(color = "black", size = 11),
        axis.text.x = element_blank(),
        axis.title = element_text(face = "bold", size = 13),
        axis.title.y.left = element_text(margin = margin(r = 0.25, unit = "cm")),
        axis.title.y.right = element_text(margin = margin(l = 0.25, unit = "cm")),
        legend.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11),
        legend.position = "right",
        legend.key.height = unit(0.5, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(0.5, "cm")) +
  annotate("text", x = 757, y = 0.13, label = paste(("lambda(F1*', '* Y) * ' = ' *"), round(loadings.y["F1"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = -0.5, label = paste(("lambda(F2*', '* Y) * ' = ' *"), round(loadings.y["F2"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = -0.35, label = paste(("lambda(F3*', '* Y) * ' = ' *"), round(loadings.y["F3"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.45, label = paste(("lambda(F4*', '* Y) * ' = ' *"), round(loadings.y["F4"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.90, label = paste(("lambda(F5*', '* Y) * ' = ' *"), round(loadings.y["F5"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  xlab(NULL) +
  annotate("text", x = 400, y = 1, label = "A", color = "white", parse = TRUE, size = 4, hjust = 0)

ggsave("plots/MegaLMM_hyper_single_date_M5.png", width = 24, height = 7.5, units = "cm")

# M = 10 =======================================================================
# Setting seed:
set.seed(1997)

# Some MegaLMM parameters:
M <- 10
burnin <- 10000
posterior <- 100000
thin <- 2

# MegaLMM config:
run_parameters <- MegaLMM::MegaLMM_control(
  drop0_tol = 1e-10,
  scale_Y = FALSE,
  h2_divisions = 20,
  h2_step_size = NULL,
  burn = 0,
  K = M,
  save_current_state = TRUE,
  thin = thin
)

priors = MegaLMM::MegaLMM_priors(
  tot_Y_var = list(V = 0.5, nu = 10),
  tot_F_var = list(V = 18/20, nu = 20),
  Lambda_prior = list(
    sampler = MegaLMM::sample_Lambda_prec_horseshoe,
    prop_0 = 0.1,
    delta = list(shape = 3, scale = 1),
    delta_iterations_factor = 100
  ),
  h2_priors_resids_fun = function(h2s, n) 1,
  h2_priors_factors_fun = function(h2s, n) 1
)

# Creating run ID:
run_ID <- "hyper/megalmm_states/MegaLMM_hyper_single_date_M10"

# Initializing MegaLMM:
MegaLMM_state = MegaLMM::setup_model_MegaLMM(d[, 2:ncol(d)],
                                             ~ 1 + (1|G),
                                             data = d,
                                             relmat = list(G = K),
                                             run_parameters = run_parameters,
                                             run_ID = run_ID)

maps = MegaLMM::make_Missing_data_map(MegaLMM_state, verbose = T)
MegaLMM_state <- MegaLMM::set_Missing_data_map(MegaLMM_state, maps$Missing_data_map)

MegaLMM_state <- MegaLMM::set_priors_MegaLMM(MegaLMM_state, priors)
MegaLMM_state <- MegaLMM::initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state <- MegaLMM::initialize_MegaLMM(MegaLMM_state, verbose = T)

MegaLMM_state$Posterior$posteriorSample_params <- c("Lambda")
MegaLMM_state$Posterior$posteriorFunctions <- list(pred = "U_R[,1] + U_F %*% Lambda[,1]")

# Clearing posterior samples if they already exist for some reason:
if (file.exists(run_ID)) {
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
}

# Burn-in, collect samples, reorder factors, clear posterior, repeat 10 times
n_iter <- burnin / 10
n_burn_in <- 10
for (i in 1:n_burn_in) {
  MegaLMM_state <- MegaLMM::reorder_factors(MegaLMM_state, drop_cor_threshold = 0.6)
  MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
}

# Clearing the burn-in samples:
MegaLMM_state <- MegaLMM::clear_Posterior(MegaLMM_state)

# Collecting the posterior samples in a single go:
n_iter <- posterior
n_sampling <- 1
for (i in 1:n_sampling) {
  MegaLMM_state <- MegaLMM::sample_MegaLMM(MegaLMM_state, n_iter, grainSize = 1)
  MegaLMM_state <- MegaLMM::save_posterior_chunk(MegaLMM_state)
}

# Reloading the saved posterior samples:
Lambda_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "Lambda")
pred_samples <- MegaLMM::load_posterior_param(MegaLMM_state, "pred")

# Saving the arrays containing all posterior predictions and loadings
save(Lambda_samples, file = "hyper/megalmm_hyper_single_date_arrays/LAMBDA_M10.RData")
save(pred_samples, file = "hyper/megalmm_hyper_single_date_arrays/PREDS_M10.RData")

# Deleting MegaLMM state files:
unlink(run_ID, recursive = TRUE)

# Calculating posterior means:
mean_Lambda <- MegaLMM::get_posterior_mean(Lambda_samples)
mean_pred <- MegaLMM::get_posterior_mean(pred_samples)

loadings <- as.data.frame(t(mean_Lambda))
names(loadings) <- paste0("F", 1:nrow(mean_Lambda))
loadings.wl <- loadings[-1,]
loadings.y <- loadings[1,]

loadings.wl$Wavelength <- as.numeric(rownames(loadings.wl))

loadings.wl.long <- tidyr::pivot_longer(loadings.wl, names_to = "Factor",
                                        cols = 1:(ncol(loadings.wl) - 1),
                                        values_to = "Loading")

loadings.wl.long$Factor <- factor(loadings.wl.long$Factor, levels = paste0("F", 1:nrow(mean_Lambda)))

# Some data for the spectral background:
conesdata <- read.csv("http://www.cvrl.org/database/data/cones/linss10e_5.csv")
names(conesdata) <- c("Wavelength", "Red", "Green", "Blue")
conesdata[is.na(conesdata)] <- 0
conesdata$colour <- rgb(conesdata$Red, conesdata$Green, conesdata$Blue, alpha = 0.6)
gradient <- t(conesdata$colour[conesdata$Wavelength >= 400 & conesdata$Wavelength <= 800])
g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# gradient <- t(photobiology::w_length2rgb(400:800))
# g <- rasterGrob(gradient, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)

# Plotting:
ggplot(data = loadings.wl.long, mapping = aes(x = Wavelength, y = Loading, color = Factor)) +
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  geom_line(linewidth = 1.5) +
  scale_color_manual(values = NatParksPalettes::natparks.pals("Acadia", 10)) +
  ylim(c(-1.1, 1.1)) +
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
  annotate("text", x = 757, y = -1.0, label = paste(("lambda(F1*', '* Y) * ' = ' *"), round(loadings.y["F1"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.75, label = paste(("lambda(F2 *', '* Y) * ' = ' *"), round(loadings.y["F2"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.90, label = paste(("lambda(F3*', '* Y) * ' = ' *"), round(loadings.y["F3"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.60, label = paste(("lambda(F4*', '* Y) * ' = ' *"), round(loadings.y["F4"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.45, label = paste(("lambda(F5*', '* Y) * ' = ' *"), round(loadings.y["F5"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = 0.30, label = paste(("lambda(F6*', '* Y) * ' = ' *"), round(loadings.y["F6"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = -0.85, label = paste(("lambda(F7*', '* Y) * ' = ' *"), round(loadings.y["F7"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = -0.70, label = paste(("lambda(F8*', '* Y) * ' = ' *"), round(loadings.y["F8"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = -0.55, label = paste(("lambda(F9*', '* Y) * ' = ' *"), round(loadings.y["F9"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  annotate("text", x = 757, y = -0.40, label = paste(("lambda(F10*', '* Y) * ' = ' *"), round(loadings.y["F10"], 2)),
           color = "white", parse = TRUE, size = 4, hjust = 0) +
  xlab("Wavelength (nm)") +
  annotate("text", x = 400, y = 1, label = "B", color = "white", parse = TRUE, size = 4, hjust = 0)

ggsave("plots/MegaLMM_hyper_single_date_M10.png", width = 24, height = 8.3, units = "cm")

# M = 5 traceplotting ==========================================================
rm(list = ls())
if (!("Lambda_samples" %in% ls() & "pred_samples" %in% ls())) {
  load("hyper/megalmm_hyper_single_date_arrays/LAMBDA_M5.RData")
  load("hyper/megalmm_hyper_single_date_arrays/PREDS_M5.RData")
}

dimnames(Lambda_samples)[[2]] <- paste0("F", 1:dim(Lambda_samples)[2])
pred_samples <- pred_samples[, , 1]
dimnames(pred_samples)[[2]] <- substr(dimnames(pred_samples)[[2]], 1, 7)

loadings.y <- as.data.frame(Lambda_samples[, , "Y"])
loadings.y$Iteration <- 1:nrow(loadings.y)
loadings.y <- tidyr::pivot_longer(loadings.y, cols = 1:(ncol(loadings.y) - 1), names_to = "Factor", values_to = "Loading")
loadings.y$Factor <- factor(loadings.y$Factor, levels = paste0("F", 1:dim(Lambda_samples)[2]))

# Plotting the significant Y loadings:
ggplot(loadings.y, aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Factor), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia")[3:7]) +
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
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/MegaLMM_hyper_single_date_traceplot_Y_M5.png", width = 24, height = 8, units = "cm")

# The important wavelengths around the "switching" point in gfBLUP:
(WL <- dimnames(Lambda_samples)[[3]][39:47])

factors <- paste0("F", 1:dim(Lambda_samples)[2])
niter <- dim(Lambda_samples)[1]
data <- as.data.frame(Lambda_samples[, factors[1], WL])
data$Factor <- rep(factors[1], niter)
data$Iteration <- 1:niter

for (f in factors[-1]) {
  subset <- as.data.frame(Lambda_samples[, f, WL])
  subset$Factor <- rep(f, niter)
  subset$Iteration <- 1:niter
  data <- rbind(data, subset)
}

data <- tidyr::pivot_longer(data, cols = 1:length(WL),
                            names_to = "Wavelength", values_to = "Loading")
data$Wavelength <- factor(data$Wavelength, levels = WL, labels = paste0(WL, " nm"))
data$Factor <- factor(data$Factor, levels = factors)

scaleFUN <- function(x) sprintf("%.3f", x)

p <- ggplot(data, aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Wavelength), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia")[1:9]) +
  theme_classic(base_size = 11) +
  facet_grid(rows = vars(Factor), scales = "free_y") +
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
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  ylab("Loading") +
  scale_y_continuous(labels=scaleFUN)

ggsave(plot = p, filename = "plots/MegaLMM_hyper_single_date_traceplot_WL_M5.png",
       width = 24, height = 31.3, units = "cm")

# M = 10 traceplotting =========================================================
rm(list = ls())
if (!("Lambda_samples" %in% ls() & "pred_samples" %in% ls())) {
  load("hyper/megalmm_hyper_single_date_arrays/LAMBDA_M10.RData")
  load("hyper/megalmm_hyper_single_date_arrays/PREDS_M10.RData")
}

dimnames(Lambda_samples)[[2]] <- paste0("F", 1:dim(Lambda_samples)[2])
pred_samples <- pred_samples[, , 1]
dimnames(pred_samples)[[2]] <- substr(dimnames(pred_samples)[[2]], 1, 7)

loadings.y <- as.data.frame(Lambda_samples[, , "Y"])
loadings.y$Iteration <- 1:nrow(loadings.y)
loadings.y <- tidyr::pivot_longer(loadings.y, cols = 1:(ncol(loadings.y) - 1), names_to = "Factor", values_to = "Loading")
loadings.y$Factor <- factor(loadings.y$Factor, levels = paste0("F", 1:dim(Lambda_samples)[2]))

# Plotting the significant Y loadings:
ggplot(loadings.y, aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Factor), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia", 10)) +
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
  guides(fill = guide_legend(override.aes = list(size = 5)))

ggsave("plots/MegaLMM_hyper_single_date_traceplot_Y_M10.png", width = 24, height = 8, units = "cm")

# The important wavelengths around the "switching" point in gfBLUP:
(WL <- dimnames(Lambda_samples)[[3]][39:47])

factors <- paste0("F", 1:dim(Lambda_samples)[2])
niter <- dim(Lambda_samples)[1]
data <- as.data.frame(Lambda_samples[, factors[1], WL])
data$Factor <- rep(factors[1], niter)
data$Iteration <- 1:niter

for (f in factors[-1]) {
  subset <- as.data.frame(Lambda_samples[, f, WL])
  subset$Factor <- rep(f, niter)
  subset$Iteration <- 1:niter
  data <- rbind(data, subset)
}

data <- tidyr::pivot_longer(data, cols = 1:length(WL),
                            names_to = "Wavelength", values_to = "Loading")
data$Wavelength <- factor(data$Wavelength, levels = WL, labels = paste0(WL, " nm"))
data$Factor <- factor(data$Factor, levels = factors)

scaleFUN <- function(x) sprintf("%.3f", x)

p1 <- ggplot(data[which(data$Factor %in% factors[1:5]),], aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Wavelength), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia")[1:9]) +
  theme_classic(base_size = 11) +
  facet_grid(rows = vars(Factor), scales = "free_y") +
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
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  ylab("Loading") +
  scale_y_continuous(labels=scaleFUN)

p2 <- ggplot(data[which(data$Factor %in% factors[6:10]),], aes(x = Iteration, y = Loading)) +
  geom_point(aes(fill = Wavelength), color = "black", pch = 21, stroke = 0.05, size = 1) +
  scale_fill_manual(values = NatParksPalettes::natparks.pals("Acadia")[1:9]) +
  theme_classic(base_size = 11) +
  facet_grid(rows = vars(Factor), scales = "free_y") +
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
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  ylab("Loading") +
  scale_y_continuous(labels=scaleFUN)

ggsave(plot = p1, filename = "plots/MegaLMM_hyper_single_date_traceplot_WL_M10_part1.png",
       width = 24, height = 31.3, units = "cm")

ggsave(plot = p2, filename = "plots/MegaLMM_hyper_single_date_traceplot_WL_M10_part2.png",
       width = 24, height = 31.3, units = "cm")



