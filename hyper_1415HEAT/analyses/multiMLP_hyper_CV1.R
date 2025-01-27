###############################################################################
### Preliminaries #############################################################
###############################################################################

# Will take about 6 hours for all 250 datasets on the CPU.

# Setting CV:
CV <- "CV1"
prep <- "splines"
# Setting seed:
# set.seed(1997)

# Loading libraries:
library(rlist)
library(tictoc)
library(keras)
library(gfBLUP)
library(tensorflow)

# Setting seeds:
tensorflow::set_random_seed(1997, disable_gpu = TRUE)

# Setting working directory:
wd <- getwd()
setwd(wd)

# Load marker data:
load("genotypes/M_hyper.RData")
M <- M_hyper; rm(M_hyper)
M <- as.data.frame(M)
M <- cbind(data.frame(G = rownames(M)), M)

###############################################################################
### Defining all 8 architectures ##############################################
###############################################################################

# Setting hyperparameter values and creating list of all combinations:
L <- c(4, 8) # Number of layers
N <- c(16, 64) # Number of neurons per layer
D <- c(0.1, 0.3) # Dropout regularization rate

hypers <- vector("list", length(L) * length(N) * length(D))
i <- 1
for (l in 1:length(L)) {
  for (n in 1:length(N)) {
    for (d in 1:length(D)) {
      hypers[[i]] <- list(layers = L[l], neurons = N[n], dropout = D[d])
      i <- i + 1
    }
  }
}

# Function to define the required architectures:
define_model <- function(L, N, D, input_shape, output_shape) {
  if (L == 1) {
    model <- keras_model_sequential()
    
    model %>%
      layer_dense(units = N, activation = "relu", input_shape = input_shape) %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = output_shape)
    
    name <- sprintf("multiMLP_L1_N%d_D%.1f_CV1", N, D)
    
  } else if (L == 2) {
    model <- keras_model_sequential()
    
    model %>%
      layer_dense(units = N, activation = "relu", input_shape = input_shape) %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = output_shape)
    
    name <- sprintf("multiMLP_L2_N%d_D%.1f_CV1", N, D)
  } else if (L == 4) {
    model <- keras_model_sequential()
    
    model %>%
      layer_dense(units = N, activation = "relu", input_shape = input_shape) %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = output_shape)
    
    name <- sprintf("multiMLP_L4_N%d_D%.1f_CV1", N, D)
  } else if (L == 8) {
    model <- keras_model_sequential()
    
    model %>%
      layer_dense(units = N, activation = "relu", input_shape = input_shape) %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = N, activation = "relu") %>%
      layer_dropout(rate = D) %>%
      layer_dense(units = output_shape)
    
    name <- sprintf("multiMLP_L8_N%d_D%.1f_CV1", N, D)
  }
  return(list(model = model, name = name))
}

###############################################################################
### Running on datasets #######################################################
###############################################################################

# Hyperspectral datasets to run:
first <- 1
last <- 250
n.datasets <- length(first:last)
      
# Creating test accuracy storage:
acc <- numeric(n.datasets)
archs <- character(n.datasets)
histories <- vector("list", n.datasets)

tic(sprintf("%s multiMLP", CV))
cat(sprintf("%s multiMLP\n\n", CV))
run <- first
for (run in first:last) {
  
  cat(sprintf("### Running on dataset %d / %d (dataset %d) ###\n\n", run, n.datasets, run))
  
  # Loading simulated datasets:
  cat(sprintf("Loading dataset %d...\n", run))
  datalist <- list.load(sprintf("hyper_1415HEAT/datasets/%s/hyper_dataset_%d.RData", prep, run))
  pred.target <- datalist$pred.target
  train.set <- datalist$train.set
  test.set <- datalist$test.set
  
  # Storing training and test data:
  d.train <- droplevels(na.omit(datalist$data))
  # d.test <- droplevels(datalist$data[which(is.na(datalist$data$Y)),])
  
  # Redundancy filtering:
  sec <- names(d.train[2:(ncol(d.train) - 1)])
  foc <- names(d.train)[ncol(d.train)]
  temp <- gfBLUP::redundancyFilter(data = d.train[c("G", sec)], tau = 0.95, verbose = TRUE)
  d.train <- cbind(temp$data.RF, d.train[foc])
  # d.test <- d.test[colnames(d.train)]
  sec.RF <- names(d.train[2:(ncol(d.train) - 1)])
  
  # Going to BLUEs:
  d.train <- gfBLUP:::genotypeMeans(d.train)
  # d.test <- gfBLUP:::genotypeMeans(d.test)
  
  # Rescaling now we only have means:
  d.train[, 2:ncol(d.train)] <- sapply(d.train[, 2:ncol(d.train)], scale)
  # d.test[, 2:ncol(d.test)] <- sapply(d.test[, 2:ncol(d.test)], scale)
  
  d.train[, 2:ncol(d.train)] <- sapply(d.train[, 2:ncol(d.train)], as.numeric)
  # d.test[, 2:ncol(d.test)] <- sapply(d.test[, 2:ncol(d.test)], as.numeric)
  
  # Subsetting marker data to train and test set and scaling:
  M.train <- M[which(M$G %in% train.set),]
  M.test <- M[which(M$G %in% test.set),]
  
  M.train[, 2:ncol(M.train)] <- sapply(M.train[, 2:ncol(M.train)], scale)
  M.test[, 2:ncol(M.test)] <- sapply(M.test[, 2:ncol(M.test)], scale)
  
  M.train[, 2:ncol(M.train)] <- sapply(M.train[, 2:ncol(M.train)], as.numeric)
  M.test[, 2:ncol(M.test)] <- sapply(M.test[, 2:ncol(M.test)], as.numeric)
  
  # We need to check whether this specific train-test division led to markers within
  # the training or test set that had 0 variance and therefore became NA after scaling:
  discard <- which(is.na(colMeans(M.train[, -1])))
  discard <- unique(names(c(discard, which(is.na(colMeans(M.test[, -1]))))))
  M.train <- M.train[, setdiff(names(M.train), discard)]
  M.test <- M.test[, setdiff(names(M.test), discard)]
  
  # Creating folds for 5-fold hyper-parameter tuning:
  tic(sprintf("Dataset %d", run))
  ngeno <- length(d.train$G)
  folds <- gfBLUP::createFolds(genos = d.train$G, folds = 5)
  
  # Creating temporary tuning results storage:
  tuning.results <- data.frame()
  
  # Doing 5-fold CV tuning:
  fold <- 1
  for (fold in 1:length(folds)) {
    
    # Making vectors of tuning and evaluation genotypes:
    eval.set <- d.train$G[d.train$G %in% folds[[fold]]]
    tune.set <- setdiff(d.train$G, eval.set)
    
    # Creating tuning and evaluation marker input layer data (CV2 so marker and secondary data in input layer):
    x_tune <- M.train[which(M.train$G %in% tune.set),]
    x_eval <- M.train[which(M.train$G %in% eval.set),]
    
    # Creating tuning and evaluation output layer data (CV2 so only yield data in output layer).
    # Yield data is the last column:
    y_tune <- d.train[which(d.train$G %in% tune.set), ]
    
    # Matching x_tune to y_tune:
    x_tune <- x_tune[match(y_tune$G, x_tune$G),]
    
    # Storing evaluation prediction targets and matching (note the selection of only yield for y_eval!):
    y_eval <- d.train[which(d.train$G %in% eval.set), c(1, ncol(d.train))]
    x_eval <- x_eval[match(y_eval$G, x_eval$G),]
    
    rownames(y_eval) <- rownames(x_eval)
    rownames(y_tune) <- rownames(x_tune)
    
    # Scaling again:
    y_eval[, 2] <- scale(y_eval[, 2])
    y_tune[, 2:ncol(y_tune)] <- sapply(y_tune[, 2:ncol(y_tune)], scale)
    
    x_eval[, 2:ncol(x_eval)] <- sapply(x_eval[, 2:ncol(x_eval)], scale)
    x_tune[, 2:ncol(x_tune)] <- sapply(x_tune[, 2:ncol(x_tune)], scale)
    
    y_eval[, 2] <- as.numeric(y_eval[, 2])
    y_tune[, 2:ncol(y_tune)] <- sapply(y_tune[, 2:ncol(y_tune)], as.numeric)
    
    x_eval[, 2:ncol(x_eval)] <- sapply(x_eval[, 2:ncol(x_eval)], as.numeric)
    x_tune[, 2:ncol(x_tune)] <- sapply(x_tune[, 2:ncol(x_tune)], as.numeric)
    
    # We need to check for 0-variance markers in the evaluation or tuning set again:
    discard <- which(is.na(colMeans(x_eval[, -1])))
    discard <- unique(names(c(discard, which(is.na(colMeans(x_tune[, -1]))))))
    x_eval <- x_eval[, setdiff(names(x_eval), discard)]
    x_tune <- x_tune[, setdiff(names(x_tune), discard)]
    
    # Creating temporary per-architecture result storage:
    temp <- numeric(length(hypers))
    
    #Evaluating accuracy of all architectures for the current fold:
    arch <- 1
    for (arch in 1:length(hypers)) {
      
      cat(sprintf("Evaluating architecture %d / %d on fold %d / 5 for dataset %d / %d (%d)...\n",
                  arch, length(hypers), fold, run, n.datasets, run))
      
      # Defining and compiling architecture (note that the output_shape = ncol(y_tune) due to
      # the secondary features being in the output layer together with yield).
      # The input shape is ncol(x_tune) - 1 (1 because "G" is omitted):
      
      model.def <- define_model(L = hypers[[arch]]$layers,
                                N = hypers[[arch]]$neurons,
                                D = hypers[[arch]]$dropout,
                                input_shape = ncol(x_tune) - 1,
                                output_shape = ncol(y_tune) - 1)
      
      model <- model.def$model
      model.name <- model.def$name
      
      model %>% compile(loss = "mse",
                        optimizer = optimizer_rmsprop(learning_rate = 0.0005),
                        metrics = list())
      
      # Training network:
      history.arch <- model %>% fit(as.matrix(x_tune[, -1]),
                                    as.matrix(y_tune[, -1]),
                                    epochs = 25, batch_size = 128,
                                    validation_split = 0.1, verbose = 0,
                                    callbacks = callback_early_stopping(monitor = "val_loss",
                                                                        patience = 10,
                                                                        mode = "auto",
                                                                        restore_best_weights = FALSE))
      
      # Making evaluation predictions and calculating accuracy:
      preds <- model %>% predict(as.matrix(x_eval[, -1]))
      temp[arch] <- cor(preds[, ncol(preds)], y_eval$Y)
      names(temp)[arch] <- model.name
      
      # Clearing Keras backend:
      k_clear_session()
    }
    # Adding results for the current fold to the overall results:
    tuning.results <- rbind(tuning.results, temp)
    names(tuning.results) <- names(temp)
  }
  mean.accs <- sapply(tuning.results, mean)
  
  # Setting any architecture with failed fits to 0:
  mean.accs[is.na(mean.accs)] <- 0
  
  best.arch <- which(mean.accs == max(mean.accs))
  
  if (mean(mean.accs) != 0) {
    
    cat(sprintf("Best architecture is %s with 5-fold mean accuracy of %.2f\n", names(best.arch), mean.accs[best.arch]))
    
    # Redefining and recompiling the best architecture:
    model.def <- define_model(L = hypers[[best.arch]]$layers,
                              N = hypers[[best.arch]]$neurons,
                              D = hypers[[best.arch]]$dropout,
                              input_shape = ncol(M.train) - 1,
                              output_shape = ncol(d.train) - 1)
    
    model <- model.def$model
    model.name <- model.def$name
    
    model %>% compile(loss = "mse",
                      optimizer = optimizer_rmsprop(learning_rate = 0.0005),
                      metrics = list())
    
    # matching data just to be sure:
    M.train <- M.train[match(d.train$G, M.train$G),]
    M.test <- M.test[match(pred.target$G, M.test$G),]
    
    cat(sprintf("Training and evaluating hypertuned CV1 multiMLP for dataset %d / %d (%d)...\n\n", run, n.datasets, run))
    
    # Training the final hypertuned network using early stopping:
    history <- model %>% fit(as.matrix(M.train[, -1]),
                             as.matrix(d.train[, -1]),
                             epochs = 100, batch_size = 128,
                             validation_split = 0.1, verbose = 0,
                             callbacks = callback_early_stopping(monitor = "val_loss",
                                                                 patience = 20,
                                                                 mode = "auto",
                                                                 restore_best_weights = FALSE))
    toc(log = TRUE)
    
    histories[[run]] <- history
    
    # Making predictions and evaluating accuracy:
    preds <- model %>% predict(as.matrix(M.test[, -1]))
    acc[run] <- cor(preds[, ncol(preds)], pred.target$pred.target)
    
    # Storing best architecture:
    archs[run] <- model.name
    
    # Clearing Keras backend:
    k_clear_session()
    
  } else if (mean(mean.accs) == 0) {
    
    cat("All fits failed!\n\n")
    acc[run] <- NA
    archs[run] <- NA
    
    # Clearing Keras backend:
    k_clear_session()
  }
}
toc()

# Retrieve computational times:
tictoc.logs <- tic.log(format = FALSE)
comptimes <- unlist(lapply(tictoc.logs, function(x) x$toc - x$tic))
tic.clearlog()

# Collect results:
results <- data.frame(acc = acc,
                      comptimes = comptimes,
                      archs = archs)

# Making correct CV label:
if (CV == "CV1") {
  lab <- "a"
} else if (CV == "CV2") {
  lab <- "b"
}

# Export results:
write.csv(results, sprintf("hyper_1415HEAT/results/%s/8%s_hyper_results_multiMLP_%s.csv",
                           prep, lab, CV))

list.save(histories, sprintf("hyper_1415HEAT/results/%s/8%s_hyper_results_multiMLP_%s.RData",
                             prep, lab, CV))


