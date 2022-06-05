############################# LOAD LIBRARIES ##################################
library(data.table)
library(foreach)
library(doParallel)
library(e1071)
library(xgboost)
library(randomForest)
library(caret)
library(mice)
library(Amelia)
library(DMwR)
library(pROC)
library(pryr)
library(rpart)
library(rpart.plot)
library(treeClust)

############################# GLOBAL OPTIONS ##################################
options(digits.secs=3)
Sys.setlocale("LC_TIME", "en_US")

############################# GLOBAL PATHS ##################################
path_root <- "/home/"
path_project <- file.path(path_root, "ygala", "UOC_TFM")
path_scripts <- file.path(path_root, "scripts", "R")
path_data <- file.path(path_project,  "output", "predictions")
input_path <- file.path(path_project, "data", "covid")

############################# SOURCE AUXILIARY SCRIPTS ##################################
source(file.path(path_scripts, "aux", "clean_data.R"))
source(file.path(path_scripts, "aux", "aux_functions.R"))

############################# GLOBAL PARAMETERS #############################
fixed_var <- c("dataset", "target")
min_transactions <- 3
max_n_var <- 300 # Maximum possible number of variables
samples_per_var <- 1 # n_variables*samples_per_var <= n_samples
min_n_categories <- 5 # Minimum number of categories per variable to set as maximum
min_size_others <- 3 # Minimun number of low frequency categories to apply grouping

############################# MAIN ######################################################D


############################# READ DATA #############################

### Predictions
dat_train <- fread(file.path(path_data, "X_train_augmented.csv"))
dat_train$pred <- as.numeric(t(fread(file.path(path_data, "predictions_train_augmented.csv"))))
dat_train$dataset <- "train"

dat_val <- fread(file.path(path_data, "X_val_augmented.csv"))
dat_val$pred <- as.numeric(t(fread(file.path(path_data, "predictions_val_augmented.csv"))))
dat_val$dataset <- "val"

dat_test <- fread(file.path(path_data, "X_test_augmented.csv"))
dat_test$pred <- as.numeric(t(fread(file.path(path_data,  "predictions_test_augmented.csv"))))
dat_test$dataset <- "test"

dat <- rbind(dat_train[, intersect(colnames(dat_train), colnames(dat_test)), with = F], dat_val[, intersect(colnames(dat_train), colnames(dat_test)), with = F], dat_test)

colnames(dat) <- tolower(colnames(dat))
dat <- dat[, c("dataset", "patientid", "filename", "pred", "survival"), with = F]

### Tabular data
tab_dat <- readRDS(file.path(input_path, "tabular_data.RDS"))
tab_dat$studytime <- as.integer(tab_dat$studytime)
dat <- merge(dat, tab_dat, by = c("patientid", "filename"))

############################# PRE-PROCESSING #############################
# Get target
target <- as.factor(dat$survival)

# Get days from admission, urg, icu
dat$image_date <- as.Date(as.character(dat$studydate), format = '%Y%m%d')
dat$fec_ing <- as.Date(dat$`f_ingreso/admission_d_ing/inpat`,format = '%d/%m/%Y')
dat$fec_urg <- as.Date(dat$`f_ingreso/admission_date_urg/emerg`,format = '%d/%m/%Y')
dat$fec_icu <- as.Date(dat$`f_entrada_uc/icu_date_in`,format = '%d/%m/%Y')
dat$days_img <- as.numeric(dat$image_date - dat$fec_ing)
dat$month_ing <- as.character(month(dat$fec_ing))
dat$weekday_ing <- weekdays(dat$fec_ing)
dat$days_urg <- as.numeric(dat$fec_urg - dat$fec_ing)
dat$month_urg <- as.character(month(dat$fec_urg))
dat$weekday_urg <- weekdays(dat$fec_urg)
dat$days_icu <- as.numeric(dat$fec_icu - dat$fec_ing)
dat$month_icu <- as.character(month(dat$fec_icu))
dat$weekday_icu <- weekdays(dat$fec_icu)

# Get hour of admission
dat$`hora/time_admision/admission_urg/emerg` <- as.character(substr(dat$`hora/time_admision/admission_urg/emerg`, 1, 1))

# Select vars
selected_vars <- setdiff(colnames(dat), c("patientid", "filename", "survival", "studytime", "studydate", "f_ingreso/admission_d_ing/inpat",
                                          "image_date", "fec_ing", "f_ingreso/admission_date_urg/emerg", "f_entrada_uc/icu_date_in",
                                          "fec_urg", "fec_icu"))
dat <- data.table(dat[, selected_vars, with = F], target = target)


# Castings
dat[, setdiff(colnames(dat)[sapply(dat, class) == "character"], fixed_var)] <- 
  data.table(sapply(dat[,setdiff(colnames(dat)[sapply(dat, class) == "character"], fixed_var), with = F], as.factor),
             stringsAsFactors = T)
dat[, setdiff(colnames(dat)[sapply(dat, class) == "integer"], fixed_var)] <- 
  data.table(sapply(dat[,setdiff(colnames(dat)[sapply(dat, class) == "integer"], fixed_var), with = F], as.numeric),
             stringsAsFactors = T)
# Clean data
dat_original <- dat
res <- clean_data(dat,
                  selected_variables = c(),
                  factor_datos = 0 ,
                  threshold_variables = 0,
                  over_factor = 0,
                  k = 3,
                  normalization = "mean", 
                  fill = "mean", 
                  anomaly_factor = 3,
                  threshold_variance = 0,
                  redundant_threshold = 0.99,
                  irrelevant_threshold = 0.001,
                  restricted_oversampling = F,
                  n_perc_varimp = 0,
                  seed=TRUE,
                  threshold_nas = 0.4,
                  fixed_var = fixed_var, 
                  target = dat$target)
dat <- res$dat

################################################## [1] GRID SEARCH ##################################################
dat_original <- dat




################################################## [1.1] LINEAR REGRESSION ##################################################
### Train
model <- glm(formula = "target ~ .", data = dat[dataset == "train", -setdiff(fixed_var, "target"), with = F], family = "binomial", control = list(maxit = 100))

### Predict
pred_train <- predict(model, newdata = dat[dataset == "train"])
pred_val<- predict(model, newdata = dat[dataset == "val"])

### Evaluation
auc_train <- AUC(target = dat[dataset == "train"]$target, prediction =  pred_train)
auc_val <- AUC(target = dat[dataset == "val"]$target, prediction =  pred_val)

### Results
grid_results_lr <- data.table(auc_train = auc_train,  auc_val = auc_val)
grid_results_lr
################################################## [1.2] RANDOM FOREST ##################################################

### Random Forest parameters grid definition
var_names <- c("ntree_values","mtry_values","nodesize_values")

# Initial grid
n_grid <- c(2, 2, 2)
lims <- list(ntree_values=c(10,500),mtry_values=c(1, 20),nodesize=c(3,50))
params_rf <- define_grid(var_names, n_grid, lims, previous_best = NA, cherkassky_point = NA, empirical_point = NA)
params_rf$ntree_values <- round(params_rf$ntree_values, 0)
params_rf$mtry_values <- round(params_rf$mtry_values, 0)
params_rf$nodesize_values <- round(params_rf$nodesize_values, 0)
params <- params_rf

### GRID
grid_results_rf <- foreach (ntree = params$ntree_values, .combine = rbind)%:%
  foreach (mtry = params$mtry_values, .combine = rbind)%:%
  foreach (nodesize = params$nodesize_values, .combine = rbind)%dopar%{
    
      # Log
      print(sprintf("Start of ntree = %s - mtry = %s - nodesize = %s", ntree, mtry, nodesize))
    
      # Train
      set.seed(14)
      model <- randomForest(x = dat[dataset == "train", -c("target", "dataset"), with = F], y = dat[dataset == "train"]$target,
                            ntree=ntree, mtry=mtry,nodesize = nodesize)
      
      # Predict
      pred_train <- predict(model, newdata = dat[dataset == "train"], type = "prob")[,2]
      pred_val<- predict(model, newdata = dat[dataset == "val"], type = "prob")[,2]
      
      # Evaluation
      auc_train <- AUC(target = dat[dataset == "train"]$target, prediction =  pred_train, quiet = TRUE)
      auc_val <- AUC(target = dat[dataset == "val"]$target, prediction =  pred_val, quiet = TRUE)
      
      # Log
      print(sprintf("End of ntree = %s - mtry = %s - nodesize = %s. AUC train = %s - AUC val = %s", 
                    ntree, mtry, nodesize, auc_train, auc_val))
      
      
      ### Results
      ret<-data.table(ntree = ntree, mtry = mtry, nodesize = nodesize,
                      auc_train = auc_train,  auc_val = auc_val)
      ret
    }
  


################################################## [1.3] XGBOOST ##################################################
### xgboost parameters grid definition
var_names<-c("eta_values","gamma_values","max_depth_values","min_child_weight_values",
             "subsample_values","colsample_bytree_values","num_parallel_tree_values","nrounds_values",
             "lambda_values", "alpha_values")

# Initial grid
n_grid<-c(1,1,1,1,1,1,1,1,1,1)
lims<-list(eta_values=c(10^-2, 10^0), gamma_values=c(10^-3, 10^3),max_depth_values=c(10, 100),min_child_weight_values=c(1, 5),
           subsample_values=c(0.5, 1),colsample_bytree_values=c(0.1, 1),num_parallel_tree_values=c(10, 100),
           nrounds_values=c(10, 30), lambda_values = c(0, 1), alpha_values = c(0, 1))

params_xgboost <- define_grid(var_names,n_grid,lims,previous_best = NA, empirical_point = NA)
params_xgboost$max_depth_values <- round(params_xgboost$max_depth_values)
params_xgboost$num_parallel_tree_values <- round(params_xgboost$num_parallel_tree_values)
params_xgboost$nrounds_values <- round(params_xgboost$nrounds_values)
params_xgboost$nthread <- detectCores()-1
params <- params_xgboost


grid_results_xgb <- data.table()
for (i in 1:length(params$eta_values)){
  for (j in 1:length(params$gamma_values)){
    for (k in 1:length(params$max_depth_values)){
      for (l in 1:length(params$min_child_weight_values)){
        for (m in 1:length(params$subsample_values)){
          for (n in 1:length(params$colsample_bytree_values)){
            for (o in 1:length(params$num_parallel_tree_values)){
              for (p in 1:length(params$nrounds_values)){
                for (q in 1:length(params$lambda_values)){
                  for (r in 1:length(params$alpha_values)){
                    eta <- params$eta_values[i]
                    gamma<-params$gamma_values[j]
                    max_depth<-params$max_depth_values[k]
                    min_child_weight<-params$min_child_weight_values[l]
                    subsample<-params$subsample_values[m]
                    colsample_bytree <- params$colsample_bytree_values[n]
                    num_parallel_tree <- params$num_parallel_tree_values[o]
                    nrounds <- params$nrounds_values[p]
                    lambda <- params$lambda_values[q]
                    alpha <- params$alpha_values[r]
                    
                    # Log
                    print(sprintf("Start of params=%s",c(eta,gamma,max_depth,min_child_weight,subsample,colsample_bytree,num_parallel_tree,nrounds, lambda, alpha)))
                    
                    # Train
                    set.seed(14)
                    flag_constant <- 0     
                    dat_train_xgboost <- xgb.DMatrix(data = as.matrix(dat[dataset == "train", -fixed_var, with = F]), 
                                                     label = as.numeric(dat[dataset == "train"]$target) - 1)
                    model <- xgboost(data = dat_train_xgboost, params=list(eta = eta, gamma = gamma,max_depth=max_depth,
                                                                             min_child_weight=min_child_weight,subsample=subsample,colsample_bytree=colsample_bytree,
                                                                             num_parallel_tree=num_parallel_tree, lambda = lambda, alpha = alpha, objective = "binary:logistic", eval_metric = "auc"),
                                       save_period=NULL,
                                       nthread = params$nthread, nrounds = nrounds, verbose = 0)
                      
                    # Predict
                    pred_train <- predict(model, newdata = as.matrix(dat[dataset == "train", -fixed_var, with = F]))
                    pred_val <- predict(model, newdata = as.matrix(dat[dataset == "val", -fixed_var, with = F]))
                      
                    if (length(unique(pred_train)) == 1){
                        print("Constant model!!!")
                        flag_constant <- 1
                        break
                    }
                    
                    # Evaluation
                    auc_train <- AUC(target = dat[dataset == "train"]$target, prediction =  pred_train, quiet = TRUE)
                    auc_val <- AUC(target = dat[dataset == "val"]$target, prediction =  pred_val, quiet = TRUE)
                    
                    # Log
                    print(paste0(paste(sprintf("End of params= %s. ",
                                  c(eta,gamma,max_depth,min_child_weight,subsample,colsample_bytree,num_parallel_tree,nrounds, lambda, alpha)), sep = "", collapse = ""),
                                 sprintf("AUC train = %s - AUC val = %s", auc_train, auc_val)))
                    
                    
                    # Results
                    ret <- data.table(eta = eta, gamma = gamma,max_depth=max_depth,
                                    min_child_weight=min_child_weight,subsample=subsample,colsample_bytree=colsample_bytree,
                                    num_parallel_tree=num_parallel_tree, nrounds = nrounds, 
                                    lambda = lambda, alpha = alpha, nthread = params$nthread,
                                    auc_train = auc_train, auc_val = auc_val)
                    grid_results_xgb <- rbind(grid_results_xgb, ret)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}


################################################## [1.4] LINEAR SVM ##################################################

### SVM parameters grid definition
var_names<-c("c_values")

# Initial grid
n_grid<-c(7)
lims<-list(c_values = 10^c(-3, 3))

params_linear_svm <- define_grid(var_names, n_grid, lims, previous_best = NA, cherkassky_point, empirical_point)
params <- params_linear_svm

### GRID
grid_results_linear_svm <- foreach (c = params$c_values, .combine = rbind)%dopar%{
    
    # Log
    print(sprintf("Start of c = %s", c))
    
    # Train
    set.seed(14)
    model <- svm(x = dat[dataset %in% c("train"), -fixed_var, with = F],
                 y = dat[dataset %in% c("train")]$target, 
                 kernel = "linear", tolerance = 0.00001, cachesize = 2000,
                 cost = c,
                 shrinking = TRUE, x.scale=F, y.scale=F, scale = FALSE,
                 probability = TRUE)
    
    # Predict
    pred_train <- predict(model, newdata=dat[dataset == "train",-fixed_var,with=F],probability = TRUE)
    pred_val <- predict(model, newdata=dat[dataset == "val",-fixed_var,with=F],probability = TRUE)
    pred_train <- attr(pred_train, "probabilities")[,1]
    pred_val <- attr(pred_val, "probabilities")[,1]
    
    # Evaluation
    auc_train <- AUC(target = dat[dataset == "train"]$target, prediction =  pred_train, quiet = TRUE)
    auc_val <- AUC(target = dat[dataset == "val"]$target, prediction =  pred_val, quiet = TRUE)
    
    # Log
    print(sprintf("End of c = %s. AUC train = %s - AUC val = %s", 
                  c, auc_train, auc_val))
    
    
    ### Results
    ret<-data.table(c = c, 
                    auc_train = auc_train,  auc_val = auc_val)
    ret
  }

################################################## [1.5] GAUSSIAN SVM ##################################################

### SVM parameters grid definition
var_names<-c("c_values", "g_values")
  
# Initial grid
n_grid<-c(1, 1)
lims<-list(c_values = 10^c(-3, 3), g_values = 10^c(-6, 0))


# Cherkassky point
ret <- cherkassky(dat[dataset %in% c("train")], as.numeric(dat[dataset %in% c("train")]$target)-1)
cherkassky_point <- data.table(c_values = ret$cost,
                               eps_values = ret$eps,
                               g_values = ret$gamma)

# Empirical point
features <- ncol(dat)
empirical_point <- data.table(c_values = 10^4,
                              eps_values = sd(as.numeric(dat[dataset %in% c("train")]$target)-1)/10,
                              g_values = 1/features)

params_svm <- define_grid(var_names, n_grid, lims, previous_best = NA, cherkassky_point, empirical_point)
params <- params_svm

### GRID
grid_results_svm <- foreach (c = params$c_values, .combine = rbind)%:%
  foreach (g = params$g_values, .combine = rbind)%dopar%{
    
    # Log
    print(sprintf("Start of c = %s - gamma = %s", c, g))
    
    # Train
    set.seed(14)
    model <- svm(x = dat[dataset %in% c("train"), -fixed_var, with = F],
                 y = dat[dataset %in% c("train")]$target, 
                 kernel = "radial", tolerance = 0.00001, cachesize = 2000,
                 gamma = g,  cost = c,
                 shrinking = TRUE, x.scale=F, y.scale=F, scale = FALSE,
                 probability = TRUE)
    
    # Predict
    pred_train <- predict(model, newdata=dat[dataset == "train",-fixed_var,with=F],probability = TRUE)
    pred_val <- predict(model, newdata=dat[dataset == "val",-fixed_var,with=F],probability = TRUE)
    pred_train <- attr(pred_train, "probabilities")[,1]
    pred_val <- attr(pred_val, "probabilities")[,1]
    
    # Evaluation
    auc_train <- AUC(target = dat[dataset == "train"]$target, prediction =  pred_train, quiet = TRUE)
    auc_val <- AUC(target = dat[dataset == "val"]$target, prediction =  pred_val, quiet = TRUE)
    
    # Log
    print(sprintf("End of c = %s - g = %s. AUC train = %s - AUC val = %s", 
                  c, g, auc_train, auc_val))
    
    
    ### Results
    ret<-data.table(c = c, g = g, 
                    auc_train = auc_train,  auc_val = auc_val)
    ret
  }


models <- c("Gaussian SVM",
            "Random Forest",
            "XGBoost",
            "Linear SVM",
            "Linear Regression")
flag_best_model <- which.max(c(max(grid_results_svm$auc_val,
                               max(grid_results_xgb$auc_val),
                               max(grid_results_rf$auc_val),
                               max(grid_results_linear_svm$auc_val,
                               max(grid_results_lr$auc_val)))))
flag_best_model <- models[flag_best_model]
auc_val <- max(c(max(grid_results_lr$auc_val),
                 max(grid_results_rf$auc_val),
                 max(grid_results_xgb$auc_val),
                 max(grid_results_linear_svm$auc_val),
                 max(grid_results_svm$auc_val)))

################################################## [2] TRAIN ##################################################
# Combine train and val
# dat[dataset == "val"]$dataset <- "train"

if (flag_best_model == "Linear Regression"){
  ### Train
  model <- glm(formula = "target ~ .", data = dat[dataset == "train", -setdiff(fixed_var, "target"), with = F], family = "binomial", control = list(maxit = 100))
  
  ### Predict
  pred_train <- predict(model, newdata = dat[dataset == "train"], type= "response")
  pred_test<- predict(model, newdata = dat[dataset == "test"], type= "response")
} else if (flag_best_model == "Random Forest"){
  # Best grid point
  grid_results_rf <- grid_results_rf[order(-auc_val, -auc_train)]
  best <- grid_results_rf[1]
  
  # Train
  set.seed(14)
  model <- randomForest(x = dat[dataset == "train", -c("target", "dataset"), with = F], y = dat[dataset == "train"]$target,
                        ntree=best$ntree, mtry=best$mtry,nodesize = best$nodesize)
  
  # Predict
  pred_train <- predict(model, newdata = dat[dataset == "train"], type = "prob")[,2]
  pred_test<- predict(model, newdata = dat[dataset == "test"], type = "prob")[,2]
} else if (flag_best_model == "XGBoost"){
  # Best grid point
  grid_results_xgb <- grid_results_xgb[order(-auc_val, -auc_train)]
  best <- grid_results_xgb[1]
  
  # Train
  set.seed(14)
  flag_constant <- 0     
  dat_train_xgboost <- xgb.DMatrix(data = as.matrix(dat[dataset == "train", -fixed_var, with = F]), 
                                   label = as.numeric(dat[dataset == "train"]$target) - 1)
  model <- xgboost(data = dat_train_xgboost, params=list(eta = best$eta, gamma = best$gamma,max_depth=best$max_depth,
                                                         min_child_weight=best$min_child_weight,subsample=best$subsample,colsample_bytree=best$colsample_bytree,
                                                         num_parallel_tree=best$num_parallel_tree, lambda = best$lambda, alpha = best$alpha, objective = "binary:logistic", eval_metric = "auc"),
                   save_period=NULL,
                   nthread = best$nthread, nrounds = best$nrounds, verbose = 0)
  
  # Predict
  pred_train <- predict(model, newdata = as.matrix(dat[dataset == "train", -fixed_var, with = F]))
  pred_test <- predict(model, newdata = as.matrix(dat[dataset == "test", -fixed_var, with = F]))
} else if (flag_best_model == "Linear SVM"){
  # Best grid point
  grid_results_linear_svm <- grid_results_linear_svm[order(-auc_val, -auc_train)]
  best <- grid_results_linear_svm[1]
  
  # Train
  set.seed(14)
  model <- svm(x = dat[dataset %in% c("train"), -fixed_var, with = F],
               y = dat[dataset %in% c("train")]$target, 
               kernel = "linear", tolerance = 0.00001, cachesize = 2000,
               cost = best$c,
               shrinking = TRUE, x.scale=F, y.scale=F, scale = FALSE,
               probability = TRUE)
  
  # Predict
  pred_train <- predict(model, newdata=dat[dataset == "train",-fixed_var,with=F],probability = TRUE)
  pred_test <- predict(model, newdata=dat[dataset == "test",-fixed_var,with=F],probability = TRUE)
  pred_train <- attr(pred_train, "probabilities")[,1]
  pred_test <- attr(pred_test, "probabilities")[,1]
} else if (flag_best_model == "Gaussian SVM"){
  # Best grid point
  grid_results_svm <- grid_results_svm[order(-auc_val, -auc_train)]
  best <- grid_results_svm[1]
  # Train
  set.seed(14)
  model <- svm(x = dat[dataset %in% c("train"), -fixed_var, with = F],
               y = dat[dataset %in% c("train")]$target, 
               kernel = "radial", tolerance = 0.00001, cachesize = 2000,
               gamma = best$g,  cost = best$c,
               shrinking = TRUE, x.scale=F, y.scale=F, scale = FALSE,
               probability = TRUE)
  
  # Predict
  pred_train <- predict(model, newdata=dat[dataset == "train",-fixed_var,with=F],probability = TRUE)
  pred_test <- predict(model, newdata=dat[dataset == "test",-fixed_var,with=F],probability = TRUE)
  pred_train <- attr(pred_train, "probabilities")[,1]
  pred_test <- attr(pred_test, "probabilities")[,1]
}



################################################## [3] EVALUATION ##################################################
auc_train <- AUC(target = dat[dataset == "train"]$target, prediction =  pred_train, quiet = TRUE)
auc_val < -auc_val
auc_test <- AUC(target = dat[dataset == "test"]$target, prediction =  pred_test, quiet = TRUE)
# Log
print(sprintf("AUC train = %s - AUC val = %s - AUC test = %s ", 
              auc_train, auc_val, auc_test))

