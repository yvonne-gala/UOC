############################### CLEAN DATA ###############################################
clean_data <- function (dat,
                        selected_variables,factor_datos,threshold_variables,over_factor,k,normalization, fill, anomaly_factor,threshold_variance,
                        redundant_threshold,irrelevant_threshold, n_perc_varimp,
                        seed=TRUE,threshold_nas,
                        fixed_var, diag_vars, proc_vars,
                        target, country){
  
  
  ############################### [1] SELECT VARIABLES #####################################
  res <- select_variables(dat, selected_variables, fixed_var);
  dat <- res$dat;
  
  # Free memory
  rm(res);
  gc();
  
  ################################ [2] OVERSAMPLING ###########################
  if (over_factor!=0){
    ################################ [2.1] Normalization (only numeric) ###########################
    # Get numerical data
    numerical<-dat[, which(sapply(dat, class) == "numeric"), with=F];
    
    # Get mean and variance
    m_smote<-apply(numerical[dat$dataset=="train"],2, get_statistic, mean, na.rm=TRUE);
    s_smote<-apply(numerical[dat$dataset=="train"],2, get_statistic, sd, na.rm=TRUE);
    
    # Apply normalization
    numerical<-data.table(foreach (v = names(s_smote),.combine=cbind)%dopar%{
      m_v <- m_smote[names(m_smote)==v]
      s_v <- s_smote[names(s_smote)==v];
      s_v[s_v==0] <- 1;
      if (normalization=="tip") {
        ret<-(numerical[,v,with=F]-m_v)/s_v;
      } else {
        ret<-numerical[,v,with=F]-m_v;
      }
      ret;
    });
    
    
    ################################ [2.2] OVERSAMPLING (SMOTE - CN) ###########################
    
    # Data formatting
    dat[, which(sapply(dat, class) == "numeric")] <- numerical;
    
    # Set train dataset used for SMOTE
    smote_dat <- dat[dataset=="train"];
    
    # Remove unnecessary variables
    smote_dat<-smote_dat[,-setdiff(fixed_var, "target"), with=F];
    
    # Apply oversampling
    if (over_factor!=0){
      smote_dat_ori<-smote_dat;
      smote_dat <- data.table(SMOTE(target ~ ., smote_dat, perc.over = over_factor, k=k));
      
      # Keep only new samples created
      smote_dat_ori$type<-"old";
      smote_dat$type<-"new";
      comb<-rbind(smote_dat_ori,smote_dat);
      rep<-duplicated(comb[,-"type",with=F]);
      new<-comb[!rep & type=="new"];
      new<-new[,-"type",with=F];
      
      # Mark as new samples
      new$dataset<-"train";
    } 
    
    # Add new samples
    dat <- rbind(dat,new);
    
    # Free memory
    rm(new)
    ;
    gc();
    
    ################################ [2.3] RECOVER ORIGINAL INFORMATION ###########################
    
    # Unapply normalization
    numerical<-dat[, which(sapply(dat, class) == "numeric"), with=F];
    numerical<-data.table(foreach (v = names(s_smote),.combine=cbind)%dopar%{
      m_v<-m_smote[names(m_smote)==v]
      s_v<-s_smote[names(s_smote)==v]
      if (normalization=="tip") {
        ret<-(numerical[,v,with=F]*s_v)+m_v;
      } else {
        ret<-numerical[,v,with=F]+m_v;
      }
      ret;
    })
    
    # Paste all variables
    dat[, which(sapply(dat, class) == "numeric")] <- numerical;
  }
  
  # Free memory
  gc();
  
  ################################# [3] UNDERSAMPLING (RANDOM) ###########################
  ### UNDERSAMPLING (RANDOM)
  if (factor_datos != 0){
    dat <- rbind(dat[dataset!="train"], dat[dataset=="train" & target==1],
                 dat[dataset=="train" & target==0][sample(1:nrow(dat[dataset=="train" & target==0]), factor_datos*nrow(dat[dataset=="train" & target==1]))])
  }
  
  # Free memory
  gc();
  
  
  ############################### [4] REMOVE CONSTANT VARIABLES #####################################
  
  constant_var <- constant_variables(dat[dataset=="train"]);
  constant_var <- setdiff(constant_var,fixed_var);
  
  if (length(constant_var) > 0){
    dat <- dat[,-constant_var,with=F];
  }
  
  # Free memory
  gc();
  
  ############################### [5] REMOVE NOT INFORMED VARIABLES ############################## 
  
  not_informed_var <- not_informed(dat[dataset=="train", 
                                       which(sapply(dat, class) == "numeric"), 
                                       with=F],
                                   threshold_nas);
  not_informed_var<-setdiff(not_informed_var,fixed_var);
  
  if (length(not_informed_var) > 0){
    dat<-dat[,-not_informed_var,with=F];
  }
  
  # Free memory
  gc();
  
  
  ############################### [6] FILL MISSING VALUES #####################################
  # Fill numerical missing values
  dat <- fill_values(dat, fill);
  
  # Free memory
  gc();
  
  ################################ [7] UNTOUCHABLE CATEGORIES ###########################
  
  untouchable_values<-data.table(t(sapply(setdiff(colnames(dat)[which(sapply(dat, class) == "factor")], fixed_var), 
                                          function(x){untouchable_values_f(x,dat[dataset=="train"][[x]], as.numeric(dat[dataset == "train"]$target) - 1,
                                                                           min_transactions, anomaly_factor)})));
  
  # Free memory
  gc();
  ################################# [8] GROUP CATEGORIES ###########################
  # Get minimum total number of categories to apply grouping
  max_n_categorical_variables <- floor(min(c(sum(dat$dataset=="train")/samples_per_var, max_n_var)) - ncol(dat[,which(sapply(dat, class) == "numeric"),with=F]));
  max_n_categories <- max(min_n_categories, round(max_n_categorical_variables/length(setdiff(colnames(dat)[which(sapply(dat, class) == "factor")], fixed_var)) - 1));
  
  # Drop empty levels
  dat[, colnames(dat)[sapply(dat, class) == "factor"]] <- data.table(sapply(dat[, colnames(dat)[sapply(dat, class) == "factor"], with = F], droplevels),
                                                                     stringsAsFactors = T)
  
  # Group categories
  res <- group_categories(setdiff(colnames(dat)[which(sapply(dat, class) == "factor")], fixed_var),
                          dat,
                          dat[dataset=="train"],
                          threshold_variables, min_size_others, max_n_categories, untouchable_values,
                          diag_vars, proc_vars);
  dat<-res$dat;
  grouped_values_list<-res$grouped_values_list;
  
  # Drop empty levels
  dat[, colnames(dat)[sapply(dat, class) == "factor"]] <- data.table(sapply(dat[, colnames(dat)[sapply(dat, class) == "factor"], with = F], droplevels),
                                                                     stringsAsFactors = T)
  
  # Free memory
  gc();
  
  
  ################################# [9] ONE HOT SECONDARY DIAGNOSIS/PROCEDURES ###########################
  # Format secondary diagnosis and procedures
  diag_vars <- intersect(diag_vars, colnames(dat))
  proc_vars <- intersect(proc_vars, colnames(dat))
  dat <- one_hit_diag_proc(dat, diag_vars[-1], proc_vars[-1])
  
  # Free memory
  gc();
  
  
  ################################# [10] ONE HOT ENCODING ###########################
  # one hot encoding
  dat <- one_hot_encoding(dat,
                          categorical = setdiff(colnames(dat)[which(sapply(dat, class) == "factor")], fixed_var),
                          no_categorical = setdiff(colnames(dat)[which(sapply(dat, class) == "numeric")], fixed_var),
                          fixed_var);
  # Free memory
  gc();
  
  # Remove cheating variables result of applying One-hot encoding
  dat <- remove_cheating_variables_afteroh(dat, target, country = country);
  
  # Free memory
  gc();
  
  ################################# [11] TIPIFY (Categorical and numerical variables) ###########################
  # Get mean and standard deviation for train data (excluding target) 
  m_tip <- sapply(dat[dataset=="train",-fixed_var,with=F], get_statistic, mean);
  s_tip <- sapply(dat[dataset=="train",-fixed_var,with=F], get_statistic, sd);
  
  # Free memory
  gc();
  
  # Tipify data
  ret<-tipify(dat[,-fixed_var,with=F],m_tip,s_tip);
  
  # Free memory
  gc();
  
  dat<-data.table(dat[,fixed_var,with=F],ret$dat);
  m_tip<-ret$m;
  s_tip<-ret$s;
  
  # Free memory
  rm(ret);
  gc();
  
  ################################# [12] REMOVE IRRELEVANT ###########################
  if (irrelevant_threshold > 0 | redundant_threshold > 0 ){
    correlations<-abs(cor(data.table(target=as.numeric(dat[dataset=="train"]$target),dat[dataset=="train",-fixed_var,with=F])));
  }
  
  # Free memory
  gc();
  
  
  if (irrelevant_threshold > 0){
    irrelevant_variables<-remove_irrelevant(correlations,irrelevant_threshold);
    dat<-dat[,setdiff(colnames(dat),irrelevant_variables),with=F];
  }
  
  # Free memory
  gc();
  
  ################################# [13] REMOVE REDUNDANT ###########################
  # Redundant variables
  if (redundant_threshold > 0 ){
    redundant_variables<-remove_redundant(correlations,redundant_threshold);
    dat<-dat[,setdiff(colnames(dat),redundant_variables),with=F];
  }
  
  # Free memory
  gc();
  
  ################################# [14] VARIABLE IMPORTANCE ###########################
  # Variable importance
  res <- select_important(dat[dataset=="train"], n_perc_varimp);
  varimp <- res$varimp;
  if (n_perc_varimp > 0){
    important_variables <- res$selected;
    dat<-dat[,c(fixed_var,important_variables),with=F];
  }
  
  
  # Free memory
  gc();
  
  ################################# SAVE FINAL VARIABLES ###########################
  selected_vars <- setdiff(colnames(dat),fixed_var);
  
  ################################# [15] PCA (TRAIN) ###########################
  if(threshold_variance > 0){
    ret<-pca(dat,dat[dataset=="train"],threshold_variance, fixed_var);
    dat<-ret$dat;
    m_tip_pca <- ret$m_tip;
    s_tip_pca <- ret$s_tip;
    prin_comp<-ret$prin_comp;
  } else{ # No PCA
    prin_comp<-NA;
    m_tip_pca <- NA;
    s_tip_pca <- NA;
  }
  
  # Free memory
  gc();
  
  
  ############################ [16] RETURN RESULTS ###############################################
  return(list(dat = dat, m_tip = m_tip, s_tip = s_tip, m_tip_pca = m_tip_pca, s_tip_pca = s_tip_pca,
              selected_vars = selected_vars, prin_comp = prin_comp,
              grouped_values_list = grouped_values_list,
              varimp = varimp));
}

############################### CLEAN DATA ONLINE ###############################################
clean_data_online <- function (dat, best_grid, ret_clean, fixed_var){
  
  ############################### [1] SELECT VARIABLES #####################################
  res <- select_variables(dat, best_grid$selected_variables[[1]], fixed_var);
  dat <- res$dat;
  
  # Free memory
  rm(res);
  gc();
  
  ############################### [2] FILL MISSING VALUES #####################################
  # Fill numerical missing values
  dat <- fill_values(dat, fill);
  
  # Cast categorical missing values to "NA"
  index_NA <- colnames(dat)[sapply(dat, function(x){sum(is.na(x))}) > 0 & sapply(dat, class) == "factor"];
  dat[, index_NA] <- data.table(sapply(dat[, index_NA, with = F],
                                       missing_to_NA_category), stringsAsFactors = TRUE);
  
  
  ################################# [3] GROUP CATEGORIES ###########################
  dat[,which(colnames(dat) %in% setdiff(colnames(dat)[which(sapply(dat, class) == "factor")], fixed_var))]<-data.table(sapply(setdiff(colnames(dat)[which(sapply(dat, class) == "factor")], fixed_var), 
                                                                                                                              function(x){group_values(x, unlist(ret_clean$grouped_values_list[name==x]$grouped_values),dat[,x,with=F])}), stringsAsFactors = TRUE);
  ################################# [4] ONE HOT ENCODING ###########################
  # one hot encoding
  dat <- one_hot_encoding(dat, 
                          categorical = setdiff(colnames(dat)[which(sapply(dat, class) == "factor")], fixed_var),
                          no_categorical = setdiff(colnames(dat)[which(sapply(dat, class) == "numeric")], fixed_var),
                          fixed_var = intersect(colnames(dat), fixed_var));
  
  # Free memory
  gc();
  
  
  
  ################################# 54] FIT TO MODEL ###########################
  ret <- fit_to_model(dat, ret_clean);
  dat <- data.table(dat[, fixed_var, with = F], ret);
  
  ################################# [6] TIPIFY (Categorical and numerical variables) ###########################
  common_names <- intersect(names(ret_clean$m_tip), colnames(dat));
  ret <- tipify(dat[,common_names,with=F],ret_clean$m_tip[common_names],
                ret_clean$s_tip[common_names]);
  dat <-  data.table(dat[, fixed_var, with = F], ret$dat);
  
  
  ################################# [7] PCA (TRAIN) ###########################
  if(best_grid$threshold_variance > 0){
    # choose number of variables
    prop_varex <- summary(ret_clean$prin_comp)$importance[3,];
    number_of_variables <- min(which(prop_varex > best_grid$threshold_variance));
    dat <- data.table(dat[, fixed_var, with = F],predict(ret_clean$prin_comp, newdata = dat)[,1:number_of_variables,with=F]);
    
    #Tipify output of PCA
    common_names <- intersect(names(ret_clean$m_tip_pca), colnames(dat));
    ret <- tipify(dat[,common_names,with=F],ret_clean$m_tip_pca[common_names],ret_clean$s_tip_pca[common_names]);
    dat <-  data.table(dat[, fixed_var, with = F], ret$dat);
  }
  
  
  # Return output
  return(dat);
}

