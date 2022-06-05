######################### MAE ##########################################
MAE <- function(real, pred){
  return(mean(abs(real-pred)));
}

######################### MAPE ##########################################
MAPE <- function(real, pred){
  return(100*mean(abs(real-pred)/abs(sapply(real, max, 1))));
}

######################### MSE ##########################################
MSE <- function(real, pred){
  return(mean((real-pred)^2));
}

######################### ROC CURVE #########################
AUC <- function(target, prediction){
  roc_curve <- roc(target, prediction);
  auc_value <- auc(roc_curve);
  return(list(roc_curve = roc_curve, auc_value = auc_value));
}




######################### f1 score ##########################################
f1_score <- function(real, pred, all = FALSE){
  f1 <- c();
  classes <- sort(unique(real))
  if (length(classes) == 2){classes <- classes[-1]}
  for (class in classes){
    TP <- sum(pred == class & real == class);
    FP <- sum(pred == class & real != class);
    FN <- sum(pred != class & real == class);
    TN <- sum(pred != class & real != class);
    precision <- TP / (TP + FP);
    recall <- (TP) / (TP + FN);
    f1 <- c(f1, 2 * (precision*recall) / (precision+recall));
  }
  if (all){
    return(list(f1=f1, mean_f1=mean(f1), precision=precision, recall=recall))
  } else {
    return(f1);
  }
}



######################### ACCURACY ##########################################
acc <- function(real, pred, type = "normal"){
  if (type == "normal"){
    acc_value <- sum(real == pred)/length(real);
  } else if (type == "sentiment") {
    acc_value <- sum(real == pred | pred == 1)/length(real);
  }
  return(acc_value);
}



######################### REPLACE ##########################################
# Replace null values with mean (first row)
replace <- function(x){
  m<-x[1];
  x<-x[-1];
  x[is.na(x)]<-m;
  return(x);
}

######################### UNIFY NULL VALUES ##########################################

unify_null <- function(x){
  na_index<-(is.na(x) | x=="null" | x=="NULL" | x=='\\N' | x=='\\n' | x=="");
  x[na_index]<-"NA";
  return(x)
}

######################### TO NUMERIC ##########################################

to_numeric <- function(x){
  if (class(x)!="numeric"){
    x<-as.numeric(x)
  }
  x <- round(x,2);
  return(x)
}

#################################### SELECT VARIABLES #################################

select_variables <- function(dat,var="ALL",variables_total,categorical_total,no_categorical_total,
                             variables,categorical,no_categorical){
  # Selected variables for the model
  if (var=='ALL'){
    dat <- dat[,variables_total[variables_total %in% colnames(dat)],with=F];
    categorical <- categorical_total;
    no_categorical <- no_categorical_total;
  }else{
    dat <- dat[,variables[variables %in% colnames(dat)],with=F];
  }
  
  return(list(dat=dat,categorical=categorical,no_categorical=no_categorical));
}




######################### REMOVE CONSTANT ##########################################
constant_variables <- function(dat){
  unique_values<-apply(dat,2,function(x){length(unique(x))})
  constant_variables<-names(unique_values)[unique_values == 1];
  return(constant_variables);
}

###############################  REMOVE NOT INFORMED VARIABLES ############################## 

not_informed <- function(dat,threshold_nas=0.4){
  ratio_nas<-apply(dat,2,function(x){sum(is.na(x) | x=="NA")/length(x)});
  not_informed_var<-names(ratio_nas)[ratio_nas>threshold_nas];
  return(not_informed_var);
}

######################### FILL_VALUES ##########################################
# Fill missing values
fill_values <- function(dat,fill,no_categorical){
  if (sum(sapply(dat[,no_categorical,with=F],function(x){sum(is.na(x))})) > 0){
    # Get original order;
    dat$ord<-1:nrow(dat);
    
    # Generalize dataset
    if (sum(colnames(dat)=="dataset")==0){
      dat$dataset <- "dummy";
    }
    
    # Fill missing values with mean
    if (fill=="mean"){
      m<-apply(dat[dataset!="original_train",no_categorical,with=F],2,mean,na.rm=TRUE);
      dat[,no_categorical]<-data.table(apply(rbind(t(m),dat[,no_categorical,with=F]),2,replace));
      
      # Fill missing using mice package
    } else if (fill=="mice"){
      n_m<-5;
      for (iter in 1:2){
        try(imputed<-mice(dat[,-"target",with=F],maxit=5,m=n_m));
        if (exists("imputed")){
          res<-complete(imputed,1)[,no_categorical];
          for (i in 2:n_m){
            res<-res+complete(imputed,i)[,no_categorical];
          }
          dat[,no_categorical]<-res/n_m;
        }
      }
      
      # Fill persisting missing values with mean
      m<-apply(dat[,no_categorical,with=F],2,mean,na.rm=TRUE);
      dat[,no_categorical]<-apply(rbind(t(m),dat[,no_categorical,with=F]),2,replace);
      
      # Fill missing using Amelia package
    } else if (fill=="amelia"){
      nominal<-setdiff(colnames(dat),c("target","index",no_categorical))
      idvars<-c("index","order",colnames(dat)[sapply(dat,class)=="character" & sapply(dat,function(x){sum(is.na(x))}) > 0]);
      noms<-setdiff(colnames(dat)[sapply(dat,class)=="character"],idvars)
      
      dat[,no_categorical]<-amelia(dat[,-"target",with=F],noms=noms,
                                   idvars=idvars,m = 1,parallel="multicore",cl=4)$imputations[[1]][,no_categorical];
      
      # Fill persisting missing values with mean
      m<-apply(dat[,no_categorical,with=F],2,mean,na.rm=TRUE);
      dat[,no_categorical]<-apply(rbind(t(m),dat[,no_categorical,with=F]),2,replace);
    }
    dat<-dat[order(ord)];
    dat<-dat[,-"ord",with=F];
  }
  
  # Generalize dataset
  if (sum(dat$dataset=="dummy")==nrow(dat)){
    dat <- dat[,-"dataset",with=F]
  }
  
  return(dat);
}

######################### UNTOUCHABLE VALUES ##########################################
# Compute "untouchable" categories (the ones with high frequency on fraud transactions)
untouchable_values_f <- function (name,x,index,index_blacklist,index_reject,index_session_id,index_suspected,metric_weights,
                                  min_transactions=2, anomaly_factor=1){
  # Format values
  values<-as.character(t(x));
  
  # List of categories
  unique_values<-unique(values);
  
  # Computed weighted target
  weighted_target<-(as.numeric(index %in% index_blacklist)*metric_weights[1]+as.numeric(index %in% index_reject)*metric_weights[2]
                    +as.numeric(index %in% index_session_id)*metric_weights[3]+as.numeric(index %in% index_suspected)*metric_weights[4])/sum(metric_weights);
  
  # Get total fraud weight per category
  fraud_frequency<-sapply(unique_values,function(x,values,weighted_target){mean(weighted_target[values==x])},
                          values,weighted_target);
  total_frequency<-sapply(unique_values,function(x,values,weighted_target){sum(values==x)},
                          values);
  m<-mean(fraud_frequency);
  s<-sd(fraud_frequency);
  
  # Set untouchable variables as the one with anomaly fraud frequency
  untouchable<-names(fraud_frequency)[((fraud_frequency > (m + anomaly_factor*s)) | (fraud_frequency < (m - anomaly_factor*s))) & total_frequency >= min_transactions];
  
  if (length(untouchable)>0){
    return(data.table(name=name,untouchable=list(untouchable)))
  } else {
    return(data.table(name=name,untouchable=NA));
  }
}

######################### GROUP CATEGORIES LIST ##########################################
# Get list of categories to group
group_categories_list_modified <- function (name,x,threshold_variables,min_size_others=2,max_n_categories=10,
                                            untouchable_values=data.table()){
  
  # Reformat values
  x<-as.character(t(x));
  
  grouped_values<-NA;
  ind_others<- (table(x)/length(x)) <  threshold_variables;
  if (length(ind_others) >= max_n_categories & sum(ind_others) > min_size_others){
    grouped_values<-names(ind_others)[ind_others];
    if(length(untouchable_values)>0){
      grouped_values<-grouped_values[!(grouped_values %in% unlist(untouchable_values$untouchable[untouchable_values$name==name]))]
    }
  }
  
  if (!is.na(grouped_values)[1] & length(grouped_values)>0){
    return(data.table(name=name,grouped_values=grouped_values)) # used to be grouped_values=list(grouped_values)
  } else {
    return(data.table(name=name,grouped_values=NA));
  }
}



######################### GROUP CATEGORIES ##########################################
# Group values with few instances
group_values <- function (name,values_to_group,x){
  
  # Reformat values
  x<-as.character(t(x));
  
  # Group categories
  x[x%in% values_to_group]<-"others"; 
  
  return(x);
}

######################### NEW TO OTHERS ##########################################
new_to_others<- function(x,x_train){
  # Reformat values
  x<-as.character(t(x));
  x_train<-as.character(t(x_train));
  
  train_values<-unique(x_train);
  x[!(x %in% train_values)]<-"others";
  return(x);
}

######################### GROUP CATEGORIES ##########################################
# Group categories
group_categories <- function(categorical,dat,dat_train,threshold_variables,min_size_others,max_n_categories,untouchable_values){
  
  # Get list of categories to group
  grouped_values_list<-data.table(t(sapply(categorical, 
                                           function(x){group_categories_list(x,dat_train[,x,with=F],threshold_variables,
                                                                             min_size_others,max_n_categories,untouchable_values)})));
  # Apply categories grouping                                                                                                                                                      min_size_others,max_n_categories,untouchable_values)})));
  dat[,categorical]<-data.table(sapply(categorical, 
                                       function(x){group_values(x,unlist(grouped_values_list[name==x]$grouped_values),dat[,x,with=F])}));
  
  # Set values of val and test not in train to others
  dat[,categorical]<-data.table(sapply(categorical,function(x){new_to_others(dat[,x,with=F],dat_train[,x,with=F])}));
  
  
  
  return(list(dat=dat,grouped_values_list=grouped_values_list));
  
}



######################### ONE HOT ENCODING ##########################################
onehot <- function(x){
  values<-unique(x);
  
  ret <- matrix(0,nrow = length(x),ncol=length(values));
  
  for (i in 1:length(values)){
    ret[, i] = as.numeric(x==values[i]);
  }
  
  colnames(ret)<-values;
  return(ret);
}

one_hot_encoding <- function(dat,categorical,no_categorical,fixed_var = NA){
  
  cat_dat <- as.matrix(dat[, categorical, with = F]);
  vars <- sample(categorical, length(categorical));
  ret <- foreach (c = vars, .combine = cbind) %dopar%{
    inner_ret <- onehot(cat_dat[,c]);
    colnames(inner_ret) <- paste0(c,colnames(inner_ret));
    inner_ret;
  }
  ret <- data.table(ret);
  
  if (length(no_categorical) > 0){
    if (!is.na(fixed_var)[1]){
      dat <- data.table(dat[,fixed_var,with=F],ret, dat[,no_categorical,with=F]);
    } else {
      dat <- data.table(ret, dat[,no_categorical,with=F]);
    }
  } else {
    if (!is.na(fixed_var)[1]){
      dat <- data.table(dat[,fixed_var,with=F],ret);
    } else {
      dat <- data.table(ret);
    }
  }
  
  return(dat);
}


#################################### TIPIFY #################################
tipify <- function(dat,m,s, remove_constant = TRUE){
  
  # Convert to data.table
  if (class(dat)[1]!="data.table"){
    dat<-data.table(dat);
  }
  
  # Remove constant variables
  constant_var<-which(s==0);
  if (length(constant_var)>0){
    if (remove_constant){
      dat<-dat[,-constant_var,with=F];
      m<-m[-constant_var];
      s<-s[-constant_var];
      print(sprintf("%s constant variables removed",length(constant_var)));
    }
    else {
      print(sprintf("%s constant variables NOT removed",length(constant_var)));
      s[constant_var] <- 1;
    }
  }
  
  dat <- as.matrix(dat);
  for (i in 1:ncol(dat)){
    dat[,i] <-  (dat[,i]-m[i])/s[i];
  }
  dat <- data.table(dat)
  
  colnames(dat)<-names(m);
  
  return(list(dat=dat,m=m,s=s));
}

#################################### UNTIPIFY #################################
untipify <- function(dat,m,s){
  
  # Convert to data.table
  if (class(dat)[1]!="data.table"){
    dat<-data.table(dat);
  }
  
  # Remove constant variables
  constant_var<-which(s==0);
  if (length(constant_var)>0){
    dat<-dat[,-constant_var,with=F];
    m<-m[-constant_var];
    s<-s[-constant_var];
  }
  print(sprintf("%s constant variables removed",length(constant_var)));
  
  dat <- as.matrix(dat);
  for (i in 1:ncol(dat)){
    dat[,i] <-  dat[,i]*s[i] + m[i];
  }
  dat <- data.table(dat)
  
  
  colnames(dat)<-names(m);
  
  return(list(dat=dat,m=m,s=s));
}

#################################### IRRELEVANT VARIABLES #################################

remove_irrelevant<-function(correlations,irrelevant_threshold){
  # Irrelevant variables
  relevance<-correlations[1,-1];
  irrelevant_variables<-names(relevance)[is.na(relevance) | relevance<irrelevant_threshold];
  return(irrelevant_variables);
}

#################################### REDUNDANT VARIABLES #################################
remove_redundant <- function(correlations,redundant_threshold){
  redundancy<-apply(correlations,2,function(x){which(x>redundant_threshold)});
  redundancy<-redundancy[which(sapply(redundancy,length)>1)]
  
  redundant_variables<-c();
  for (i in redundancy){
    imp<-sort(correlations[1,i],decreasing = TRUE);
    redundant_variables<-c(redundant_variables,names(imp)[2:length(i)])
  }
  redundant_variables<-unique(redundant_variables);
  return(redundant_variables)
} 

#################################### IMPORTANT VARIABLES #################################
select_important<-function(dat,n_perc_varimp,priority_class=2,y=NA){
  if (is.na(y[1])){
    varimp<-filterVarImp(x=dat[,-fixed_var,with=F],y=dat$target,nonpara=TRUE);
  } else {
    varimp<-filterVarImp(x=dat,y=y,nonpara=TRUE);
  }
  varimp<-data.table(variable=rownames(varimp),imp=varimp[,priority_class])
  varimp<-varimp[order(-imp)];
  selected<-varimp$variable[1:ceiling(n_perc_varimp*nrow(varimp))];
  return(selected);
}

#################################### PCA #################################
pca <- function(dat,dat_train,threshold_variance, fixed_var){
  # PCA
  prin_comp <- prcomp(dat_train[,-fixed_var,with=F],center = FALSE,scale. = FALSE); 
  
  # choose number of variables
  prop_varex <- summary(prin_comp)$importance[3,];
  number_of_variables <- min(which(prop_varex > threshold_variance))
  
  # Apply PCA
  dat<-data.table(dat[,fixed_var,with=F],
                  data.table(predict(prin_comp, newdata = dat[,-fixed_var,with=F]))[,1:number_of_variables,with=F]);
  
  # Get mean and standard deviation for train data (excluding target) 
  m_tip<-apply(dat[dataset=="train",-fixed_var,with=F],2,mean);
  s_tip<-apply(dat[dataset=="train",-fixed_var,with=F],2,sd);
  
  # Tipify data
  ret<-tipify(dat[,-fixed_var,with=F],m_tip,s_tip);
  dat<-data.table(dat[,fixed_var,with=F],ret$dat);
  m_tip<-ret$m;
  s_tip<-ret$s;
  
  return(list(dat=dat,prin_comp=prin_comp,m_tip=m_tip,s_tip=s_tip));
}






