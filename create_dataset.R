library(data.table);


input_path <- "/home/ygala/TFM_UOC/"


#################################################### READ DATA ################################################################
# IMAGE METADATA
dat_metadata <- fread("/home/jesusprada/proyecto_python/x-ray/data/covid_images/clean_data.csv") # 1814 IDs
colnames(dat_metadata) <- tolower(colnames(dat_metadata));
vars_metadata <- c("patientid",           
                   "filename", "studydate", "studytime", "modality", "viewposition", "seriesdescription",
                   "patientorientation", "photometricinterpretation",
                   "rows", "columns", "pixelspacing");
dat_metadata <- dat_metadata[, vars_metadata, with = F];
patient_ids <- unique(dat_metadata$patientid)

# CMBD
dat_1 <- fread(file.path(input_path, "01.csv"))
colnames(dat_1) <- tolower(colnames(dat_1));
vars_1 <- c("patient id", "edad/age", "sexo/sex", "diag ing/inpat", "f_ingreso/admission_d_ing/inpat", "f_entrada_uc/icu_date_in",
            "f_ingreso/admission_date_urg/emerg", "hora/time_admision/admission_urg/emerg", "especialidad/department_urg/emerg",
            "diag_urg/emerg");
dat_1 <- dat_1[dat_1$`patient id` %in% patient_ids, vars_1, with = F];

# FARMACIA (SACAR LOS MAS COMUNES EN COVID. EJ: PAR = 0/1, valor_PAR = 0/0.7)
dat_2 <- fread(file.path(input_path, "02.csv"))
colnames(dat_2) <- tolower(colnames(dat_2));
vars_2 <- c("patient id", "farmaco/drug_nombre_comercial/comercial_name", "dosis_media_diaria/daily_avrg_dose",
            "inicio_trat/drug_start_date", "fin_trat/drug_end_date", "id_atc5", "id_atc7")
dat_2 <- dat_2[dat_2$`patient id` %in% patient_ids, vars_2, with = F];
  

# CONSTANTES (SACAR PRIMERO, ULTIMO, MAXIMO, MINIMO, MEDIA)
dat_3 <- fread(file.path(input_path, "03.csv"))
colnames(dat_3) <- tolower(colnames(dat_3));
vars_3 <- c("patient id",           
            "fc/hr_ing/inpat", "ta_max_ing/inpat", "ta_min_ing/inpat",  "temp_ing/inpat",
            "constants_ing/inpat_fecha/date", "constants_ing/inpat_hora/time");
dat_3 <- dat_3[dat_3$`patient id` %in% patient_ids, vars_3, with = F];

# LABORATORIO (SACAR LOS MAS COMUNES EN COVID. EJ: CO2 (determinacion) = 0/1 resultado_CO2 = 0/0.7, ref_CO2_min, ref_CO2_max)
# dat_4 <- fread(file.path(input_path, "04.csv"))
# colnames(dat_4) <- tolower(colnames(dat_4));
# vars_4 <- c("patient id",           
#             "determinacion/item_lab", "resultado/val_result", "valores_referencia/ref_values");

# DIAGNOSTICOS Y PROCEDIMIENTOS (URG)
dat_5 <- fread(file.path(input_path, "05.csv"))
colnames(dat_5) <- tolower(colnames(dat_5));
vars_5 <- c("patient id",           
            "dia_ppal",   "dia_02",
            "proc_01",    "proc_02");
dat_5 <- dat_5[dat_5$`patient id` %in% patient_ids, vars_5, with = F];

# DIAGNOSTICOS Y PROCEDIMIENTOS (HOSP)
dat_6 <- fread(file.path(input_path, "06.csv"))
colnames(dat_6) <- tolower(colnames(dat_6));
vars_6 <- c("patient id",           
            "dia_ppal",   "dia_02",     "dia_03",
            "poad_ppal",   "poad_02",     "poad_03", 
            "proc_01",    "proc_02", "proc_03");
dat_6 <- dat_6[dat_6$`patient id` %in% patient_ids, vars_6, with = F];


#################################################### COMBINE DATA ################################################################

### [1] Merge with CMBD
dat <- merge(dat_metadata, dat_1,  by.x = "patientid", by.y = "patient id", all.x = TRUE);

### [2] Merge with farmacia

# Only data in the first 4 days
dat_2 <- merge(dat_2, dat_1[, c('patient id', 'f_ingreso/admission_d_ing/inpat'), with = F], by = "patient id");
dat_2$fec_ing <-  as.Date(dat_2$`f_ingreso/admission_d_ing/inpat`,format = '%d/%m/%Y');
dat_2$fec_toma <-  as.Date(dat_2$`inicio_trat/drug_start_date`,format = '%d/%m/%Y');
dat_2$days <- as.numeric(dat_2$fec_toma - dat_2$fec_ing);
dat_2 <- dat_2[days <= 4];

# Get only last dose for each patient
dat_2 <- dat_2[order(fec_toma)];
dat_2 <- dat_2[, .SD[1], by = c('patient id', 'farmaco/drug_nombre_comercial/comercial_name')]

# Most frequent drugs for exitus/non-exitus
dat_2 <- merge(dat_2,  fread(file.path(input_path, "01.csv"))[, c('PATIENT ID', 'MOTIVO_ALTA/DESTINY_DISCHARGE_ING'), with = F], 
               by.x = "patient id", by.y = 'PATIENT ID');
drugs_0 <- head(names(sort(table(dat_2[`MOTIVO_ALTA/DESTINY_DISCHARGE_ING` == "Domicilio"]$`farmaco/drug_nombre_comercial/comercial_name`), decreasing = T)), 10)
drugs_1 <- head(names(sort(table(dat_2[`MOTIVO_ALTA/DESTINY_DISCHARGE_ING` == "Fallecimiento"]$`farmaco/drug_nombre_comercial/comercial_name`), decreasing = T)), 10)
drugs <- unique(c(drugs_0, drugs_1));
dat_2 <- dat_2[`farmaco/drug_nombre_comercial/comercial_name` %in% drugs];

# Get one row per patient
get_drugs <- function(dat, drugs){
  res <- data.table(`patient id` = dat$`patient id`[1]);
  for (d in drugs){
    if (d %in% dat$`farmaco/drug_nombre_comercial/comercial_name`){
      inner_res <- data.table(1, dat$`dosis_media_diaria/daily_avrg_dose`[d == dat$`farmaco/drug_nombre_comercial/comercial_name`])
    } else {
      inner_res <- data.table(0, 0L)
    }
    colnames(inner_res) <- c(d, paste0(d, "_dose"))
    res <- cbind(res, inner_res)
  }
  return(res);
}

dat_2 <- dat_2[,get_drugs(.SD, drugs) , by = "patient id"];

# Merge
dat <- merge(dat, dat_2,  by.x = "patientid", by.y = "patient id", all.x = TRUE);

### [3] Merge with vital signs

# Only data in the first 4 days
dat_3 <- merge(dat_3, dat_1[, c('patient id', 'f_ingreso/admission_d_ing/inpat'), with = F], by = "patient id");
dat_3$fec_ing <-  as.Date(dat_3$`f_ingreso/admission_d_ing/inpat`,format = '%d/%m/%Y');
dat_3$fec_toma <-  as.Date(dat_3$`constants_ing/inpat_fecha/date`,format = '%d/%m/%Y');
dat_3$days <- as.numeric(dat_3$fec_toma - dat_3$fec_ing);
dat_3 <- dat_3[days <= 4];

# Get one row per patient
get_constants <- function(dat){
  res <- data.table(`patient id` = dat$`patient id`[1]);
  for (s in setdiff(colnames(dat), "patient id")){
    values <- dat[[s]];
    inner_res <- data.table(values[1], values[length(values)], min(values), max(values), mean(values));
    colnames(inner_res) <- paste0(c("first_", "last_", "min_", "max_", "mean_"),s);
    res <- cbind(res, inner_res)
  }
  return(res);
}

dat_3 <- dat_3[order(fec_toma, `constants_ing/inpat_hora/time`)];
remove_vars <- c("constants_ing/inpat_fecha/date", "constants_ing/inpat_hora/time", "f_ingreso/admission_d_ing/inpat",
                 "fec_ing", "fec_toma", "days");
dat_3 <- dat_3[, -remove_vars, with = F];
dat_3$`temp_ing/inpat`<- as.numeric(gsub(pattern = ",", replacement = ".", dat_3$`temp_ing/inpat`));
dat_3 <- dat_3[,get_constants(.SD) , by = "patient id"];

# Merge
dat <- merge(dat, dat_3,  by.x = "patientid", by.y = "patient id", all.x = TRUE);


### [5] Merge with diag and proc from urgency
colnames(dat_5)[2:ncol(dat_5)] <- paste0("urg_", colnames(dat_5)[2:ncol(dat_5)])
dat <- merge(dat, dat_5,  by.x = "patientid", by.y = "patient id", all.x = TRUE);

### [6] Merge with diag and proc from hospitalization
dat <- merge(dat, dat_6,  by.x = "patientid", by.y = "patient id", all.x = TRUE);


#################################################### SAVE DATA ################################################################
dat <- unique(dat);
saveRDS(dat, file.path(input_path, "tabular_data.RDS"))
