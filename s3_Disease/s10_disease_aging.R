rm(list = ls())
# Set your working directory and data path here

library(survival)
packageVersion("survival")
CoxReg <- function(status, x, co, data) {

  cat(paste("Running Cox regression for", status, x, "\n"))
  time <- paste0(status, "_duration")
  data <- data[!is.na(data[[x]]) & !is.na(data[[status]]) & !is.na(data[[time]]), ]
  formula <- paste("Surv(", time, ",", status, ") ~", paste(c(x, co), collapse = '+'))

  cox_model <- coxph(as.formula(formula), data = data)

  cox_summary <- summary(cox_model)
  cox_ci <- exp(confint(cox_model))

  hr <- paste0(exp(cox_summary$coefficients[1, 1]))
  lower <- paste0(cox_ci[1, 1])
  upper <- paste0(cox_ci[1, 2])


  cox.Coeff_B <- paste0(
    round(cox_summary$coefficients[1, 2], 2),
    ' (', round(cox_ci[1, 1], 2), ', ', round(cox_ci[1, 2], 2), ')'
  )

  cox.p_B <- cox_summary$coefficients[1, 5]

  cox.coeff <- data.frame(
    'Characteristics' = x,
    'y' = status,
    'HR' = c(hr),
    'Lower' = c(lower),
    'Upper' = c(upper),
    'Coef[95% CI]' = c(cox.Coeff_B),
    'P value' = c(cox.p_B)
  )

  return(cox.coeff)
}

res_cox <- function(status, Var.analys, co, data) {
  res <- data.frame()

  for (y in status) {
    for (i in 1) {
      combos <- combn(Var.analys, i)
      for (j in 1:ncol(combos)) {
        predictors <- combos[, j]
        res <- rbind(res, CoxReg(status, predictors, co, data))
      }
    }
  }

  return(res)
}

CoxReg_4class <- function(time, status, x, co, data) {

  data[[x]] <- as.factor(data[[x]])

  formula <- paste("Surv(", time, ",", status, ") ~", paste(c(x, co), collapse = '+'))

  cox_model <- coxph(as.formula(formula), data = data)

  cox_summary <- summary(cox_model)
  cox_ci <- exp(confint(cox_model))

  coeffs <- c()
  p_vals <- c()
  HRs <- c()
  Low_CIs <- c()
  Up_CIs <- c()

  for (i in 2:4) {
    cox.Coeff_B <- paste0(
      round(cox_summary$coefficients[i-1, 2], 2),
      ' (', round(cox_ci[i-1, 1], 2), ', ', round(cox_ci[i-1, 2], 2), ')'
    )
    cox.p_B <- cox_summary$coefficients[i-1, 5]
    HRs <- c(HRs, exp(cox_summary$coefficients[i-1, 1]))
    Low_CIs <- c(Low_CIs, cox_ci[i-1, 1])
    Up_CIs <- c(Up_CIs, cox_ci[i-1, 2])
    coeffs <- c(coeffs, cox.Coeff_B)
    p_vals <- c(p_vals, cox.p_B)
  }
  cases <- list()
  totals <- list()
  for (i in 1:4) {
    cases[[i]] <- sum(data[[status]][data[[x]] == levels(data[[x]])[i]], na.rm = TRUE)
    totals[[i]] <- sum(!is.na(data[[status]][data[[x]] == levels(data[[x]])[i]]))
  }

  cox.coeff <- data.frame(
    'Characteristics' = x,
    'Predictors' = paste(x, levels(data[[x]])[2:4], sep = '_'),
    'y' = status,
    'Cases' = unlist(cases[2:4]),
    'Totals' = unlist(totals[2:4]),
    'HR' = HRs,
    'Lower' = Low_CIs,
    'Upper' = Up_CIs,
    'Coef[95% CI]' = coeffs,
    'P value' = p_vals
  )

  return(cox.coeff)
}

res_4class_cox <- function(time, status, Var.analys, co, data) {
  res <- data.frame()

  for (y in status) {
    for (i in 1) {
      combos <- combn(Var.analys, i)
      for (j in 1:ncol(combos)) {
        predictors <- combos[, j]

        res <- rbind(res, CoxReg_4class(time, status, predictors, co, data))
      }
    }
  }

  return(res)
}
library(data.table)
Myt <- fread("data/Food/MYT_score_either.csv")

Myt <- Myt[, .(eid, MYT)]
Myt$MYT_scale <- scale(Myt$MYT)

df <- fread("data/Y_total.csv")
df <- merge(df, Myt, by = 'eid', all = F)

disease_df <- data.frame()
disease_dir <- 'data/Disease'
disease_lst_df <- readxl::read_xlsx('data/Disease/disease_code_Aging.xlsx')
colnames(disease_lst_df)
disease_lst <- disease_lst_df$Disease
disease <- disease_lst[1]
low_case <- c()
nb_case <- c()
for (disease in disease_lst) {
  print(disease)
  subfolder <- disease_lst_df$Category[disease_lst_df$Disease == disease]
  subfolder <- gsub(' ', '_', subfolder)
  subfolder <- paste0('Aging_', subfolder)
  disease <- gsub(' ', '_', disease)
  disease <- gsub(',', '', disease)
  file_name <- paste0(disease_dir, '/', subfolder, '/', disease, '.csv')
  temp <- fread(file_name)
  duration_col <- paste0(disease,'_duration')
  temp <- temp[get(duration_col) > 0, .(eid, get(disease), get(duration_col))]
  colnames(temp) <- c('eid', disease, duration_col)
  case_nb <- nrow(temp[temp[[disease]] == 1, ])
  if (case_nb < 500) {
    cat(paste(disease, 'has fewer than 500 cases, skipping\n'))
    low_case <- c(low_case, disease)
    nb_case <- c(nb_case, case_nb)
    next
  }
  if (nrow(disease_df) == 0) {
    disease_df <- temp
  } else {
    disease_df <- merge(disease_df, temp, by = 'eid', all = T)
  }
}
low_dis <- data.frame('Disease' = low_case, 'Case_nb' = nb_case)
fwrite(low_dis, 'results/s7_disease/low_case_disease_Aging.csv')

cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ','energy','Vitamin','fast_time')
df <- merge(df, disease_df, by = 'eid', all.x = T)
status <- 'Oesophagus'


disease_lst <- gsub(' ', '_', disease_lst)
disease_lst <- gsub(',', '', disease_lst)
disease_lst <- setdiff(disease_lst, low_case)

colnames(df) <- enc2utf8(colnames(df))
colnames(df) <- gsub(" ", " ", colnames(df), useBytes = TRUE)
colnames(df) <- trimws(colnames(df))
colnames(df) <- gsub("\\s*_\\s*", "_", colnames(df), perl = TRUE)
colnames(df) <- gsub("-", "_", colnames(df), perl = TRUE)
disease_lst <- enc2utf8(disease_lst)
disease_lst <- gsub(" ", " ", disease_lst, useBytes = TRUE)
disease_lst <- trimws(disease_lst)
disease_lst <- gsub("\\s*_\\s*", "_", disease_lst, perl = TRUE)
disease_lst <- gsub("-", "_", disease_lst, perl = TRUE)


disease_lst <- setdiff(disease_lst, 'Hyperlipidemia')
results <- data.frame()
for (disease in disease_lst) {
  temp <- res_cox(
    status = disease,
    Var.analys = c('MYT_scale'),
    co = cov,
    data = df
  )

  if (nrow(results) == 0) {
    results <- temp
  } else {
    results <- rbind(results, temp)
  }
}

results$P_fdr<- p.adjust(results$P.value, method = 'fdr')
fwrite(results,'Results/s7_disease/Disease_Aging_MYT_continuous.csv', row.names = F)

sum(results$P_fdr < 0.05)

Myt <- df
Myt$MYT_4cat <- cut(Myt$MYT,
                    breaks = quantile(Myt$MYT, probs = seq(0, 1, by = 0.25), na.rm = TRUE),
                    labels = c("0-2", "3", "4", "5+"),include.lowest = TRUE)
Myt$MYT_4cat <- factor(Myt$MYT_4cat, labels = c("0-2", "3", "4", "5+"))
df <- Myt

results_4class <- data.frame()
for (disease in disease_lst) {
  cat(paste0("Processing ", disease, "\n"))
  temp <- res_4class_cox(
    time = paste0(disease, '_duration'),
    status = disease,
    Var.analys = c('MYT_4cat'),
    co = cov,
    data = df
  )

  if (nrow(results_4class) == 0) {
    results_4class <- temp
  } else {
    results_4class <- rbind(results_4class, temp)
  }
}
results_4class$P_fdr<- p.adjust(results_4class$P.value, method = 'fdr')
fwrite(results_4class,'Results/s7_disease/Disease_Aging_MYT_4class.csv', row.names = F)
