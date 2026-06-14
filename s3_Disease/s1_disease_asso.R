rm(list = ls())
# Set your working directory and data path here

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

CoxReg_4class <- function(status, x, co, data) {
  print(paste("Processing variable:",status, x))
  data[[x]] <- as.factor(data[[x]])

  time <- paste0(status, "_duration")
  formula <- paste("Surv(", time, ",", status, ") ~", paste(c(x, co), collapse = '+'))

  cox_model <- coxph(as.formula(formula), data = data)

  cox_summary <- summary(cox_model)
  cox_ci <- exp(confint(cox_model))

  coeffs <- c()
  p_vals <- c()

  for (i in 2:5) {
    cox.Coeff_B <- paste0(
      round(cox_summary$coefficients[i-1, 2], 2),
      ' (', round(cox_ci[i-1, 1], 2), ', ', round(cox_ci[i-1, 2], 2), ')'
    )
    cox.p_B <- round(cox_summary$coefficients[i-1, 5], 3)

    coeffs <- c(coeffs, cox.Coeff_B)
    p_vals <- c(p_vals, cox.p_B)
  }

  cox.coeff <- data.frame(
    'Characteristics' = x,
    'Predictors' = paste(x, levels(data[[x]])[2:5], sep = '_'),
    'HR' = cox_summary$coefficients[1:4, 2],
    'Lower CI' = cox_ci[1:4, 1],
    'Upper CI' = cox_ci[1:4, 2],
    'Coef[95% CI]' = coeffs,
    'P value' = p_vals
  )

  return(cox.coeff)
}

res_4class_cox <- function(status, Var.analys, co, data) {
  res <- data.frame()

  for (y in status) {
    for (i in 1) {
      combos <- combn(Var.analys, i)
      for (j in 1:ncol(combos)) {
        predictors <- combos[, j]

        res <- rbind(res, CoxReg_4class(status, predictors, co, data))
      }
    }
  }

  return(res)
}


library(data.table)
Myt <- fread("data/Food/MYT_score_noBeer.csv")
Myt <- Myt[, .(eid, MYT)]
Myt$MYT_5cat <- cut(Myt$MYT,
                    breaks = quantile(Myt$MYT, probs = seq(0, 1, by = 0.20), na.rm = TRUE),
                    labels = c("Q1", "Q2", "Q3", "Q4",'Q5'),include.lowest = TRUE)
table(Myt$MYT_5cat)
Q1 <- Myt[MYT_5cat == "Q1",]
Q2 <- Myt[MYT_5cat == "Q2",]
Q3 <- Myt[MYT_5cat == "Q3",]
Q4 <- Myt[MYT_5cat == "Q4",]
Q5 <- Myt[MYT_5cat == "Q5",]
table(Q1$MYT)
table(Q2$MYT)
table(Q3$MYT)
table(Q4$MYT)
table(Q5$MYT)
Myt$MYT_5cat <- factor(Myt$MYT_5cat, labels = c("0-2", "3", "4", "5",'6+'))
table(Myt$MYT_5cat)

df <- fread("data/Y_total.csv")
df <- merge(df, Myt, by = 'eid', all = F)

disease_df <- data.frame()
disease_dir <- 'data/Disease'
disease_lst_df <- readxl::read_xlsx('data/Disease/disease_code.xlsx')
colnames(disease_lst_df)
disease_lst_df <- subset(disease_lst_df, Category != 'Mortality')
disease_lst <- disease_lst_df$Disease
disease <- disease_lst[1]
for (disease in disease_lst) {
  subfolder <- disease_lst_df$Category[disease_lst_df$Disease == disease]
  subfolder <- gsub(' ', '_', subfolder)
  disease <- gsub(' ', '_', disease)
  file_name <- paste0(disease_dir, '/', subfolder, '/', disease, '.csv')
  temp <- fread(file_name)
  duration_col <- paste0(disease,'_duration')
  temp <- temp[get(duration_col) > 0, .(eid, get(disease), get(duration_col))]
  colnames(temp) <- c('eid', disease, duration_col)
  if (nrow(disease_df) == 0) {
    disease_df <- temp
  } else {
    disease_df <- merge(disease_df, temp, by = 'eid', all = T)
  }
}

cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ','energy','Vitamin','fast_time')
df <- merge(df, disease_df, by = 'eid', all.x = T)
library(survival)
disease_lst <- gsub(' ', '_', disease_lst)
disease_lst1 <- disease_lst[1]
results <- res_cox(
  status = disease_lst[2],
  Var.analys = c('MYT'),
  co = cov,
  data = df
)

for (disease in disease_lst) {
  temp <- res_cox(
    status = disease,
    Var.analys = c('MYT'),
    co = cov,
    data = df
  )
  results <- rbind(results, temp)
}
