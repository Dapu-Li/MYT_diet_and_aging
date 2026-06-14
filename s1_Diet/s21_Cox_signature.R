rm(list = ls())
# Set your working directory and data path here

library(survival)

CoxReg <- function(time, status, x, co, data) {

  cat(paste("Running Cox regression for", status, x, "\n"))

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

res_cox <- function(time, status, Var.analys, co, data) {
  res <- data.frame()

  for (y in status) {
    for (i in 1) {
      combos <- combn(Var.analys, i)
      for (j in 1:ncol(combos)) {
        predictors <- combos[, j]
        res <- rbind(res, CoxReg(time, status, predictors, co, data))
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
df <- fread('data/Y_total.csv')
MYT_sig <- fread('results/s5_Lasso/MYT_Pro_signature.csv')
df <- merge(df, MYT_sig, by = 'eid', all.y = TRUE)
df$MYT_Pro_signature <- scale(df$MYT_Pro_signature, center = TRUE, scale = TRUE)
df$MYT_Pro_signature_Group <- cut(df$MYT_Pro_signature, breaks = c(-Inf, quantile(df$MYT_Pro_signature, probs = c(0.25, 0.5, 0.75)), Inf),
                                        labels = c('Q1', 'Q2', 'Q3', 'Q4'))
colnames(df)
table(df$death)
cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ','energy','Vitamin','fast_time')
x <- c('MYT_Pro_signature','MYT_Pro_signature_Group')
y <- 'death';duration <- 'death_duration'

df$death <- as.numeric(df$death)
df <- na.omit(df)
cox_num <- res_cox(duration, y, 'MYT_Pro_signature', cov, df)
cox_4class <- res_4class_cox(duration, y, 'MYT_Pro_signature_Group', cov, df)


cox_num$HR    <- as.numeric(cox_num$HR)
cox_num$Lower <- as.numeric(cox_num$Lower)
cox_num$Upper <- as.numeric(cox_num$Upper)
cox_4class$HR    <- as.numeric(cox_4class$HR)
cox_4class$Lower <- as.numeric(cox_4class$Lower)
cox_4class$Upper <- as.numeric(cox_4class$Upper)
cox_num$HR_Lower <- cox_num$HR - cox_num$Lower
cox_num$HR_Upper <- cox_num$Upper - cox_num$HR
cox_4class$HR_Lower <- cox_4class$HR - cox_4class$Lower
cox_4class$HR_Upper <- cox_4class$Upper - cox_4class$HR

fwrite(cox_num, 'Results/s5_Lasso/cox_num_signature.csv', row.names = F)

fwrite(cox_4class, 'Results/s5_Lasso/cox_4class_signature.csv', row.names = F)
