rm(list = ls())
# Set your working directory and data path here

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

Y_total <- fread('data/Y_total.csv')
CDAI_DASH <- fread('data/Food/CDAI_DASH_scores.csv')
HEI <- fread('data/Food/HEI_2010.csv')
MYT <- fread('data/Food/MYT_score_either.csv')
cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ','energy','Vitamin','fast_time')

df <- merge(Y_total, HEI, by='eid', all.x=TRUE)
df <- merge(df, MYT, by='eid', all.x=TRUE)
df <- merge(df, CDAI_DASH, by='eid', all.x=TRUE)
colnames(df)
df$MYT <- scale(df$MYT)
df$HEI_2010 <- scale(df$HEI_2010)
df$CDAI <- scale(df$CDAI_score)
df$DASH <- scale(df$DASH_score)
colnames(df)
library(survival)
diet <- c('HEI_2010','MYT','CDAI','DASH')
res_val <- res_cox(time = 'death_duration', status = c('death'),
                Var.analys = c('HEI_2010','MYT','CDAI','DASH'), co = cov, data = df)

for (d in diet) {
  cat(paste0("Processing ", d, "\n"))
  quantiles <- quantile(df[[d]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  df[[paste0(d, '_Category')]] <- cut(df[[d]],
                                      breaks = unique(quantiles),
                                      include.lowest = TRUE,
                                      labels = c('Q1', 'Q2', 'Q3', 'Q4'))
}
res_4class_val <- res_4class_cox(time = 'death_duration', status = c('death'),
                        Var.analys = c('HEI_2010_Category','MYT_Category','CDAI_Category','DASH_Category'), co = cov, data = df)
fwrite(res_val, 'results/s11_val/Cox_validation_diet_scores_continuous.csv', row.names = FALSE)
fwrite(res_4class_val, 'results/s11_val/Cox_validation_diet_scores_4class.csv', row.names = FALSE)
