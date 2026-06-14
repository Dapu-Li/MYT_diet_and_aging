rm(list = ls())
# Set your working directory and data path here

CoxReg <- function(time, status, x, co, data) {

  cat(paste("Running Cox regression for", status, x, "\n"))

  vars <- c(x, co)
  vars <- paste0("`", vars, "`")

  formula <- as.formula(
    paste("Surv(", time, ",", status, ") ~", paste(vars, collapse = "+"))
  )
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
library(data.table)
food_code <- fread(paste0(data_dir, "/Food/raw/Food_code.csv"))
Food_Group_lst <- unique(food_code$Food_Group)
Food_Group_lst_cat <- paste0(Food_Group_lst,'_Category')

fg <- fread("data/Food/Food_group_category.csv")[, c('eid', Food_Group_lst_cat), with = F]
colnames(fg)

Food <- Food_Group_lst_cat[1]
for (Food in Food_Group_lst_cat) {
  print(Food)
  print(table(fg[[Food]], useNA = 'ifany'))
  labels <- names(table(fg[[Food]], useNA = "ifany"))
  start_values <- as.numeric(sub("([0-9.]+).*", "\\1", labels))
  fg[[Food]] <- factor(fg[[Food]], labels = start_values)
  fg[[Food]] <- as.numeric(as.character(fg[[Food]]))
  print(table(fg[[Food]], useNA = 'ifany'))
}

df <- fread("data/Y_total.csv")
colnames(df)
df <- merge(df, fg, by = 'eid', all = F)
component <- Food_Group_lst_cat[2]
ref = 2
y <- c('death','death_duration')
library(survival)
library(dplyr)
cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ','energy','Vitamin','fast_time')
res <- res_cox('death_duration','death', Food_Group_lst_cat, cov, df)

res$Liner <- ifelse(res$P.value < 0.05, '*', ' ')

fwrite(res, 'results/s0_Association/Cox_Food_group_trend_P.csv')
