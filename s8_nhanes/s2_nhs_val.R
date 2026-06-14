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
library(data.table)
fg <- fread('data/nhanes/nhanes_food3.csv')
colnames(fg)

lst <- c('Poultry','Red_meat','Beer','Refined grains','Processed meats','Nuts',
         'Coffee','Legume','Butter','SSB','Wine')

score_col <- paste0(lst, '_score')
score_col <- gsub(' ', '_', score_col)
fg$MYT <- rowSums(fg[, ..score_col])
table(fg$MYT)

load('data/nhanes/nhanes_1999_2018.RData')
load('data/nhanes/zhirong.rdata')
load('data/nhanes/pattern.rdata')
colnames(zhirong)
pattern <- subset(zhirong, select = c('seqn','DASH.Mellen','dii','CDAI'))
colnames(pattern) <- c('seqn','DASH','DII','CDAI')

need_col <- c('seqn','ageyr','sex','eth1','BMI_kg.m2','poverty','edu','energy_kcal',
              'hei2015_total_score','mortstat','permth_exm')
d <- d[, need_col]
need_col <- c('seqn','smoke','PA_total_MET')
dd <- pa[, need_col]

dd <- merge(d,dd, by='seqn', all.y=TRUE)
colnames(fg)[1] <- 'seqn'
df <- merge(dd, fg, by='seqn', all.y=TRUE)
df <- df[!is.na(df$mortstat), ]
cov <- c('ageyr','sex','eth1','BMI_kg.m2','poverty','edu','energy_kcal')
na_df <- data.frame(
  Variable = names(df),
  Missing_Count = sapply(df, function(x) sum(is.na(x))),
  Missing_Percentage = sapply(df, function(x) sum(is.na(x)) / length(x) * 100)
)
na_df <- na_df[na_df$Missing_Count > 0, ]
na_df$Variable
cata_col <- c('PA_total_MET')
cata_col1 <- c('smoke','edu')
num_col <- c('BMI_kg.m2','poverty','BMI')

for (col in cata_col) {
  df[[col]][is.na(df[[col]])] <- "Unknown"
}
for (col in cata_col1) {
  mode_value <- names(sort(table(df[[col]]), decreasing = TRUE))[1]
  df[[col]][is.na(df[[col]])] <- mode_value
}
for (col in num_col) {
  df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
}
colnames(df)
table(df$mortstat)
df$death <- ifelse(df$mortstat == 'Assumed deceased', 1, 0)
df$MYT_scale <- scale(df$MYT)
df$HEI_scale <- scale(df$hei2015_total_score)
df$edu <- trimws(tolower(df$edu))
table(df$edu)
df$edu2 <- ifelse(df$edu %in% c("more than high school", "high school graduate",
                                "some college or aa degree",'college graduate or above'),
                  ">= High school",'< High school')
df <- merge(df, pattern, by='seqn', all.x=TRUE)
df$DASH_scale <- scale(df$DASH)
df$DII_scale <- scale(df$DII)
df$CDAI_scale <- scale(df$CDAI)
library(survival)
cov <- c('ageyr','sex','eth1','BMI_kg.m2','poverty','edu2','energy_kcal')
x <- c('HEI_scale','MYT_scale','DASH_scale','DII_scale','CDAI_scale')
res_val <- res_cox(time = 'permth_exm', status = c('death'),
                   Var.analys = x, co = cov, data = df)
