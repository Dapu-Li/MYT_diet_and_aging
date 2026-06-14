rm(list = ls())
# Set your working directory and data path here

LinearReg_4class <- function(x, y, covariates, data) {

  data[[x]] <- as.factor(data[[x]])

  formula <- paste(y, '~', paste(c(x, covariates), collapse = '+'))

  lr.t0 <- lm(as.formula(formula), data = data)

  lr.t1 <- summary(lr.t0)
  lr.t2 <- confint(lr.t0)

  coeffs <- c()
  p_vals <- c()

  for (i in 2:4) {
    lr.Coeff_B <- paste0(
      round(lr.t1$coefficients[i, 1], 2),
      ' (', round(lr.t2[i, 1], 2), ', ', round(lr.t2[i, 2], 2), ')'
    )
    lr.p_B <- round(lr.t1$coefficients[i, 4], 3)

    coeffs <- c(coeffs, lr.Coeff_B)
    p_vals <- c(p_vals, lr.p_B)
  }

  lr.coeff <- data.frame(
    'Characteristics' = y,
    'Predictors' = paste(x, levels(data[[x]])[2:4], sep = '_'),
    'Beta' = lr.t1$coefficients[2:4, 1],
    'Lower_CI' = lr.t2[2:4, 1],
    'Upper_CI' = lr.t2[2:4, 2],
    'Coef[95% CI]' = coeffs,
    'P value' = p_vals
  )

  return(lr.coeff)
}

res_4class <- function(Var.analys, Y, co, data) {
  res <- data.frame()

  for (y in Y) {
    for (i in 1) {
      combos <- combn(Var.analys, i)
      for (j in 1:ncol(combos)) {
        predictors <- combos[, j]
        res <- rbind(res, LinearReg_4class(predictors, y, co, data))
      }
    }
  }

  return(res)
}

LinearReg <- function(x, y, covariates, data) {

  print(paste("Running linear regression for", y, x))
  data <- data[!is.na(data[[x]]), ]
  data <- data[!is.na(data[[y]]), ]
  n <- nrow(data)
  formula <- paste(y, '~', paste(c(x, covariates), collapse = '+'))
  lr.t0 <- lm(as.formula(formula), data = data)

  lr.t1 <- summary(lr.t0)
  lr.t2 <- confint(lr.t0)

  beta <- paste0(lr.t1$coefficients[2, 1])

  lr.Coeff_B <- paste0(
    round(lr.t1$coefficients[2, 1], 2),
    ' (', round(lr.t2[2, 1], 2), ', ', round(lr.t2[2, 2], 2), ')'
  )
  lr.p_B <- lr.t1$coefficients[2, 4]

  lr.coeff <- data.frame(
    'Characteristics' = y,
    'Predictors' = c(x),
    'Beta' = c(beta),
    'Lower_CI' = c(lr.t2[2, 1]),
    'Upper_CI' = c(lr.t2[2, 2]),
    'Coef[95% CI]' = c(lr.Coeff_B),
    'P value' = c(lr.p_B),
    'N' = c(n)
  )

  return(lr.coeff)
}

res <- function(Var.analys, Y, co, data) {
  res <- data.frame()
  for (y in Y) {
    for (i in 1) {
      combos <- combn(Var.analys, i)
      for (j in 1:ncol(combos)) {
        predictors <- combos[, j]
        res <- rbind(res, LinearReg(predictors, y, co, data))
      }
    }
  }

  return(res)
}

library(data.table)
Myt <- fread("data/Food/MYT_score_either.csv")
Myt <- Myt[, .(eid, MYT)]
Myt$MYT_4cat <- cut(Myt$MYT,
                    breaks = quantile(Myt$MYT, probs = seq(0, 1, by = 0.25), na.rm = TRUE),
                    labels = c("Q1", "Q2", "Q3", "Q4"),include.lowest = TRUE)

table(Myt$MYT_4cat)
Q1 <- Myt[MYT_4cat == "Q1",]
Q2 <- Myt[MYT_4cat == "Q2",]
Q3 <- Myt[MYT_4cat == "Q3",]
Q4 <- Myt[MYT_4cat == "Q4",]
table(Q1$MYT)
table(Q2$MYT)
table(Q3$MYT)
table(Q4$MYT)

Myt$MYT_4cat <- factor(Myt$MYT_4cat, labels = c("0-2", "3", "4", "5+"))
table(Myt$MYT_4cat)

HEI <- fread("data/Food/HEI_2010.csv")
HEI$HEI_5cat <- cut(HEI$HEI,
                    breaks = quantile(HEI$HEI, probs = seq(0, 1, by = 0.25), na.rm = TRUE),
                    labels = c("Q1", "Q2", "Q3", "Q4"),include.lowest = TRUE)
HEI$HEI <- scale(HEI$HEI_2010)

df <- fread("data/Y_total.csv")
df <- merge(df, Myt, by = 'eid', all = F)
df <- merge(df, HEI, by = 'eid', all = F)
organ_df <- data.frame()
organ_dir <- 'data/Organ_age_NMED'
organ_lst <- c('Brain', 'Artery', 'Liver', 'Immune', 'Intestine',
            'Lung', 'Heart', 'Pancreas', 'Muscle', 'Adipose', 'Kidney','Organismal','Conventional')
organ <- 'Brain'
for (organ in organ_lst) {
  file_name <- paste0(organ_dir, '/', organ, '/', organ, '_age_total.csv')
  temp <- fread(file_name)
  temp[, paste0(organ, '_PAD') := get(paste0(organ, '_age')) - get('age')]
  temp <- subset(temp, select = c('eid', paste0(organ, '_PAD')))
  if (nrow(organ_df) == 0) {
    organ_df <- temp
  } else {
    organ_df <- merge(organ_df, temp, by = 'eid', all = F)
  }
}
Y <- paste0(organ_lst, '_PAD')
cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ','energy','Vitamin','fast_time')
df <- merge(df, organ_df, by = 'eid', all = F)

res_continuous <- res(
  Var.analys = c('MYT'),
  Y = Y,
  co = cov,
  data = df
)
res_categorical <- res_4class(
  Var.analys = c('MYT_4cat'),
  Y = Y,
  co = cov,
  data = df
)

brain_age <- fread("data/brain_age_MRI/Brain_age_total.csv")
brain_age <- subset(brain_age, select = c('eid', 'Brain_PAD'))
colnames(brain_age) <- c('eid', 'Brain_PAD_MRI')
df <- merge(df, brain_age, by = 'eid', all = F)
Y <- c('Brain_PAD_MRI')
res_brain <- res(
  Var.analys = c('MYT'),
  Y = Y,
  co = cov,
  data = df
)
res_brain_categorical <- res_4class(
  Var.analys = c('MYT_4cat'),
  Y = Y,
  co = cov,
  data = df
)

res_con <- rbind(res_continuous, res_brain)
res_cat <- rbind(res_categorical, res_brain_categorical)

fwrite(res_con, 'results/s6_Organ_asso/Organ_age_MYT_continuous.csv', row.names = FALSE)
fwrite(res_cat, 'results/s6_Organ_asso/Organ_age_MYT_4class.csv', row.names = FALSE)
