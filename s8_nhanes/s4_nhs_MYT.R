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

fg$'Refined grains' <- fg$'Refined grains'*31.103
fg$'Processed meats' <- fg$'Processed meats'*31.103
fg$Nuts <- fg$Nuts*31.103
fg$Poultry <- fg$Poultry*31.103
fg$Red_meat <- fg$Red_meat *31.103

lst <- c('Poultry','Red_meat','Beer','Refined grains','Processed meats','Nuts',
         'Coffee','Legume','Butter','SSB','Wine')
group <- lst[1]
for (group in lst) {
  print(group)
  group_col <- group
  fg[[group_col]] <- as.numeric(fg[[group_col]])
  new_col <- paste0(group, "_Category")

  nb_zero <- sum(fg[[group_col]] == 0, na.rm = TRUE)
  nb_total <- nrow(fg)
  prec <- nb_zero / nb_total
  print(paste0("Zero: ", nb_zero, "/", nb_total, " (", round(prec, 3), ")"))

  fg[[new_col]] <- NA

  if (prec < 1/5) {
    fg[[new_col]] <- cut( fg[[group_col]],
                          breaks = unique(quantile(fg[[group_col]], probs = seq(0, 1, length.out = 6), na.rm = TRUE)),
                          include.lowest = TRUE, labels = FALSE )
  } else if (prec < 2/5) {
    fg[[new_col]] <- ifelse(fg[[group_col]] == 0, 0, NA)
    fg[[new_col]][fg[[group_col]] > 0] <- cut(
      fg[[group_col]][fg[[group_col]] > 0],
      breaks = unique(quantile(fg[[group_col]][fg[[group_col]] > 0], probs = seq(0, 1, length.out = 5), na.rm = TRUE)),
      include.lowest = TRUE, labels = FALSE )

  } else if (prec < 3/5) {
    fg[[new_col]] <- ifelse(fg[[group_col]] == 0, 0, NA)
    fg[[new_col]][fg[[group_col]] > 0] <- cut(
      fg[[group_col]][fg[[group_col]] > 0],
      breaks = unique(quantile(fg[[group_col]][fg[[group_col]] > 0], probs = seq(0, 1, length.out = 4), na.rm = TRUE)),
      include.lowest = TRUE, labels = FALSE )

  } else if (prec < 4/5) {
    fg[[new_col]] <- ifelse(fg[[group_col]] == 0, 0, NA)
    fg[[new_col]][fg[[group_col]] > 0] <- cut(
      fg[[group_col]][fg[[group_col]] > 0],
      breaks = unique(quantile(fg[[group_col]][fg[[group_col]] > 0], probs = seq(0, 1, length.out = 3), na.rm = TRUE)),
      include.lowest = TRUE, labels = FALSE )
  } else {
    fg[[new_col]] <- ifelse(fg[[group_col]] == 0, 0, 1)
  }

  if (prec < 1/5) {
    nb_groups <- 5
  } else if (prec < 2/5) {
    nb_groups <- 4
  } else if (prec < 3/5) {
    nb_groups <- 3
  } else if (prec < 4/5) {
    nb_groups <- 2
  } else {
    nb_groups <- 1
  }

  if (prec < 1/5) {
    grp_vals <- split(fg[[group_col]], fg[[new_col]])
  } else {
    grp_vals <- split(fg[[group_col]][fg[[group_col]] > 0], fg[[new_col]][fg[[group_col]] > 0])
  }

  labels <- sapply(grp_vals, function(x) {
    paste0(round(median(x, na.rm = TRUE), 1), " (",
           round(min(x, na.rm = TRUE), 1), ", ",
           round(max(x, na.rm = TRUE), 1), ")")
  })
  labels[1] <- sub("\\(", "[", labels[1])
  labels[1] <- sub("\\)", "]", labels[1])
  for (i in 2:(length(labels))) {
    labels[i] <- sub("\\)", "]", labels[i])
  }

  for (i in 2:length(labels)) {
    last_upper <- sub(".*,(.*)\\]", "\\1", labels[i - 1])
    last_upper <- trimws(last_upper)

    labels[i] <- sub("\\(([^,]+),", paste0("(", last_upper, ","), labels[i])
  }

  if (prec >= 1/5) {
    labels <- c("0" = "0", labels)
  }

  fg[[new_col]] <- factor(fg[[new_col]], labels = labels[as.character(sort(unique(fg[[new_col]])))])
  print(table(fg[[new_col]], useNA = "ifany"))
}

score_col <- paste0(lst, '_score')
score_col <- gsub(' ', '_', score_col)
score_col <- c(score_col, 'Alcohol_score')

fg$MYT <- rowSums(fg[, ..score_col])
table(fg$MYT)
fg$MYT_scale <- scale(fg$MYT)

colnames(fg)
fg <- subset(fg,select = c('SEQN','MYT','MYT_scale'))

fwrite(fg, 'data/nhanes/nhanes_MYT.csv')
