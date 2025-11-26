rm(list = ls())
setwd("~/Desktop/Diet_aging")

# Cox
CoxReg <- function(time, status, x, co, data) {
  
  cat(paste("正在进行", status, x, "的Cox回归分析\n"))
  
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
  
  # 提取 p 值
  cox.p_B <- cox_summary$coefficients[1, 5]
  
  # 创建结果数据框
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
fg <- fread("data/Food/Food_group_category.csv")
imp <- fread('results/s1_Prediction/Food_Importance.csv')
imp <- imp[order(imp$TotalGain_cv, decreasing = T),]
colnames(imp)
imp_lst <- imp$Food
imp_lst <- gsub("_"," ",imp_lst)
#imp_lst <- imp_lst[imp_lst != 'Beer']

#/Users/mrli/Desktop/Diet_aging/results/s0_Association/food_group_paint1.pdf
proper <- fread('results/s0_Association/Cox_Food_group_pro.csv')
proper <- subset(proper,HR == 1)
proper <- subset(proper,exposure %in% imp_lst)
colnames(proper)
proper <- subset(proper,select = c("exposure","class"))
fwrite(proper,'data/Food/proper.csv',row.names = F)
diet_df <- data.frame()
table(fg$`Butter_Category`)
table(diet_df$Butter)
for (food in imp_lst){
  print(food)
  food_cat <- paste0(food,'_Category')
  tmp_df <- subset(fg,select = c("eid",food_cat))
  food_proper <- proper$class[which(proper$exposure == food)]
  tmp_df[[food]] <- ifelse(tmp_df[[food_cat]] %in% food_proper,1,0)
  tmp_df <- subset(tmp_df,select = c("eid",food))
  if (ncol(diet_df) == 0){
    diet_df <- tmp_df
  } else {
    diet_df[[food]] <- tmp_df[[food]]
  }
}


pattern_df <- subset(diet_df,select = c('eid',imp_lst[1]))
colnames(pattern_df)[2] <- 'Imp1'
for (i in 2:length(imp_lst)){
  print(i)
  col <- imp_lst[1:i]
  tmp_df <- subset(diet_df,select = c('eid',col))
  name <- paste0('Imp',i)
  pattern_df[[name]] <- rowSums(subset(tmp_df,select = col))
}
table(pattern_df$Imp11)
# scale
pattern_df[, (2:ncol(pattern_df)) := lapply(.SD, scale), .SDcols = 2:ncol(pattern_df)]

Y_total <- fread('data/Y_total.csv')
data <- merge(pattern_df,Y_total,by = 'eid',all = F)
cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ','energy','Vitamin','fast_time')
colnames(data)
library(survival)
res <- res_cox(time = 'death_duration', status = c('death'), Var.analys = colnames(pattern_df)[2:ncol(pattern_df)], co = cov, data = data)

fwrite(res,'results/s0_Association/Cox_Imp_food_pattern.csv',row.names = F)
