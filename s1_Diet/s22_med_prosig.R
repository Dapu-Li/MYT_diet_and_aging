rm(list = ls())
# Set your working directory and data path here

library(data.table)
df <- fread('data/Y_total.csv')
MYT <- fread('data/Food/MYT_score_either.csv')
MYT$MYT_scale <- scale(MYT$MYT)
df <- merge(df,MYT[,c('eid','MYT_scale')],by='eid')
MYT_sig <- fread('Results/s5_Lasso/MYT_Pro_signature.csv')
df <- merge(df, MYT_sig, by = 'eid', all.y = TRUE)
colnames(df)
df <- na.omit(df)
cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ','energy','Vitamin','fast_time')

x <- 'MYT_scale';m <- 'MYT_Pro_signature';y <- 'death'

library(broom)

var <- x
get_result <- function(model, var){
  tidy_mod <- tidy(model, conf.int = TRUE, exponentiate = inherits(model, "glm"))
  res <- tidy_mod[tidy_mod$term == var, ]
  data.frame(
    Variable = var,
    'OR,Beta' = round(res$estimate, 3),
    CI = paste0("(", round(res$conf.low, 3), ", ", round(res$conf.high, 3), ")"),
    Pvalue = signif(res$p.value, 3)
  )
}

fit1 <- lm(as.formula(paste(m, "~", x, "+", paste(cov, collapse = "+"))), data = df)
summary(fit1)
res1 <- get_result(fit1, x)

fit2 <- glm(as.formula(paste(y, "~", m, "+", paste(cov, collapse = "+"))),
            data = df, family = binomial)
summary(fit2)
exp(coef(fit2))
res2 <- get_result(fit2, m)

fit3 <- glm(as.formula(paste(y, "~", x, "+", paste(cov, collapse = "+"))),
            data = df, family = binomial)
summary(fit3)
res3 <- get_result(fit3, x)

final_res <- rbind(
  cbind(Path = "x to m", res1),
  cbind(Path = "m to y", res2),
  cbind(Path = "x to y", res3)
)
fwrite(final_res, 'Results/s5_Lasso/MYT_Pro_mediation_results.csv')
