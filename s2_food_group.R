rm(list = ls())
path1 <- "~/Desktop/Diet_aging"
path2 <- "C:/Users/Administrator/Desktop/Diet_aging"
if (dir.exists(path1)) {
  setwd(path1)
} else if (dir.exists(path2)) {
  setwd(path2)
} else {
  stop("路径不存在，请检查路径")
}

data_dir1 <- '/Users/mrli/Desktop/UKB_data'
data_dir2 <- 'C:/Users/Administrator/Desktop/UKB_data'
if (dir.exists(data_dir1)) {
  data_dir <- data_dir1
} else if (dir.exists(data_dir2)) {
  data_dir <- data_dir2
} else {
  stop("数据目录不存在，请检查路径")
}

library(data.table)
food_code <- fread(paste0(data_dir, "/Food/raw/Food_code.csv"))
Food_Group_lst <- unique(food_code$Food_Group)
Food_Group_lst_cat <- paste0(Food_Group_lst,'_Category')

# 'data/Food/Food_group_category.csv')
fg <- fread("data/Food/Food_group_category.csv")
colnames(fg)

# "data/Y_total.csv",
df <- fread("data/Y_total.csv")
colnames(df)

df <- merge(df, fg, by = 'eid', all = F)
component <- Food_Group_lst_cat[2]
ref = 2
y <- c('death','death_duration')
library(survival)
library(dplyr)

func_cox <- function(component, ref, cov){
  cat("Processing ", component, " with ref =", ref, "\n")
  tmp_data <- select(df, eid, y, all_of(component), all_of(cov))
  setnames(tmp_data, component, "pheno")

  levs <- levels(factor(tmp_data$pheno))
  lower_bound <- as.numeric(sub("\\s*\\[?\\(?([0-9\\.]+).*", "\\1", levs))
  levs_sorted <- levs[order(lower_bound)]
  levs <- levs_sorted
  
  n <- length(levs)
  ref_store <- ref
  if(ref > n){
    message(paste("ref =", ref, "over the number of levels =", n, "in", component))
    ref <- n
  }

  new_levels <- c(levs[ref], levs[-ref])
  tmp_data$pheno <- factor(tmp_data$pheno, levels = new_levels)
  events <- list()
  totals <- list()
  for (i in 1:n) {
    events[[i]] <- sum(tmp_data$death[tmp_data$pheno == new_levels[i]], na.rm = TRUE)
    totals[[i]] <- sum(!is.na(tmp_data$death[tmp_data$pheno == new_levels[i]]))
  }

  fmla <- as.formula(paste("Surv(death_duration, death) ~ pheno +", paste(cov, collapse = "+")))
  
  fit <- coxph(fmla, data = tmp_data)
  res <- cox.zph(fit)
  fit0 <- summary(fit)
  
  coef <- as.data.frame(fit0$coefficients)
  confint <- as.data.frame(fit0$conf.int)
  exposure <- sub("_Category", "", component)
  
  output <- data.frame(
    Ref = rep(ref_store, n),
    exposure = rep(exposure, n),
    outcome = rep("Death", n),
    events = unlist(events),
    totals = unlist(totals),
    class = levels(tmp_data$pheno),
    HR = c(1, coef$`exp(coef)`[1:(n-1)]),
    LCI = c(1, confint$`lower .95`[1:(n-1)]),
    HCI = c(1, confint$`upper .95`[1:(n-1)]),
    'HR(95%CI) '= c("Reference", paste0(round(coef$`exp(coef)`[1:(n-1)], 2), " (", 
                                        round(confint$`lower .95`[1:(n-1)], 2), "-", 
                                        round(confint$`upper .95`[1:(n-1)], 2), ")")),
    P = c(1, coef$`Pr(>|z|)`[1:(n-1)]),
    res.zph = rep(res$table[1,3], n)
  )
  
  return(output)
}


output <- data.frame()
cov1 <- c('Age','Sex','Ethnic_group')
cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ','energy','Vitamin','fast_time')
for (ref in 1:5) {
  for (i in Food_Group_lst_cat) {
    output0 <- func_cox(i,ref,cov)
    output <- rbind(output,output0)
  }
}
dir.create("results/s0_Association", showWarnings = FALSE, recursive = TRUE)


o1 <- subset(output, Ref == 1)[, -3]
o2 <- subset(output, Ref == 2)[, -3]
o3 <- subset(output, Ref == 3)[, -3]
o4 <- subset(output, Ref == 4)[, -3]
o5 <- subset(output, Ref == 5)[, -3]
out <- cbind(o1, o2, o3, o4, o5)
fwrite(out, "results/s0_Association/Cox_Food_group.csv")
fwrite(output, "results/s0_Association/output.csv")

output <- fread("results/s0_Association/output.csv")
colnames(output)
output$P <- p.adjust(output$P, method = "fdr")
output1 <- subset(output, P < 0.05) # /length(Food_Group_lst)
need_lst <- unique(output1$exposure)
paste0(need_lst, collapse = ", ")
writeLines(need_lst, "results/s0_Association/Need_Food_lst.txt")

cor_mat <- cor(fg[,..need_lst], method = "spearman", use = "pairwise.complete.obs")

dist_mat <- as.dist(1 - abs(cor_mat))  
hc <- hclust(dist_mat, method = "average")  
plot(hc, main = "Hierarchical clustering of food groups (Spearman)",
     xlab = "", sub = "", cex = 0.8)

threshold <- 0.4
abline(h = 1 - threshold, col = "red", lty = 2)
# save the plot

clusters <- cutree(hc, h = 1 - threshold)

need_lst1 <- tapply(names(clusters), clusters, function(x) x[1]) |> unname()

need_lst1
paste0(need_lst1, collapse = ", ")
writeLines(need_lst1, "results/s0_Association/Need_Food_lst.txt")
