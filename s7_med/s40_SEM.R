rm(list = ls())
# Set your working directory and data path here

library(data.table)

df <- fread('data/Y_total.csv')
MYT <- fread('data/Food/MYT_score_either.csv')
MYT$MYT_scale <- scale(MYT$MYT)
df <- merge(df,MYT[,c('eid','MYT_scale')],by='eid')
colnames(df)

x_Meta <- fread('Results/s4_med_either/Meta_summary.csv')
colnames(x_Meta)[1] <- 'Meta_code'
x_Meta <- x_Meta[order(x_Meta$p_fdr), ]
Meta_lst <- x_Meta$Meta_code[x_Meta$p_fdr < 0.05]

x_Infla <- fread('Results/s4_med_either/Infla_summary.csv')
colnames(x_Infla)[1] <- 'Infla_code'
x_Infla <- x_Infla[order(x_Infla$p_fdr), ]
Infla_lst <- x_Infla$Infla_code[x_Infla$p_fdr < 0.05]

Meta <- fread(paste0(data_dir, '/Metabolite/processed/meta.csv'))[, c("eid", Meta_lst), with = FALSE]
Meta[, 2:ncol(Meta) := lapply(.SD, log1p), .SDcols = 2:ncol(Meta)]
Meta[, 2:ncol(Meta) := lapply(.SD, scale), .SDcols = 2:ncol(Meta)]
Infla <- fread(paste0(data_dir, '/Inflammation/processed/Inflammation.csv'))[, c("eid", Infla_lst), with = FALSE]
Infla[, names(Infla) := lapply(.SD, function(x) {
  x[is.infinite(x)] <- NA
  return(x)
})]
Infla[, 2:ncol(Infla) := lapply(.SD, scale), .SDcols = 2:ncol(Infla)]

dd <- merge(df[, c('eid', 'MYT_scale', 'death')], Meta, by = 'eid', all.x = TRUE)
dd <- merge(dd, Infla, by = 'eid', all.x = TRUE)


names(dd) <- gsub("-", "_", names(dd))
names(dd) <- gsub(" ", "_", names(dd))
names(dd) <- gsub("%", "_", names(dd))
dd <- na.omit(dd)
Meta_lst <- gsub("-", "_", Meta_lst)
Meta_lst <- gsub("%", "_", Meta_lst)

library(lavaan)
model_meta <- paste0("Metabolites =~ ", paste(Meta_lst, collapse = " + "))
model_infla <- paste0("Inflammation =~ ", paste(Infla_lst, collapse = " + "))
model_struct <- "
  Metabolites ~ MYT_scale
  Inflammation ~ MYT_scale + Metabolites
  death ~ MYT_scale + Metabolites + Inflammation
"
sem_model <- paste(model_meta, model_infla, model_struct, sep="\n")

fit <- sem(sem_model, data = dd, optim.method = "BFGS")
lavInspect(fit, 'optim.gradient')


getOption("max.print")
options(max.print = 1e7)
summary(fit, standardized = TRUE, fit.measures = TRUE)
