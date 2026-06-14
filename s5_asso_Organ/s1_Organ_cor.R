rm(list = ls())
# Set your working directory and data path here

library(data.table)

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


brain_age <- fread("data/brain_age_MRI/Brain_age_total.csv")
brain_age <- subset(brain_age, select = c('eid', 'Brain_PAD'))
brain_age$Brain_PAD <- scale(brain_age$Brain_PAD)
colnames(brain_age)[2] <- 'Brain_MRI_PAD'
organ_df1 <- merge(brain_age, organ_df, by = 'eid', all = F)

library(corrplot)
for (col in colnames(organ_df)[-1]) {
  organ_df[[col]] <- scale(organ_df[[col]])
}
cor_mat <- cor(subset(organ_df, select = -eid), use = "pairwise.complete.obs", method = "pearson")

cor_mat1 <- cor(subset(organ_df1, select = -eid), use = "pairwise.complete.obs", method = "pearson")
cor_mat1 <- as.data.frame(cor_mat1)
cor_mat1 <- subset(cor_mat1, select = Brain_MRI_PAD)
cor_mat <- as.data.frame(cor_mat)

fwrite(cor_mat, file = 'results/s3_Organ_age_NMED/Correlation_Plasma.csv',row.names = TRUE)
fwrite(cor_mat1, file = 'results/s3_Organ_age_NMED/Correlation_Brain_MRI_vs_Plasma.csv',row.names = TRUE)
