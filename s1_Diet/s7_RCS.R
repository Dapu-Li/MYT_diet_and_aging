rm(list = ls())
# Set your working directory and data path here

library(data.table)
library(rcssci)
packageVersion("rcssci")
food_code <- fread(paste0(data_dir, "/Food/raw/Food_code.csv"))
Food_Group_lst <- unique(food_code$Food_Group)
Food_Group_lst_cat <- paste0(Food_Group_lst,'_Category')

fg <- fread("data/Food/Food_group_category.csv")[, c('eid', Food_Group_lst), with = F]
colnames(fg)


df <- fread("data/Y_total.csv")
colnames(df)
df <- merge(df, fg, by = 'eid', all = F)
cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ',
         'energy','Vitamin','fast_time','PHQ4')

dir.create('results/s0_Association/RCS', recursive = T, showWarnings = F)

for (food in Food_Group_lst) {
  cat(food, '\n')
  path_dir <- paste0('results/s0_Association/RCS/', food)
  if (!dir.exists(path_dir)) {
    dir.create(path_dir, recursive = T)
  }
  rcssci_cox(data = df, time = 'death_duration', y = 'death',
             covariates = cov,
             x = food,
             filepath = path_dir)
}
