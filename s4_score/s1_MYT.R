rm(list = ls())
# Set your working directory and data path here

library(data.table)

food_code <- fread(paste0(data_dir, "/Food/raw/Food_code.csv"))
Food_Group_lst <- unique(food_code$Food_Group)
fg <- fread("data/Food/Food_group_category.csv")
fg <- subset(fg, select = c('eid', Food_Group_lst))

lst <- c('Poultry','Red meats','Beer','Refined grains','Processed meats','Nuts',
         'Coffee','Legumes','Butter','Sweetened beverages','Wine')

score_col <- paste0(lst, '_score')
score_col <- gsub(' ', '_', score_col)
score_col <- c(score_col, 'Alcohol_score')

fg$MYT <- rowSums(fg[, ..score_col])

table(fg$MYT)
fg <- subset(fg, select = c('eid',lst,score_col, 'MYT'))
fwrite(fg, "data/Food/MYT_score_either.csv", row.names = F)
