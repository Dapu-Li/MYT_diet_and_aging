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


fg <- subset(fg, select = c('eid',score_col, 'MYT'))

Y <- fread("data/Y_total.csv")
df <- merge(Y, fg, by = 'eid',all = F)
cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level','Smoke','IPAQ',
         'energy','Vitamin','fast_time','PHQ4')

score_col

colnames(df)



suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(survival)
  library(survminer)
  library(magick)
})

y          <- "death"
y_duration <- "death_duration"
score_col  <- c("Poultry_score","Red_meats_score","Refined_grains_score",
                "Processed_meats_score","Nuts_score","Coffee_score",
                "Legumes_score","Butter_score","Sweetened_beverages_score",
                "Alcohol_score")


outdir <- "results/KM_plots_score_col"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df$death_duration <- df$death_duration/365.25
df$death <- as.numeric(df$death)
df$death_duration <- as.numeric(df$death_duration)

for (s in score_col) {
  df[[s]] <- as.numeric(df[[s]])
}


var <- score_col[1]
for (var in score_col) {
  new_var <- gsub("_score$", "", var)
  df[[new_var]] <- df[[var]]
  df[['score']] <- df[[var]]
  new_var <- 'score'
  formula_str <- sprintf("Surv(%s, %s) ~ %s", y_duration, y, new_var)
  fit <- survfit(as.formula(formula_str), data = df)

  pretty_title <- var |>
    stringr::str_replace("_score$", "") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_title()

  p <- ggsurvplot(
    fit, data = df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    legend.title = pretty_title,
    xlab = "Time",
    ylab = "Survival probability",
    title = pretty_title,
    ylim = c(0.85, 1)
  )
  p
  outfile <- file.path(outdir, paste0(var, ".png"))
  png(outfile, width = 2000, height = 2000, res = 300)
  print(p)
  dev.off()
}

files <- file.path(outdir, paste0(score_col, ".png"))
files <- files[file.exists(files)]

if (length(files) == 0) {
  stop("No single image files available. Please check score_col or data.")
}

imgs <- image_read(files)

montage <- image_montage(imgs, tile = "5x2", geometry = "900x900+1+1")

combined_file <- file.path(outdir, "KM_combined.png")
image_write(montage, path = combined_file, format = "png")

message("Individual plots saved in: ", outdir)
message("Combined plot saved: ", combined_file)
