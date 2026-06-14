rm(list = ls())
# Set your working directory and data path here

library(data.table)
df <- fread("results/s0_Association/Cox_Food_group.csv")
food_code <- fread(paste0(data_dir, "/Food/raw/Food_code.csv"))
Food_Group_lst <- unique(food_code$Food_Group)
colnames(df)

total_df <- data.frame()

for (food in Food_Group_lst){
  print(food)
  df_sub <- df[which(df$exposure == food), ]
  sub1 <- df_sub[, c(1:11)]
  sub2 <- df_sub[, c(12:22)]
  sub3 <- df_sub[, c(23:33)]
  sub4 <- df_sub[, c(34:44)]
  sub5 <- df_sub[, c(45:55)]

  hr_lst1 <- sub1$HR
  hr_lst2 <- sub2$HR
  hr_lst3 <- sub3$HR
  hr_lst4 <- sub4$HR
  hr_lst5 <- sub5$HR

  if (all(hr_lst1 >= 1)){
    df_sub <- sub1
  } else if (all(hr_lst2 >= 1)){
    df_sub <- sub2
  } else if (all(hr_lst3 >= 1)){
    df_sub <- sub3
  } else if (all(hr_lst4 >= 1)){
    df_sub <- sub4
  } else if (all(hr_lst5 >= 1)){
    df_sub <- sub5
  } else {
    df_sub <- sub1
    print(paste0(food, " has mixed directions!"))
  }

  if (nrow(total_df) == 0){
    total_df <- df_sub
  } else {
    total_df <- rbind(total_df, df_sub)
  }
}

fwrite(total_df,'results/s0_Association/Cox_Food_group_pro.csv')
library(ggplot2)
library(forestploter)
library(grid)

df <- total_df
colnames(df)
df$event_total <- paste0(df$events, "/", df$totals)
df$events <- NULL
df$totals <- NULL
df$res.zph <- NULL
df$Ref <- NULL
food <- Food_Group_lst[1]
paint_df <- data.frame()
for (food in Food_Group_lst){
  df_sub <- df[which(df$exposure == food), ]
  colnames(df_sub)

  tmp_df <- data.frame(matrix(ncol = ncol(df_sub), nrow = 0))
  colnames(tmp_df) <- colnames(df_sub)
  tmp_df <- tmp_df[,c(1,3:ncol(tmp_df))]
  tmp_df[1,] <- c(food, '', '', '', '', '', '')
  df_sub <- subset(df_sub, select = -c(exposure))
  df_sub[[1]] <- paste0("   ", df_sub[[1]])

  colnames(df_sub)[1] <- 'exposure'
  df_sub <- rbind(tmp_df, df_sub)
  if (nrow(paint_df) == 0){
    paint_df <- df_sub
  } else {
    paint_df <- rbind(paint_df, df_sub)
  }
}
library(forestploter)

dt <- paint_df

dt$` ` <- paste(rep(" ", 20), collapse = " ")

dt$HR  <- as.numeric(dt$HR)
dt$LCI <- as.numeric(dt$LCI)
dt$HCI <- as.numeric(dt$HCI)
dt$P   <- as.numeric(dt$P)

dt$P_label <- ifelse(dt$P < 0.001, "<0.001", round(dt$P, 3))

dt1 <- dt

dt1$group <- ifelse(dt1$P < 0.05, "p<0.05", "p>=0.05")

dt1[is.na(dt1)] <- ""
tm <- forest_theme(
  base_size = 11,
  ci_pch = 15,
  ci_col = c("p<0.05" = "red", "p>=0.05" = "black"),
  ci_lwd = 1.2,
  refline_col = "grey40",
  vert_line_col = "grey70",
  vert_line_lty = 2,
  legend_name = "P value",
  legend_value = c("p<0.05", "p>=0.05"),
  legend_col = c("red", "black")
)

dt1$HR  <- as.numeric(dt1$HR)
dt1$LCI <- as.numeric(dt1$LCI)
dt1$HCI <- as.numeric(dt1$HCI)
dt1$P   <- as.numeric(dt1$P)
p <- forest(
  dt1[, c(1, 7, 5, 8, 9)],
  est       = dt1$HR,
  lower     = dt1$LCI,
  upper     = dt1$HCI,
  ci_column = 4,
  ref_line  = 1,
  xlim      = c(0.94, 1.3),
  ticks_at  = c(0.94,1, 1.3),
  footnote  = "",
  theme     = tm,
  nudge_y   = 0.2,
  group     = dt1$group
)

ggsave(
  filename = "results/s0_Association/food_group_paint1.pdf",
  plot     = p,
  width    = 8,
  height   = 45
)
