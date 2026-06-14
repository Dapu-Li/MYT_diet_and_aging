rm(list = ls())
# Set your working directory and data path here

library(data.table)

fg <- fread('data/Food/Food_group.csv')
food_code <- fread(paste0(data_dir, "/Food/raw/Food_code.csv"))

colnames(food_code)
Food_Group_lst <- unique(food_code$Food_Group)

for (group in Food_Group_lst) {
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

dir.create('data/Food/', showWarnings = FALSE, recursive = TRUE)
fwrite(fg, 'data/Food/Food_group_category_new.csv')
