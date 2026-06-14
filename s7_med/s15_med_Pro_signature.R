rm(list = ls())
# Set your working directory and data path here

library(data.table)
library(mediation)

df <- fread('data/Y_total.csv')
Pro <- fread('results/s5_Lasso/SSB_Pro_signature.csv')
df <- merge(df, Pro, by = 'eid', all.x = TRUE)

cov2 <- c("Age", "Sex", "Ethnic_group", "TDI", "education_level",'bmi_catagory',
          "Smoke", "Alcohol_freq", "IPAQ", "PHQ4", "energy")
xy <- c("SSB_scale", "depression")

x_m <- fread('results/s0_Association/SSB_pro.csv')
pro_lst <- x_m$Pro_code[x_m$pval_fdr < 0.05]

cat_vars <- c("Sex", "Ethnic_group", "education_level", "Smoke", "Alcohol_freq", "IPAQ",'bmi_catagory')
df[, (cat_vars) := lapply(.SD, as.factor), .SDcols = cat_vars]

wd <- getwd()
dir.create(paste0(wd, "/results/s4_med_either/chunks"), showWarnings = FALSE)
dir.create(paste0(wd, "/results/s4_med_either/errors"), showWarnings = FALSE)

start_time <- Sys.time()
total_pro <- length(pro_lst)
mediation_results <- data.frame()
chunk_index <- 1
count <- 0
colnames(df)
time_log <- c()
pro_lst <- 'SSB_Pro_signature'
pro <- pro_lst[1]
set.seed(123)
for (pro in pro_lst) {
  count <- count + 1
  single_start <- Sys.time()

  tmp_df <- subset(df, select = c('eid', pro, cov2, xy))
  tmp_df <- tmp_df[!is.na(tmp_df[[pro]]), ]
  if (nrow(tmp_df) < 100) next

  colnames(tmp_df)[colnames(tmp_df) == pro] <- "Mediator_Pro"
  form_m <- as.formula(paste("Mediator_Pro ~ SSB_scale +", paste(cov2, collapse = " + ")))
  form_y <- as.formula(paste("depression ~ Mediator_Pro + SSB_scale +", paste(cov2, collapse = " + ")))

  model.m <- tryCatch(
    lm(form_m, data = tmp_df),
    error = function(e) {
      writeLines(e$message, paste0("results/s4_med_either/errors/Error_", pro, ".txt"))
      return(NULL)
    }
  )
  if (is.null(model.m)) next

  model.y <- tryCatch(
    glm(form_y, family = binomial(link = "logit"), data = tmp_df),
    error = function(e) {
      writeLines(e$message, paste0("results/s4_med_either/errors/Error_", pro, ".txt"))
      return(NULL)
    }
  )
  if (is.null(model.y)) next

  med <- tryCatch(
    mediate(
      model.m = model.m,
      model.y = model.y,
      treat = "SSB_scale",
      mediator = "Mediator_Pro",
      sims = 1000
    ),
    error = function(e) {
      writeLines(e$message, paste0("results/s4_med_either/errors/Error_", pro, ".txt"))
      return(NULL)
    }
  )
  if (is.null(med)) next
  summary(med)
  mediation_results <- rbind(mediation_results, data.frame(
    Protein = pro,
    ACME = med$d0,
    ACME_p = med$d0.p,
    ADE = med$z0,
    ADE_p = med$z0.p,
    Total_Effect = med$tau.coef,
    Total_Effect_p = med$tau.p,
    Prop_Mediated = med$n0,
    Prop_Mediated_p = med$n0.p,
    stringsAsFactors = FALSE
  ))

  single_end <- Sys.time()
  time_taken <- as.numeric(difftime(single_end, single_start, units = "secs"))
  time_log <- c(time_log, time_taken)

  avg_time <- mean(time_log)
  remaining_pro <- total_pro - count
  remaining_time <- remaining_pro * avg_time

  format_time <- function(sec) {
    if (sec > 3600) {
      sprintf("%.2f hours", sec / 3600)
    } else if (sec > 60) {
      sprintf("%.2f minutes", sec / 60)
    } else {
      sprintf("%.2f seconds", sec)
    }
  }

  cat(sprintf("Finished: %4d/%d | Protein: %-20s | Time: %.2f sec | Avg: %.2f sec | ETA: %s\n",
              count, total_pro, pro, time_taken, avg_time, format_time(remaining_time)))

  if (count %% 100 == 0) {
    chunk_path <- paste0("results/s4_med_either/chunks/Pro_signature_mediation_chunk_", chunk_index, ".csv")
    fwrite(mediation_results, chunk_path)
    mediation_results <- data.frame()
    chunk_index <- chunk_index + 1
    cat("Saved chunk", chunk_index - 1, "to", chunk_path, "\n")
  }
}


if (nrow(mediation_results) > 0) {
  chunk_path <- paste0("results/s4_med_either/chunks/Pro_signature_mediation_chunk_", chunk_index, ".csv")
  fwrite(mediation_results, chunk_path)
  cat("Saved final chunk", chunk_index, "\n")
}

all_chunks <- list.files(
  path = "results/s4_med_either/chunks",
  pattern = "^Pro_signature_mediation_chunk_.*\\.csv$",
  full.names = TRUE
)
all_results <- rbindlist(lapply(all_chunks, fread), use.names = TRUE, fill = TRUE)

fwrite(all_results, "results/s4_med_either/Pro_signature_mediation.csv", row.names = FALSE)
cat("All chunks combined and saved to Pro_mediation")
