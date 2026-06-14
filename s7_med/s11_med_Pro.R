rm(list = ls())
# Set your working directory and data path here

library(data.table)
library(mediation)

df <- fread('data/Y_total.csv')
Pro <- fread(paste0(data_dir, '/Protein/processed/Protein.csv'))
df <- merge(df, Pro, by = 'eid', all.x = TRUE)
myt <- fread('data/Food/MYT_score.csv')[,c('eid','MYT')]
myt$MYT_score <- scale(myt$MYT)
df <- merge(df,myt,by = 'eid',all.y = TRUE)

cov <- c('Age','Sex','Ethnic_group','BMI','TDI','education_level',
         'Smoke','IPAQ','energy','Vitamin','fast_time')
xy <- c("MYT_score", "death")

x_m <- fread('Results/s4_med/Score_Pro.csv')
m_y <- fread('Results/s4_med/Pro_death.csv')
xmy <- merge(x_m, m_y, by = 'Pro_code', all.x = TRUE)
colnames(xmy)
xmy_sig <- subset(xmy, pval_fdr.x < 0.05 & pval_fdr.y < 0.05)
pro_lst <- xmy_sig$Pro_code

cat_vars <- c("Sex", "Ethnic_group", "education_level", "Smoke", "IPAQ",'Vitamin')
df[, (cat_vars) := lapply(.SD, as.factor), .SDcols = cat_vars]

wd <- getwd()
dir.create(paste0(wd, "/results/s4_med/chunks"), showWarnings = FALSE)
dir.create(paste0(wd, "/results/s4_med/errors"), showWarnings = FALSE)

start_time <- Sys.time()
total_pro <- length(pro_lst)
mediation_results <- data.frame()
chunk_index <- 1
count <- 0

time_log <- c()
pro_lst1 <- pro_lst[1:10]
set.seed(123)
pro <- pro_lst1[1]
for (pro in pro_lst) {
  count <- count + 1
  single_start <- Sys.time()

  tmp_df <- subset(df, select = c('eid', pro, cov, xy))
  tmp_df <- tmp_df[!is.na(tmp_df[[pro]]), ]
  if (nrow(tmp_df) < 100) next

  colnames(tmp_df)[colnames(tmp_df) == pro] <- "Mediator_Pro"
  form_m <- as.formula(paste("Mediator_Pro ~ MYT_score +", paste(cov, collapse = " + ")))
  form_y <- as.formula(paste("death ~ Mediator_Pro + MYT_score +", paste(cov, collapse = " + ")))

  model.m <- tryCatch(
    lm(form_m, data = tmp_df),
    error = function(e) {
      writeLines(e$message, paste0("results/s4_med/errors/Error_", pro, ".txt"))
      return(NULL)
    }
  )
  if (is.null(model.m)) next

  model.y <- tryCatch(
    glm(form_y, family = binomial(link = "logit"), data = tmp_df),
    error = function(e) {
      writeLines(e$message, paste0("results/s4_med/errors/Error_", pro, ".txt"))
      return(NULL)
    }
  )
  if (is.null(model.y)) next

  med <- tryCatch(
    mediate(
      model.m = model.m,
      model.y = model.y,
      treat = "MYT_score",
      mediator = "Mediator_Pro",
      sims = 1000
    ),
    error = function(e) {
      writeLines(e$message, paste0("results/s4_med/errors/Error_", pro, ".txt"))
      return(NULL)
    }
  )
  if (is.null(med)) next

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

  if (count %% 1 == 0) {
    chunk_path <- paste0("results/s4_med/chunks/Pro_mediation_chunk_", pro, ".csv")
    fwrite(mediation_results, chunk_path)
    mediation_results <- data.frame()
    chunk_index <- chunk_index + 1
    cat("Saved chunk", chunk_index - 1, "to", chunk_path, "\n")
  }
}


if (nrow(mediation_results) > 0) {
  chunk_path <- paste0("results/s4_med/chunks/Pro_mediation_chunk_", pro, ".csv")
  fwrite(mediation_results, chunk_path)
  cat("Saved final chunk", chunk_index, "\n")
}

all_chunks <- list.files(
  path = "results/s4_med/chunks",
  pattern = "^Pro_mediation_chunk_.*\\.csv$",
  full.names = TRUE
)
all_results <- rbindlist(lapply(all_chunks, fread), use.names = TRUE, fill = TRUE)
colnames(all_results)
all_results$p_fdr <- p.adjust(all_results$Prop_Mediated_p, method = "fdr")
fwrite(all_results, "results/s4_med/Pro_mediation.csv", row.names = FALSE)
cat("All chunks combined and saved to Pro_mediation")
