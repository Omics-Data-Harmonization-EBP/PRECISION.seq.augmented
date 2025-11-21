rm(list = ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(umap)
library(tidyverse)
library(gt)
library(viridis)
library(bstfun)
library(cowplot)
library(RUVSeq)
set.seed(42)

## Figure 2A
load("MSKpair_60000_count.RData")
load("precision.data_v2.RData")
augmented_benchmark_data <- do.call(cbind, lapply(benchmark[1:30], function(x) x$data))
augmented_benchmark_labe <- do.call(c, lapply(benchmark[1:30], function(x) x$label))
real_benchmark_data <- data.benchmark
real_benchmark_labe <- data.group
benchmark_eval_data <- log2(t(cbind(augmented_benchmark_data, real_benchmark_data)) + 1)
benchmark_eval_labe <- data.frame(
  label = c(augmented_benchmark_labe, real_benchmark_labe),
  group = c(rep("Augmented", length(augmented_benchmark_labe)), rep("Real", length(real_benchmark_labe)))
)
umap_result <- umap(benchmark_eval_data)
umap_df <- data.frame(
  UMAP1 = umap_result$layout[,1],
  UMAP2 = umap_result$layout[,2],
  label = benchmark_eval_labe$label,
  group = benchmark_eval_labe$group
)
fig2A <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = label, alpha = group)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("MXF" = "blue", "PMFH" = "red")) +
  scale_alpha_manual(values = c("Augmented" = 0.1, "Real" = 1)) +
  labs(title = "",
       x = "UMAP1", y = "UMAP2") +
  theme_minimal() +
  guides(color = guide_legend(title = ""),
         alpha = guide_legend(title = ""))

## Figure 2B
load("MSKpair_60000_count.RData")
load("precision.data_v2.RData")
calculate_mean_sd <- function(data) {
  data <- log2(data + 1)
  means <- apply(data, 1, mean)
  sds <- apply(data, 1, sd)
  summary <- data.frame(gene = rownames(data), mean = means, sd = sds)
  return(summary)
}
true_summary <- calculate_mean_sd(data.benchmark) %>%
  mutate(type = 'Real')
simulated_summaries <- lapply(seq_along(benchmark[1:10]), function(i) {
  data <- benchmark[[i]]$data
  rownames(data) <- rownames(data.benchmark)
  calculate_mean_sd(data) %>%
    mutate(type = "Simulated")
})
all_data <- bind_rows(true_summary, do.call(rbind, simulated_summaries))
fig2B <- ggplot(all_data, aes(x = mean, y = sd, color = type, alpha = type)) +
  geom_point() +
  scale_alpha_manual(values = c('Real' = 1, 'Simulated' = 0.2)) +
  scale_color_manual(values = c('Real' = "#fb8d62", 'Simulated' = "#66c2a5")) +
  labs(
    title = '',
    x = 'Mean',
    y = 'Standard Deviation'
  ) +
  theme_minimal() +
  theme(legend.title = element_blank())


## Figure 3
result_output <- function(scenario){
  array_data <- abind::abind(scenario, along = 3)
  average_df <- apply(array_data, c(1, 2), mean)
  rownames(average_df) <- c("HC-E", "HC-P", "HC-S", "MNM", "Kmeans",
                            "SOM", "PAM-E", "PAM-P", "PAM-S")
  colnames(average_df) <- c('Raw','TC','UQ','Med','TMM','DESeq','PoissonSeq','QN','RUVr','RUVs','RUVg')
  average_df <- round(average_df, 2)
  average_df[average_df < 0] <- 0
  # Reorder the rows (clustering methods)
  row_order <- c("HC-E", "HC-P", "HC-S", "PAM-E", "PAM-P", "PAM-S", "MNM", "Kmeans", "SOM")
  average_df <- average_df[row_order, ]
  # Reorder the RUV columns (normalization methods)
  ruv_cols <- c('RUVg', 'RUVr', 'RUVs')
  non_ruv_cols <- setdiff(colnames(average_df), c('RUVr', 'RUVs', 'RUVg'))
  average_df <- average_df[, c(non_ruv_cols, ruv_cols)]
  fig <- data.frame(average_df) %>%
    tibble::rownames_to_column("Cluster") %>%
    gt() %>%
    cols_label(Cluster = "") %>%
    data_color(
      columns = -Cluster,
      colors = scales::col_numeric(
        palette = RColorBrewer::brewer.pal(9, "Blues"),
        domain = c(0, 1)
      )
    ) %>%
    cols_align(align = "center") %>%
    bstfun::as_ggplot()
  return(fig)
}
load("clustering_sceanrio1.RData")
clean_res_c0.2_ari <- lapply(1:300, function(i) clean_results_list$clean_res_c0.2$indexessummary[[i]]$ari)
gt1 <- result_output(clean_res_c0.2_ari)
clean_res_c0.6_ari <- lapply(1:300, function(i) clean_results_list$clean_res_c0.6$indexessummary[[i]]$ari)
gt2 <- result_output(clean_res_c0.6_ari)
clean_res_c1_ari <- lapply(1:300, function(i) clean_results_list$clean_res_c1$indexessummary[[i]]$ari)
gt3 <- result_output(clean_res_c1_ari)
fig3 <- plot_grid(
  gt3, gt2, gt1,
  labels = c("(A)", "(B)", "(C)"),
  label_size = 8,
  nrow = 1
)

## Figure 4
load("clustering_sceanrio2.RData")
dirty_res_c0.6d1.2_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c0.6d1.2$indexessummary[[i]]$ari)
dirty_res_c0.6d1.2_fig <- result_output(dirty_res_c0.6d1.2_ari)
dirty_res_c0.6d1.5_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c0.6d1.5$indexessummary[[i]]$ari)
dirty_res_c0.6d1.5_fig <- result_output(dirty_res_c0.6d1.5_ari)
dirty_res_c0.6d2_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c0.6d2$indexessummary[[i]]$ari)
dirty_res_c0.6d2_fig <- result_output(dirty_res_c0.6d2_ari)
dirty_res_c1d1.2_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c1d1.2$indexessummary[[i]]$ari)
dirty_res_c1d1.2_fig <- result_output(dirty_res_c1d1.2_ari)
dirty_res_c1d1.5_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c1d1.5$indexessummary[[i]]$ari)
dirty_res_c1d1.5_fig <- result_output(dirty_res_c1d1.5_ari)
dirty_res_c1d2_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c1d2$indexessummary[[i]]$ari)
dirty_res_c1d2_fig <- result_output(dirty_res_c1d2_ari)
dirty_res_c1.5d1.2_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c1.5d1.2$indexessummary[[i]]$ari)
dirty_res_c1.5d1.2_fig <- result_output(dirty_res_c1.5d1.2_ari)
dirty_res_c1.5d1.5_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c1.5d1.5$indexessummary[[i]]$ari)
dirty_res_c1.5d1.5_fig <- result_output(dirty_res_c1.5d1.5_ari)
dirty_res_c1.5d2_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c1.5d2$indexessummary[[i]]$ari)
dirty_res_c1.5d2_fig <- result_output(dirty_res_c1.5d2_ari)
dirty_res_c2d1.2_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c2d1.2$indexessummary[[i]]$ari)
dirty_res_c2d1.2_fig <- result_output(dirty_res_c2d1.2_ari)
dirty_res_c2d1.5_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c2d1.5$indexessummary[[i]]$ari)
dirty_res_c2d1.5_fig <- result_output(dirty_res_c2d1.5_ari)
dirty_res_c2d2_ari <- lapply(1:300, function(i) dirty_results_list$dirty_res_c2d2$indexessummary[[i]]$ari)
dirty_res_c2d2_fig <- result_output(dirty_res_c2d2_ari)
fig4 <- plot_grid(
  dirty_res_c2d2_fig,
  dirty_res_c2d1.5_fig,
  dirty_res_c2d1.2_fig,
  
  dirty_res_c1.5d2_fig,
  dirty_res_c1.5d1.5_fig,
  dirty_res_c1.5d1.2_fig,
  
  dirty_res_c1d2_fig,
  dirty_res_c1d1.5_fig,
  dirty_res_c1d1.2_fig,
  
  dirty_res_c0.6d2_fig,
  dirty_res_c0.6d1.5_fig,
  dirty_res_c0.6d1.2_fig,
  
  labels = c("(A-1)", "(A-2)", "(A-3)", 
             "(B-1)", "(B-2)", "(B-3)",
             "(C-1)", "(C-2)", "(C-3)",
             "(D-1)", "(D-2)", "(D-3)"),
  label_size = 8,
  nrow = 4
)

# Figure 5
load("clustering_sceanrio1.RData")
load("clustering_sceanrio2.RData")
load("clustering_sceanrio3.RData")
load("clustering_sceanrio4.RData")
load("clustering_sceanrio5.RData")
extract_ari_matrix <- function(result_list, scenario_name, param_name) {
  ari_data   <- result_list[[param_name]]
  array_data <- abind::abind(
    lapply(seq_len(300), \(i) ari_data$indexessummary[[i]]$ari),
    along = 3
  )
  average_df <- apply(array_data, c(1, 2), mean)
  rownames(average_df) <- c("HC-E","HC-P","HC-S","MNM",
                            "Kmeans","SOM","PAM-E","PAM-P","PAM-S")
  colnames(average_df) <- c("Raw","TC","UQ","Med","TMM",
                            "DESeq","PoissonSeq","QN",
                            "RUVr","RUVs","RUVg")
  average_df |>
    as.data.frame() |>
    rownames_to_column("clustering_method") |>
    pivot_longer(-clustering_method,
                 names_to  = "normalization",
                 values_to = "ari") |>
    mutate(scenario   = scenario_name,
           parameters = param_name)
}
clean_c0.6     <- extract_ari_matrix(clean_results_list,
                                     "Biological",   "clean_res_c0.6")
dirty_c1d1.5   <- extract_ari_matrix(dirty_results_list,
                                     "Real_world",   "dirty_res_c1.5d2")
misspecify_k6  <- extract_ari_matrix(misspecify_results_list,
                                     "Mis_k6",       "misspecify_res_k6")
imbalance_r0.2 <- extract_ari_matrix(balanceness_results_list,
                                     "Imbalance",    "balanceness_res_r0.2")
all_data <- bind_rows(clean_c0.6,
                      dirty_c1d1.5,
                      misspecify_k6,
                      imbalance_r0.2)
heatmap_data <- all_data |>
  unite("scenario_method", scenario, clustering_method, sep = "_") |>
  select(normalization, scenario_method, ari) |>
  pivot_wider(names_from = scenario_method, values_from = ari) |>
  rename(id = normalization)
clustering_methods <- c("HC-E","HC-P","HC-S",
                        "PAM-E","PAM-P","PAM-S",
                        "MNM","Kmeans","SOM")
anchor_hex <- c(
  "HC-E"  = "#F2B8B5",
  "HC-P"  = "#CFE4D0",
  "HC-S"  = "#D8D1E9",
  "PAM-E" = "#F5CFCB",
  "PAM-P" = "#D0DEE0",
  "PAM-S" = "#F9E6D3",
  "MNM"   = "#F8DEC3",
  "Kmeans"= "#C7DFF3",
  "SOM"   = "#E4D7C9"
)
cell_anchor <- colorspace::darken(anchor_hex, 0.35)
cell_palettes <- purrr::map(cell_anchor,
                            ~ colorRampPalette(c("#FFFFFF", .x))(6)) |>
  set_names(tolower(names(cell_anchor)))
header_anchor <- colorspace::darken(anchor_hex, 0.50)
header_palettes <- purrr::map(header_anchor, \(x) rep(x, 2)) |>
  set_names(paste0(tolower(names(header_anchor)), "_hdr"))
palettes <- c(cell_palettes, header_palettes)
scenarios       <- c("Biological", "Real_world", "Imbalance", "Mis_k6")
scenario_labels <- c("A",          "B",          "C",         "D")
column_info <- tibble(
  id      = "id",
  group   = "normalization",
  name    = "",
  geom    = "text",
  palette = NA_character_,
  options = list(list(fontface = "plain", angle = 0))
)
for (method in clustering_methods) {
  for (i in seq_along(scenarios)) {
    column_info <- add_row(
      column_info,
      id      = paste(scenarios[i], method, sep = "_"),
      group   = method,
      name    = scenario_labels[i],      # shows A/B/C/D
      geom    = "funkyrect",
      palette = tolower(method),
      # make the A/B/C/D plain, not italic, and non-angled
      options = list(list(fontface = "plain", angle = 0))
    )
  }
}
column_groups <- tibble(
  group   = c("normalization", clustering_methods),
  level1  = c("Normalization", clustering_methods),
  palette = c(NA_character_,   paste0(tolower(clustering_methods), "_hdr"))
)
fig5 <- funkyheatmap::funky_heatmap(
  data          = heatmap_data,
  column_info   = column_info,
  column_groups = column_groups,
  palettes      = palettes,
  scale_column  = FALSE,   # ARI already 0-1
  add_abc       = FALSE
)

# Figure 6
options(chromote.headless = "new")
load("brca_subsampled_results_with_reps_v4.RData")
result_output <- function(scenario){
  scenario <- scenario[c(1,2,3,7,8,9,4,5,6),c(1:8, 11, 9, 10)]
  
  rownames(scenario) <-  c("HC-E", "HC-P", "HC-S", "PAM-E", "PAM-P", "PAM-S",
                           "MNM", "Kmeans", "SOM")
  colnames(scenario) <- c("None", "TC", "UQ", "Med", "TMM", "DESeq", "PoissonSeq", 
                          "QN", "RUVg", "RUVr", "RUVs")
  scenario <- round(scenario, 2)
  scenario[scenario < 0] <- 0
  
  fig <- data.frame(scenario) %>%
    tibble::rownames_to_column("Cluster") %>%
    gt() %>%
    cols_label(Cluster = "") %>%
    data_color(
      columns = -Cluster,
      colors = scales::col_numeric(
        palette = RColorBrewer::brewer.pal(9, "Blues"),
        domain = c(0, 1)
      )
    ) %>%
    cols_align(align = "center") %>%
    bstfun::as_ggplot()
  return(fig)
}
avg_reps_matrix <- function(result_list, base_key, metric = c("ari","silhouette"), reps = 1:10) {
  metric <- match.arg(metric)
  rep_names <- sprintf("%s__rep%03d", base_key, reps)
  rep_names <- rep_names[rep_names %in% names(result_list)]
  if (length(rep_names) == 0L)
    stop("No matching replicates found for base_key = ", base_key)
  mats <- lapply(rep_names, function(nm) {
    m <- as.matrix(result_list[[nm]]$res[[metric]])
    storage.mode(m) <- "numeric"
    m
  })
  rn <- rownames(mats[[1]])
  cn <- colnames(mats[[1]])
  mats <- lapply(mats, function(m) m[rn, cn, drop = FALSE]) 
  arr <- array(NA_real_, dim = c(length(rn), length(cn), length(mats)),
               dimnames = list(rn, cn, rep_names))
  for (i in seq_along(mats)) arr[,,i] <- mats[[i]]
  
  avg <- apply(arr, c(1,2), function(x) mean(x, na.rm = TRUE))
  dimnames(avg) <- list(rn, cn)
  attr(avg, "replicates_used") <- rep_names
  avg
}
avg_s2_LB_late <- avg_reps_matrix(result_list, "S2_allLate__LumA_vs_Basal", metric = "ari", reps = 1:10)
avg_s2_LB_balanced <- avg_reps_matrix(result_list, "S2_eraBalanced__LumA_vs_Basal", metric = "ari", reps = 1:10)
avg_s3_LB_confA <- avg_reps_matrix(result_list, "S3_confA_L1Early_L2Late__LumA_vs_Basal", metric = "ari", reps = 1:10)
avg_s3_LB_confB <- avg_reps_matrix(result_list, "S3_confB_L1Late_L2Early__LumA_vs_Basal", metric = "ari", reps = 1:10)
fig6 <- plot_grid(
  result_output(avg_s2_LB_late),     result_output(avg_s2_LB_balanced),
  result_output(avg_s3_LB_confA),    result_output(avg_s3_LB_confB),
  labels = c("(A)", "(B)", "(C)", "(D)"),
  label_size = 12,
  ncol = 2,
  align = "hv"
)




