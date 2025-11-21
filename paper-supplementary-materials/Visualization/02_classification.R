options(chromote.headless = "new")
library(tidyverse)
library(gt)
library(viridis)
library(bstfun)
library(cowplot)
library(PRECISION.seq.augmented)
library(edgeR)
library(umap)

## Figure 1
load("classification_scenario1_analysis.RData")
result_output <- function(scenario){
  array_data <- abind::abind(scenario, along = 3)
  average_df <- apply(array_data, c(1, 2), mean, na.rm = TRUE)
  average_df2 <- t(average_df)[,-c(6,10)]
  
  rownames(average_df2) <- c("LASSO", "PAM", "KNN", "SVM", "SVM-FS", "RF", "RF-FS")
  colnames(average_df2) <- c("None", "TC", "UQ", "Med", "TMM", "DESeq", "PoissonSeq", 
                             "QN", "RUVg", "RUVr", "RUVs")
  average_df2 <- round(average_df2, 2)
  average_df2[average_df2 < 0] <- 0
  
  fig <- data.frame(average_df2) %>%
    tibble::rownames_to_column("Classification") %>%
    gt() %>%
    cols_label(Classification = "") %>%
    data_color(
      columns = -Classification,
      colors = scales::col_numeric(
        palette = RColorBrewer::brewer.pal(9, "Blues"),
        domain = c(0, 1)
      )
    ) %>%
    cols_align(align = "center") %>%
    bstfun::as_ggplot()
  return(fig)
}
c_values <- c(1.0, 0.6, 0.2)
gt_plots <- list()
for (i in seq_along(c_values)) {
  name_pattern <- if(c_values[i] < 0) {
    sprintf("c-%.1f", abs(c_values[i]))
  } else {
    sprintf("c%.1f", c_values[i])
  }
  accuracies <- lapply(1:300, function(x) results[[name_pattern]][[x]]$test2$accuracy)
  gt_plots[[length(gt_plots) + 1]] <- result_output(accuracies)
}
fig1 <- plot_grid(
  gt_plots[[1]], 
  gt_plots[[2]],
  gt_plots[[3]],
  labels = c("(A)", 
             "(B)", 
             "(C)"),
  label_size = 12,
  nrow = 3
)

## Figure 2
load("classification_scenario2_analysis.RData")
result_output <- function(scenario){
  array_data <- abind::abind(scenario, along = 3)
  average_df <- apply(array_data, c(1, 2), mean)
  average_df2 <- t(average_df)[,-c(6,10)]
  rownames(average_df2) <- c("LASSO", "PAM", "KNN", "SVM", "SVM-FS", "RF", "RF-FS")
  colnames(average_df2) <- c("None", "TC", "UQ", "Med", "TMM", "DESeq", "PoissonSeq", 
                             "QN", "RUVg", "RUVr", "RUVs")
  average_df2 <- round(average_df2, 2)
  average_df2[average_df2 < 0] <- 0
  fig <- data.frame(average_df2) %>%
    tibble::rownames_to_column("Classification") %>%
    gt() %>%
    cols_label(Classification = "") %>%
    data_color(
      columns = -Classification,
      colors = scales::col_numeric(
        palette = RColorBrewer::brewer.pal(9, "Blues"),
        domain = c(0, 1)
      )
    ) %>%
    cols_align(align = "center") %>%
    bstfun::as_ggplot()
  return(fig)
}
c_values <- c(1.0, 0.6, 0.2)
d_values <- c(1.2, 1.5, 2.0)
gt_plots <- list()
for (i in seq_along(c_values)) {
  for (j in seq_along(d_values)) {
    name_pattern <- if(c_values[i] < 0) {
      sprintf("c-%.1f_d%.1f", abs(c_values[i]), d_values[j])
    } else {
      sprintf("c%.1f_d%.1f", c_values[i], d_values[j])
    }
    
    # Extract accuracies
    accuracies <- lapply(1:300, function(x) results[[name_pattern]][[x]]$test2$accuracy)
    gt_plots[[length(gt_plots) + 1]] <- result_output(accuracies)
  }
}
fig2 <- plot_grid(
  plotlist = gt_plots,
  labels = c("(A-1)", "(A-2)", "(A-3)", 
             "(B-1)", "(B-2)", "(B-3)",
             "(C-1)", "(C-2)", "(C-3)"),
  label_size = 10,
  nrow = 3
)

## Figure 3
load("classification_scenario2_analysis.RData")
result_output <- function(scenario){
  array_data <- abind::abind(scenario, along = 3)
  average_df <- apply(array_data, c(1, 2), mean)
  average_df2 <- t(average_df)[,c(5,6,9,10)]
  rownames(average_df2) <- c("LASSO", "PAM", "KNN", "SVM", "SVM-FS", "RF", "RF-FS")
  colnames(average_df2) <- c("TMM", "fTMM", "QN", "fQN")
  average_df2 <- round(average_df2, 2)
  average_df2[average_df2 < 0] <- 0
  fig <- data.frame(average_df2) %>%
    tibble::rownames_to_column("Classification") %>%
    gt() %>%
    cols_label(Classification = "") %>%
    data_color(
      columns = -Classification,
      colors = scales::col_numeric(
        palette = RColorBrewer::brewer.pal(9, "Blues"),
        domain = c(0, 1)
      )
    ) %>%
    cols_align(align = "center") %>%
    bstfun::as_ggplot()
  return(fig)
}
c_values <- c(1.0, 0.6, 0.2)
d_values <- c(1.2, 1.5, 2.0)
gt_plots <- list()
for (i in seq_along(c_values)) {
  for (j in seq_along(d_values)) {
    # Create the name pattern as in results$cv
    name_pattern <- if(c_values[i] < 0) {
      sprintf("c-%.1f_d%.1f", abs(c_values[i]), d_values[j])
    } else {
      sprintf("c%.1f_d%.1f", c_values[i], d_values[j])
    }
    
    # Extract accuracies
    accuracies <- lapply(1:300, function(x) results[[name_pattern]][[x]]$test2$accuracy)
    gt_plots[[length(gt_plots) + 1]] <- result_output(accuracies)
  }
}
fig3 <- plot_grid(
  plotlist = gt_plots,
  labels = c("(A-1)", "(A-2)", "(A-3)", 
             "(B-1)", "(B-2)", "(B-3)",
             "(C-1)", "(C-2)", "(C-3)"),
  label_size = 10,
  nrow = 3
)

## Figure 4
load("classification_scenario1_analysis.RData")
clean_result = results
load("classification_scenario2_analysis.RData")
dirty_result = results
extract_accuracy_matrix <- function(result_list, scenario_name, param_name) {
  accuracies <- lapply(1:300, function(x) result_list[[param_name]][[x]]$test2$accuracy)
  array_data <- abind::abind(accuracies, along = 3)
  average_df <- apply(array_data, c(1, 2), mean)
  average_df2 <- t(average_df)[,-c(6,10)]
  rownames(average_df2) <- c("LASSO", "PAM", "KNN", "SVM", "SVM-FS", "RF", "RF-FS")
  colnames(average_df2) <- c("None", "TC", "UQ", "Med", "TMM", "DESeq", "PoissonSeq", 
                             "QN", "RUVg", "RUVr", "RUVs")
  average_df2 |>
    as.data.frame() |>
    rownames_to_column("classification_method") |>
    pivot_longer(-classification_method,
                 names_to  = "normalization",
                 values_to = "accuracy") |>
    mutate(scenario   = scenario_name,
           parameters = param_name)
}
clean_c0.2     <- extract_accuracy_matrix(clean_result,
                                          "Biological",   "c0.2")
dirty_c0.2d1.5 <- extract_accuracy_matrix(dirty_result,
                                          "Real_world",   "c0.2_d1.5")
all_data <- bind_rows(clean_c0.2,
                      dirty_c0.2d1.5)
heatmap_data <- all_data |>
  unite("scenario_method", scenario, classification_method, sep = "_") |>
  select(normalization, scenario_method, accuracy) |>
  pivot_wider(names_from = scenario_method, values_from = accuracy) |>
  rename(id = normalization)
anchor_hex <- c(
  "LASSO"  = "#F2B8B5",  
  "PAM"    = "#CFE4D0",  
  "KNN"    = "#D8D1E9",  
  "SVM"    = "#F8DEC3",  
  "SVM-FS" = "#C7DFF3",  
  "RF"     = "#E4D7C9", 
  "RF-FS"  = "#F5CFCB"   
)
cell_anchor <- colorspace::darken(anchor_hex, 0.35)
cell_palettes <- purrr::map(cell_anchor,
                            ~ colorRampPalette(c("#FFFFFF", .x))(6)) |>
  set_names(tolower(names(cell_anchor)))
header_anchor <- colorspace::darken(anchor_hex, 0.50)
header_palettes <- purrr::map(header_anchor, \(x) rep(x, 2)) |>
  set_names(paste0(tolower(names(header_anchor)), "_hdr"))
palettes <- c(cell_palettes, header_palettes)
classification_methods <- names(anchor_hex)
scenarios          <- c("Biological", "Real_world")
scenario_labels    <- c("Bio",        "Real")
column_info <- tibble(
  id      = "id",
  group   = "normalization",
  name    = "",
  geom    = "text",
  palette = NA_character_,
  options = list(list())
)
for (method in classification_methods) {
  for (i in seq_along(scenarios)) {
    column_info <- add_row(
      column_info,
      id      = paste(scenarios[i], method, sep = "_"),
      group   = method,
      name    = scenario_labels[i],
      geom    = "funkyrect",
      palette = tolower(method),  # DATA cells use gradient palettes
      options = list(list())
    )
  }
}
column_groups <- tibble(
  group   = c("normalization", classification_methods),
  level1  = c("Normalization", classification_methods),
  palette = c(NA_character_,   paste0(tolower(classification_methods), "_hdr"))
)
fig4 <- funkyheatmap::funky_heatmap(
  data          = heatmap_data,
  column_info   = column_info,
  column_groups = column_groups,
  palettes      = palettes,
  scale_column  = FALSE,   
  add_abc       = FALSE
)

## Figure 5
load("classification_results_simulated_all_v5_1027.RData")
load("brca_LumA_vs_LumB_trainEarlyAB_testLate_with_reps.RData")  # loads `datasets`
stopifnot(
  exists("results"), 
  exists("predictions_by_dataset"),
  exists("datasets"),
  is.list(predictions_by_dataset),
  is.list(results),
  is.list(datasets)
)
harmon_methods <- c(
  "Raw", "TC", "UQ", "med", "TMM", "TMM.frozen", "DESeq", "PoissonSeq",
  "QN", "QN.frozen", "RUVg", "RUVr", "RUVs"
)
base_methods   <- c("lasso", "pam", "knn", "svm", "svm_with_fs", "ranfor", "ranfor_with_fs")
pretty_base <- c("LASSO", "PAM", "KNN", "SVM", "SVM-FS", "RF", "RF-FS")
pretty_harmon <- c("None", "TC", "UQ", "Med", "TMM", "fTMM",
                   "DESeq", "PoissonSeq", "QN", "fQN", "RUVg", "RUVr", "RUVs")
calculate_metrics <- function(predicted, actual, positive = NULL) {
  stopifnot(length(predicted) == length(actual))
  predicted <- as.character(predicted)
  actual    <- as.character(actual)
  if (is.null(positive)) {
    lvl <- unique(actual)
    positive <- lvl[1]  
  }
  negative <- setdiff(unique(actual), positive)
  if (length(negative) != 1L) {
    acc <- mean(predicted == actual)
    return(c(accuracy = acc, sensitivity = NA_real_, specificity = NA_real_))
  }
  negative <- negative[1]
  ok <- !is.na(predicted) & !is.na(actual)
  pred <- predicted[ok]
  act  <- actual[ok]
  TP <- sum(pred == positive & act == positive)
  TN <- sum(pred == negative & act == negative)
  FP <- sum(pred == positive & act == negative)
  FN <- sum(pred == negative & act == positive)
  denom_acc <- TP + TN + FP + FN
  acc <- if (denom_acc > 0) (TP + TN) / denom_acc else NA_real_
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  
  c(accuracy = acc, sensitivity = sens, specificity = spec)
}
.evaluate_one_dataset <- function(ds_name) {
  if (!ds_name %in% names(datasets)) {
    stop("Dataset '", ds_name, "' not in `datasets`.")
  }
  if (!ds_name %in% names(results)) {
    stop("Dataset '", ds_name, "' not in `results` (for split info).")
  }
  preds <- predictions_by_dataset[[ds_name]]
  split <- results[[ds_name]]$split
  ds    <- datasets[[ds_name]]
  stopifnot(is.list(preds), is.list(split))
  tr <- split$train; te <- split$test
  ytr <- factor(ds$meta$Labels[tr])
  yte <- factor(ds$meta$Labels[te], levels = levels(ytr))
  pos_class <- levels(ytr)[1]
  H <- length(preds[[ base_methods[1] ]]) 
  if (!H %in% c(13L)) {
    message("Note: ", ds_name, " has H = ", H, " harmonizations.")
  }
  empty_mat <- function() {
    matrix(NA_real_, nrow = length(harmon_methods), ncol = length(base_methods),
           dimnames = list(harmon_methods, base_methods))
  }
  mats <- list(
    train = list(accuracy = empty_mat(), sensitivity = empty_mat(), specificity = empty_mat()),
    test1 = list(accuracy = empty_mat(), sensitivity = empty_mat(), specificity = empty_mat()),
    test2 = list(accuracy = empty_mat(), sensitivity = empty_mat(), specificity = empty_mat())
  )
  for (m in base_methods) {
    if (is.null(preds[[m]])) next
    for (h in seq_along(preds[[m]])) {
      pr <- preds[[m]][[h]]
      if (is.null(pr)) next
      tr_pred  <- factor(as.character(pr$train), levels = levels(ytr))
      te1_pred <- factor(as.character(pr$test1), levels = levels(ytr))
      te2_pred <- factor(as.character(pr$test2), levels = levels(ytr))
      tr_met  <- calculate_metrics(tr_pred,  ytr, positive = pos_class)
      te1_met <- calculate_metrics(te1_pred, yte, positive = pos_class)
      te2_met <- calculate_metrics(te2_pred, yte, positive = pos_class)
      mats$train$accuracy   [h, m] <- tr_met["accuracy"]
      mats$train$sensitivity[h, m] <- tr_met["sensitivity"]
      mats$train$specificity[h, m] <- tr_met["specificity"]
      mats$test1$accuracy   [h, m] <- te1_met["accuracy"]
      mats$test1$sensitivity[h, m] <- te1_met["sensitivity"]
      mats$test1$specificity[h, m] <- te1_met["specificity"]
      mats$test2$accuracy   [h, m] <- te2_met["accuracy"]
      mats$test2$sensitivity[h, m] <- te2_met["sensitivity"]
      mats$test2$specificity[h, m] <- te2_met["specificity"]
    }
  }
  out_tbl <- t(round(mats$test2$accuracy, 2))
  rownames(out_tbl) <- pretty_base
  colnames(out_tbl) <- pretty_harmon
  list(
    dataset      = ds_name,
    metrics_mats = mats,                 
    table_test2_accuracy = out_tbl       
  )
}
ds_names <- intersect(names(predictions_by_dataset), names(results))
eval_per_dataset <- lapply(ds_names, .evaluate_one_dataset)
names(eval_per_dataset) <- ds_names
tables_test2_acc <- lapply(eval_per_dataset, `[[`, "table_test2_accuracy")
avg_reps_table <- function(table_list, base_key, reps = 1:10) {
  rep_names <- sprintf("%s__rep%03d", base_key, reps)
  rep_names <- rep_names[rep_names %in% names(table_list)]
  if (length(rep_names) == 0L) stop("No matching replicates found for base_key = ", base_key)
  mats <- lapply(rep_names, function(nm) {
    m <- as.matrix(table_list[[nm]])
    storage.mode(m) <- "numeric"
    m
  })
  rn <- rownames(mats[[1]]); cn <- colnames(mats[[1]])
  mats <- lapply(mats, function(m) m[rn, cn, drop = FALSE])
  
  arr <- array(NA_real_, dim = c(length(rn), length(cn), length(mats)),
               dimnames = list(rn, cn, rep_names))
  for (i in seq_along(mats)) arr[,,i] <- mats[[i]]
  
  avg <- apply(arr, c(1,2), function(x) mean(x, na.rm = TRUE))
  dimnames(avg) <- list(rn, cn)
  attr(avg, "replicates_used") <- rep_names
  avg
}
plot_classif_table <- function(tbl6x13) {
  m <- as.matrix(tbl6x13)
  m[m < 0.5] <- 0.5
  m[m > 1] <- 1
  m <- round(m, 2)
  data.frame(m) %>%
    tibble::rownames_to_column("Classifier") %>%
    gt() %>%
    cols_label(Classifier = "") %>%
    data_color(
      columns = -Classifier,
      colors = scales::col_numeric(
        palette = RColorBrewer::brewer.pal(9, "Blues"),
        domain = c(0.5, 1)
      )
    ) %>%
    cols_align(align = "center") %>%
    bstfun::as_ggplot()
}

avg_balanced <- avg_reps_table(tables_test2_acc, "T_balanced__LumA_vs_LumB", reps = 1:10)
avg_partB    <- avg_reps_table(tables_test2_acc, "T_partB_only__LumA_vs_LumB", reps = 1:10)
avg_confA <- avg_reps_table(tables_test2_acc, "T_confA_L1_A_L2_B__LumA_vs_LumB", reps = 1:10)
avg_confB <- avg_reps_table(tables_test2_acc, "T_confB_L1_B_L2_A__LumA_vs_LumB", reps = 1:10)
p_list <- list(
  plot_classif_table(avg_partB), plot_classif_table(avg_balanced), 
  plot_classif_table(avg_confA),    plot_classif_table(avg_confB)
)
labels <- c("(A)", "(B)",
            "(C)", "(D)")
fig5 <- cowplot::plot_grid(plotlist = p_list, labels = labels, label_size = 10, nrow = 2)









