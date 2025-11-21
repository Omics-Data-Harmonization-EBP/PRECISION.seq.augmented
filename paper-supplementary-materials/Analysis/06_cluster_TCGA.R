rm(list = ls())
library(parallel)
library(sva)
library(EDASeq)
library(edgeR)
library(RUVSeq)
library(Biobase)
library(BiocGenerics)
library(tidyverse)
library(mclust)
library(aricode)
library(RSKC)
library(cluster)
library(factoextra)
library(caret)
library(e1071)
library(sparseSVM)
library(caret)
library(TCGAbiolinks)
library(PRECISION.seq.augmented)
library(gt)
library(abind)
library(gridExtra)
load("TCGA_BRCA_miRNA_matched_batch_included.RData")
flag_dates <- as.Date(c("2011-05-17", "2011-05-31", "2011-06-14", "2011-07-12", "2011-08-31"))
sequence_depth <- data.frame(total_count = colSums(brca_expr),
                             date = brca_meta$ShipDate) %>%
  mutate(date        = as.Date(date),               
         era         = if_else(date < as.Date("2011-08-01"),
                               "pre‑2011‑08",   
                               "post‑2011‑08"),  
         date_factor = factor(date, levels = sort(unique(date))),
         is_flagged  = date %in% flag_dates,
         fill_group  = if_else(is_flagged, "flagged", era))
era_cols <- c("pre‑2011‑08"  = "#A8D2F0",   
              "post‑2011‑08" = "#F4B6B3")
p <- ggplot(sequence_depth,
            aes(x = date_factor,
                y = log2(total_count),
                fill = fill_group)) +         
  geom_boxplot(outlier.shape = NA,
               width = 0.55,
               size  = 0.4,
               colour = "black") +
  scale_fill_manual(values = c(flagged = "#cecece", era_cols), guide = "none") +  
  labs(x = "Sequencing Date",
       y = "Read Depth") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_blank(),  
    axis.ticks.x       = element_blank()
  )
ship_dates_as_date <- as.Date(brca_meta$ShipDate)
exclude_idx <- which(ship_dates_as_date %in% flag_dates)
if (length(exclude_idx) > 0) {
  brca_meta <- brca_meta[-exclude_idx, , drop = FALSE]
  brca_expr <- brca_expr[, -exclude_idx, drop = FALSE]
}
rownames(brca_meta) <- brca_meta$Sample
brca_meta <- as.data.frame(brca_meta, stringsAsFactors = FALSE)
stopifnot(ncol(brca_expr) == nrow(brca_meta))
if (!("Subtype" %in% names(brca_meta))) {
  label_map <- c(BRCA_1="LumA", BRCA_2="LumB", BRCA_3="Basal", BRCA_4="HER2")
  brca_meta$Subtype <- label_map[brca_meta$Labels]
}
if (!("era" %in% names(brca_meta))) stop("brca_meta must have an 'era' column ('Early'/'Late').")
brca_meta$Subtype <- factor(brca_meta$Subtype, levels = c("LumA","LumB","Basal","HER2"))
brca_meta$era     <- factor(brca_meta$era,     levels = c("Early","Late"))
PAIRS <- list(
  LumA_vs_Basal = c("LumA","Basal")
)
N_REPS          <- 10      
SEED_BASE       <- 42
REPLACE_WITHIN  <- FALSE  
STRICT_COUNTS   <- FALSE 
.pool_idx <- function(meta, subtype, era) {
  which(meta$Subtype == subtype & meta$era == era)
}
.take_or_scale_layout <- function(meta, pair, layout, strict = STRICT_COUNTS) {
  stopifnot(all(rownames(layout) == c("Early","Late")))
  stopifnot(all(colnames(layout) == pair))
  cap <- matrix(0L, nrow=2, ncol=2, dimnames=dimnames(layout))
  for (er in rownames(layout)) for (lab in colnames(layout)) {
    cap[er, lab] <- length(.pool_idx(meta, lab, er))
  }
  need <- layout
  ok <- all(need <= cap)
  if (ok) return(list(layout=need, scaled=FALSE, scale_factor=1.0, cap=cap))
  if (strict) {
    msg <- paste0("Insufficient capacity for layout\nNeeded:\n", capture.output(print(need)), 
                  "\nCaps:\n", capture.output(print(cap)), collapse="\n")
    stop(msg)
  }
  ratio_mat <- ifelse(need > 0, cap / pmax(1, need), Inf)
  s <- min(ratio_mat[is.finite(ratio_mat) & ratio_mat > -Inf], 1.0)
  s <- max(0, s)  # guard
  adj <- floor(need * s)
  list(layout=adj, scaled=TRUE, scale_factor=s, cap=cap)
}
.sample_from_layout <- function(meta, expr, pair, layout, replace = REPLACE_WITHIN) {
  idx <- integer(0)
  for (er in rownames(layout)) for (lab in colnames(layout)) {
    n <- as.integer(layout[er, lab])
    if (n > 0) {
      pool <- .pool_idx(meta, lab, er)
      idx  <- c(idx, sample(pool, n, replace = replace))
    }
  }
  list(
    expr = expr[, idx, drop = FALSE],
    meta = meta[idx, , drop = FALSE],
    idx  = idx
  )
}
.make_note <- function(scen, pair, layout, scaled, sf) {
  paste0(
    sprintf("%s (%s vs %s) — target cell counts [Early/Late x %s/%s] = ",
            scen, pair[1], pair[2], pair[1], pair[2]),
    paste(c(layout["Early",pair[1]], layout["Early",pair[2]],
            layout["Late",pair[1]],  layout["Late",pair[2]]), collapse="/"),
    if (scaled) sprintf("; NOTE: scaled by %.3f to meet capacity.", sf) else ""
  )
}
layout_s1 <- function(pair) { 
  m <- matrix(c(120, 120,
                50,  50),
              nrow=2, byrow=TRUE,
              dimnames=list(c("Early","Late"), pair))
  m
}
layout_s2_allLate <- function(pair) {
  matrix(c(  0,   0,
             50,  50),
         nrow=2, byrow=TRUE,
         dimnames=list(c("Early","Late"), pair))
}
layout_s2_balanced <- function(pair) {  # 50:50 era for same totals
  matrix(c( 25, 25,
            25, 25),
         nrow=2, byrow=TRUE,
         dimnames=list(c("Early","Late"), pair))
}
layout_s3_confA <- function(pair) {
  matrix(c(50,  0,
           0,  50),
         nrow=2, byrow=TRUE,
         dimnames=list(c("Early","Late"), pair))
}
layout_s3_confB <- function(pair) { 
  matrix(c(  0, 50,
             50,  0),
         nrow=2, byrow=TRUE,
         dimnames=list(c("Early","Late"), pair))
}
SCENARIOS <- list(
  S1_signal_strength      = layout_s1,
  S2_allLate              = layout_s2_allLate,
  S2_eraBalanced          = layout_s2_balanced,
  S3_confA_L1Early_L2Late = layout_s3_confA,
  S3_confB_L1Late_L2Early = layout_s3_confB
)

datasets <- list()
ledger_rows <- list()
rep_id_fmt <- function(k) sprintf("rep%03d", k)

for (pair_name in names(PAIRS)) {
  pair <- PAIRS[[pair_name]]
  
  for (scen_name in names(SCENARIOS)) {
    layout_fn <- SCENARIOS[[scen_name]]
    
    for (k in seq_len(N_REPS)) {
      set.seed(SEED_BASE + k + 1000 * match(scen_name, names(SCENARIOS)) + 100000 * match(pair_name, names(PAIRS)))
      
      target <- layout_fn(pair)
      feas   <- .take_or_scale_layout(brca_meta, pair, target, strict = STRICT_COUNTS)
      draw   <- .sample_from_layout(brca_meta, brca_expr, pair, feas$layout, replace = REPLACE_WITHIN)
      
      ds_name <- paste(scen_name, pair_name, rep_id_fmt(k), sep = "__")
      note    <- .make_note(scen_name, pair, feas$layout, feas$scaled, feas$scale_factor)
      
      datasets[[ds_name]] <- list(expr = draw$expr, meta = draw$meta, note = note)
      
      # ledger row
      tabL <- table(factor(draw$meta$Subtype, levels=pair))
      tabE <- table(factor(draw$meta$era,     levels=c("Early","Late")))
      cell_counts <- c(
        Early_L1  = sum(draw$meta$Subtype==pair[1] & draw$meta$era=="Early"),
        Early_L2  = sum(draw$meta$Subtype==pair[2] & draw$meta$era=="Early"),
        Late_L1   = sum(draw$meta$Subtype==pair[1] & draw$meta$era=="Late"),
        Late_L2   = sum(draw$meta$Subtype==pair[2] & draw$meta$era=="Late")
      )
      ledger_rows[[length(ledger_rows)+1L]] <- data.frame(
        dataset = ds_name,
        scenario = scen_name,
        pair = paste(pair, collapse=" vs "),
        replicate = k,
        N = ncol(draw$expr),
        note = note,
        LumA = as.integer(tabL["LumA"]), LumB = as.integer(tabL["LumB"]),
        Basal = as.integer(tabL["Basal"]), HER2 = as.integer(tabL["HER2"]),
        Early = as.integer(tabE["Early"]), Late = as.integer(tabE["Late"]),
        Early_L1 = cell_counts["Early_L1"],
        Early_L2 = cell_counts["Early_L2"],
        Late_L1  = cell_counts["Late_L1"],
        Late_L2  = cell_counts["Late_L2"],
        scaled = feas$scaled,
        scale_factor = feas$scale_factor,
        stringsAsFactors = FALSE, check.names = FALSE
      )
    }
  }
}
scenario_ledger <- dplyr::bind_rows(ledger_rows)
save(datasets, file = "brca_subsampled_simulation_with_reps_v4.RData")


## =========================================================
## Run all subsampled analyses (replicate-aware, mclapply) & save results
## =========================================================
set.seed(42)
library(parallel)
library(edgeR)
library(tidyverse)
library(EDASeq)
library(RUVSeq)
library(mclust)
library(cluster)
library(sva)
library(PRECISION.seq.augmented)
source("RUV_without_effects.R")
load("brca_subsampled_simulation_with_reps_v4.RData") 
MC_CORES <- as.integer(Sys.getenv("MC_CORES",
                                  unset = max(1L, parallel::detectCores(logical = TRUE) - 1L)))
message("Using mclapply with ", MC_CORES, " cores")
harmon.all <- function(object){
  object@harmon.train.data$Raw$dat.harmonized <- object@raw.train.data$data
  object <- harmon.TC(object)
  object <- harmon.UQ(object)
  object <- harmon.med(object)
  object <- harmon.TMM(object)
  object <- harmon.DESeq(object)
  object <- harmon.PoissonSeq(object)
  object <- harmon.QN(object)
  object <- harmon.RUVr(object)
  object <- harmon.RUVs(object)
  object <- harmon.RUVg(object)
  return(object)
}
cluster.all <- function(object){
  object <- cluster.hc(object, distance = "euclidean")
  object <- cluster.hc(object, distance = "pearson")
  object <- cluster.hc(object, distance = "spearman")
  object <- cluster.mnm(object)
  object <- cluster.kmeans(object)
  object <- cluster.som(object)
  object <- cluster.pam(object, distance = "euclidean")
  object <- cluster.pam(object, distance = "pearson")
  object <- cluster.pam(object, distance = "spearman")
  return(object)
}
measure <- function(data_expr, data_label, data_batch){
  data_expr <- as.matrix(data_expr)
  keep <-  rowSums(data_expr > 0) >= ceiling(0.01 * ncol(data_expr))
  data_expr <- data_expr[keep, , drop = FALSE]
  analysis <- create.precision.cluster(data = data_expr, label = data_label)
  analysis <- harmon.all(analysis)
  analysis <- cluster.all(analysis)
  cluster_methods <- c('hc_euclidean','hc_pearson','hc_spearman',
                       'mnm','kmeans','som',
                       'pam_euclidean','pam_pearson','pam_spearman')
  harmon_methods  <- c('Raw','TC','UQ','med','TMM','DESeq','PoissonSeq','QN','RUVr','RUVs','RUVg')
  ari_indexes <- matrix(NA_real_, nrow=length(cluster_methods), ncol=length(harmon_methods))
  silhouette_indexes <- matrix(NA_real_, nrow=length(cluster_methods), ncol=length(harmon_methods))
  true_label <- as.factor(data_label)
  for (i in seq_along(cluster_methods)){
    clm <- cluster_methods[i]
    for (j in seq_along(harmon_methods)){
      hm <- harmon_methods[j]
      est_cluster <- try({
        if (startsWith(clm, "hc_")) {
          distance <- sub("hc_", "", clm)
          analysis@cluster.result$hc[[distance]][[hm]]
        } else if (startsWith(clm, "pam_")) {
          distance <- sub("pam_", "", clm)
          analysis@cluster.result$pam[[distance]][[hm]]
        } else {
          analysis@cluster.result[[clm]][[hm]]
        }
      }, silent = TRUE)
      if (inherits(est_cluster, "try-error") || is.null(est_cluster)) {
        ari_indexes[i,j] <- NA_real_
        silhouette_indexes[i,j] <- NA_real_
        next
      }
      est_cluster <- as.factor(est_cluster)
      ari_indexes[i,j] <- tryCatch(
        mclust::adjustedRandIndex(true_label, est_cluster),
        error = function(e) NA_real_
      )
      sil_val <- NA_real_
      if (length(unique(est_cluster)) >= 2L && length(unique(est_cluster)) <= (ncol(data_expr) - 1L)) {
        dist_matrix <- try({
          Xh <- analysis@harmon.train.data[[hm]]$dat.harmonized
          Xh <- as.matrix(Xh)
          if (startsWith(clm, "hc_") || startsWith(clm, "pam_")) {
            switch(sub(".*_", "", clm),
                   "euclidean" = dist(t(Xh)),
                   "pearson"   = as.dist(1 - cor(Xh, method = "pearson")),
                   "spearman"  = as.dist(1 - cor(Xh, method = "spearman")))
          } else {
            dist(t(as.matrix(Xh)))
          }
        }, silent = TRUE)
        
        if (!inherits(dist_matrix, "try-error")) {
          sil <- try(suppressWarnings(cluster::silhouette(as.integer(est_cluster), dist_matrix)), silent = TRUE)
          if (!inherits(sil, "try-error")) {
            sil_val <- mean(sil[, "sil_width"])
          }
        }
      }
      silhouette_indexes[i,j] <- sil_val
    }
  }
  rownames(ari_indexes) <- cluster_methods
  colnames(ari_indexes) <- harmon_methods
  rownames(silhouette_indexes) <- cluster_methods
  colnames(silhouette_indexes) <- harmon_methods
  list(ari = as.data.frame(ari_indexes),
       silhouette = as.data.frame(silhouette_indexes))
}
parse_name <- function(nm) {
  if (grepl("__rep\\d+$", nm)) {
    parts <- strsplit(nm, "__", fixed = TRUE)[[1]]
    scen_full <- parts[1]
    pair_key  <- parts[2]
    rep_str   <- parts[3]
    scenario_id <- sub("^(S\\d+).*", "\\1", scen_full)
    variant_raw <- tolower(sub("^S\\d+_", "", scen_full))
    variant <- dplyr::case_when(
      variant_raw %in% c("alllate", "all_late", "all-late") ~ "All-Late",
      variant_raw %in% c("erabalanced","era_balanced","era-balanced") ~ "Era-balanced",
      variant_raw %in% c("duplicatebotheras","dupboth","duplicateboth") ~ "DupBothEras",
      grepl("^signal", variant_raw) ~ "Signal-strength",
      grepl("^confa", variant_raw) ~ "confA",
      grepl("^confb", variant_raw) ~ "confB",
      TRUE ~ variant_raw
    )
    pair <- gsub("_vs_", " vs ", pair_key)
    replicate <- as.integer(sub("^rep", "", rep_str))
    tibble(
      scenario_full = scen_full,
      scenario = scenario_id,
      pair = pair,
      variant = variant,
      replicate = replicate,
      dataset = nm
    )
  } else {
    scen <- sub("^s(\\d+).*", "S\\1", nm)
    pair <- if (grepl("LumA_vs_LumB", nm)) "LumA vs LumB"
    else if (grepl("LumA_vs_Basal", nm)) "LumA vs Basal"
    else if (grepl("fourLabels", nm, ignore.case = TRUE)) "4-labels"
    else "other"
    variant <- if (grepl("allLate", nm, ignore.case = TRUE)) "All-Late"
    else if (grepl("eraBalanced", nm, ignore.case = TRUE)) "Era-balanced"
    else if (grepl("confA", nm, ignore.case = TRUE)) "confA"
    else if (grepl("confB", nm, ignore.case = TRUE)) "confB"
    else if (grepl("balanced", nm, ignore.case = TRUE) && !grepl("eraBalanced", nm, ignore.case = TRUE)) "Balanced"
    else if (grepl("imbalanced", nm, ignore.case = TRUE)) "Imbalanced"
    else "—"
    tibble(scenario_full = scen, scenario = scen, pair = pair, variant = variant,
           replicate = 1L, dataset = nm)
  }
}
cluster_levels <- c('hc_euclidean','hc_pearson','hc_spearman','mnm','kmeans','som',
                    'pam_euclidean','pam_pearson','pam_spearman')
harmon_levels  <- c('Raw','TC','UQ','med','TMM','DESeq','PoissonSeq','QN','RUVr','RUVs','RUVg')
nm_vec <- names(datasets)
message("Running analyses across ", length(nm_vec), " datasets with mclapply...")
run_one <- function(nm){
  d <- datasets[[nm]]
  stopifnot(is.matrix(d$expr) || is.data.frame(d$expr))
  y_label <- d$meta$Labels
  y_batch <- d$meta$era
  res <- measure(as.matrix(d$expr), y_label, y_batch)
  list(name = nm, note = d$note, res = res)
}
placeholder_result <- function() {
  empty_mat <- matrix(NA_real_,
                      nrow = length(cluster_levels),
                      ncol = length(harmon_levels),
                      dimnames = list(cluster_levels, harmon_levels))
  list(
    ari = as.data.frame(empty_mat),
    silhouette = as.data.frame(empty_mat)
  )
}
result_list <- parallel::mclapply(
  nm_vec,
  function(nm) {
    out <- try(run_one(nm), silent = TRUE)
    if (inherits(out, "try-error")) {
      warning("Failed on dataset: ", nm, "\n  Reason: ", as.character(out))
      list(name = nm, note = datasets[[nm]]$note, res = placeholder_result())
    } else {
      out
    }
  },
  mc.cores = MC_CORES,
  mc.set.seed = TRUE,
  mc.preschedule = TRUE
)
names(result_list) <- nm_vec
tidy_one <- function(name, res){
  meta_row <- parse_name(name)
  ari <- res$ari %>% rownames_to_column("cluster") %>%
    pivot_longer(-cluster, names_to = "harmon", values_to = "ARI")
  sil <- res$silhouette %>% rownames_to_column("cluster") %>%
    pivot_longer(-cluster, names_to = "harmon", values_to = "Silhouette")
  full_join(ari, sil, by = c("cluster","harmon")) %>%
    mutate(dataset = name) %>%
    left_join(meta_row, by = "dataset")
}
results_tidy <- bind_rows(lapply(result_list, function(x) tidy_one(x$name, x$res))) %>%
  mutate(
    cluster  = factor(cluster, levels = cluster_levels),
    harmon   = factor(harmon,  levels = harmon_levels),
    scenario = factor(scenario, levels = sort(unique(scenario)))
  )
results_summary <- results_tidy %>%
  group_by(scenario, pair, variant, cluster, harmon) %>%
  summarise(
    n_reps = n_distinct(replicate),
    ARI_mean = mean(ARI, na.rm = TRUE),
    ARI_sd   = sd(ARI, na.rm = TRUE),
    ARI_se   = ARI_sd/sqrt(n_reps),
    ARI_q25  = quantile(ARI, 0.25, na.rm = TRUE),
    ARI_q50  = quantile(ARI, 0.50, na.rm = TRUE),
    ARI_q75  = quantile(ARI, 0.75, na.rm = TRUE),
    Sil_mean = mean(Silhouette, na.rm = TRUE),
    Sil_sd   = sd(Silhouette, na.rm = TRUE),
    Sil_se   = Sil_sd/sqrt(n_reps),
    Sil_q25  = quantile(Silhouette, 0.25, na.rm = TRUE),
    Sil_q50  = quantile(Silhouette, 0.50, na.rm = TRUE),
    Sil_q75  = quantile(Silhouette, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    cluster  = factor(cluster, levels = cluster_levels),
    harmon   = factor(harmon,  levels = harmon_levels),
    scenario = factor(scenario, levels = sort(unique(scenario)))
  )
save(result_list, results_tidy, results_summary,
     file = "brca_subsampled_results_with_reps_v4.RData")