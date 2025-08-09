library(PRECISION.seq.augmented)
augmented.data <- load.augmented.data(temp = TRUE)
benchmark <- augmented.data$benchmark
test <- augmented.data$test
rm(augmented.data)

# load("./data/data_source/MSKpair_60000_count.RData")

# For reproducibility:
set.seed(42)

# --- Step 1: Randomly select which 200 elements to keep ---
selected_indices <- sample(seq_len(400), size = 300)

# Initialize new lists to store the subsampled data
benchmark_sub <- vector("list", length = 300)
test_sub <- vector("list", length = 300)

# --- Step 2: For each selected element, subsample 100 MXF and 100 PMFH ---
for (i in seq_along(selected_indices)) {
  idx <- selected_indices[i]

  # Extract the original benchmark element
  bench_elem <- benchmark[[idx]]
  # bench_elem$data: genes in rows, 300 samples in columns
  # bench_elem$label: length-300 vector with 150 "MXF" and 150 "PMFH"

  # Identify the columns for each label
  mxf_indices <- which(bench_elem$label == "MXF")
  pmfh_indices <- which(bench_elem$label == "PMFH")

  # Randomly select 100 columns from each label
  sub_mxf <- sample(mxf_indices, 100)
  sub_pmfh <- sample(pmfh_indices, 100)

  # Combine these column indices
  sub_all <- c(sub_mxf, sub_pmfh)

  # Subset benchmark data by columns (since samples are in columns)
  new_bench_data <- bench_elem$data[, sub_all, drop = FALSE]
  new_bench_label <- bench_elem$label[sub_all]

  # Now handle the paired test element
  test_elem <- test[[idx]]
  # Subset test data with the same column indices
  new_test_data <- test_elem$data[, sub_all, drop = FALSE]
  new_test_label <- test_elem$label[sub_all]

  # Store the subsets in new lists
  benchmark_sub[[i]] <- list(data = new_bench_data, label = new_bench_label)
  test_sub[[i]] <- list(data = new_test_data, label = new_test_label)
}

# benchmark_sub and test_sub now each have 300 elements.
# In each element, data has 200 columns (100 MXF + 100 PMFH).
# The pairing is preserved both at the element level (same index)
# and at the sample level (same columns in each paired element).

save(benchmark_sub, test_sub, file = "./data/data_source/MSKpair_300_cluster.RData")
