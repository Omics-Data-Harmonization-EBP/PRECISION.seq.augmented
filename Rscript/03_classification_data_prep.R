load("./data/data_source/MSKpair_60000_count.RData")

#-----------------------------
# 1) Merge all "benchmark" elements
#-----------------------------
big_bench_data_list  <- lapply(benchmark, function(x) x$data)
big_bench_data       <- do.call(cbind, big_bench_data_list)
# big_bench_data is now (genes x 120,000)

big_bench_label_list <- lapply(benchmark, function(x) x$label)
big_bench_label      <- unlist(big_bench_label_list)
# big_bench_label is length 120,000, with 60,000 "MXF" and 60,000 "PMFH"

#-----------------------------
# 2) Merge all "test" elements in the same order
#-----------------------------
big_test_data_list  <- lapply(test, function(x) x$data)
big_test_data       <- do.call(cbind, big_test_data_list)
# (genes x 120,000)

big_test_label_list <- lapply(test, function(x) x$label)
big_test_label      <- unlist(big_test_label_list)
# length 120,000, also 60,000 "MXF" + 60,000 "PMFH"

# IMPORTANT ASSUMPTION: The j-th column in big_bench_data
# aligns with the j-th column in big_test_data (same sample pairing).

#-----------------------------
# 3) Identify all MXF and PMFH columns
#-----------------------------
mxf_indices  <- which(big_bench_label == "MXF")   # length 60,000
pmfh_indices <- which(big_bench_label == "PMFH")  # length 60,000

# We can form up to 300 disjoint subsets of size 200 from each group:
n_subset <- min(length(mxf_indices), length(pmfh_indices)) / 200  # 300

#-----------------------------
# 4) Split each label group into 300 disjoint chunks
#-----------------------------
# Sort them in ascending order so we can chunk consistently
mxf_indices_sorted  <- sort(mxf_indices)
pmfh_indices_sorted <- sort(pmfh_indices)

# Now chunk them: chunk 1 will get the first 200, chunk 2 the next 200, etc.
mxf_chunks  <- split(mxf_indices_sorted,  rep(1:n_subset, each = 200))
pmfh_chunks <- split(pmfh_indices_sorted, rep(1:n_subset, each = 200))

#-----------------------------
# 5) For each chunk, build a 400-column subset for benchmark & test
#-----------------------------
benchmark_subsets <- vector("list", length = n_subset)
test_subsets      <- vector("list", length = n_subset)

for (i in seq_len(n_subset)) {
  # Gather the columns for the i-th subset
  sub_mxf  <- mxf_chunks[[i]]   # 200 MXF columns
  sub_pmfh <- pmfh_chunks[[i]]  # 200 PMFH columns
  sub_cols <- c(sub_mxf, sub_pmfh)  # total 400 columns
  
  # Subset from the big benchmark data
  tmp_bench_data  <- big_bench_data[, sub_cols, drop = FALSE]
  tmp_bench_label <- big_bench_label[sub_cols]
  
  # Subset from the big test data
  tmp_test_data  <- big_test_data[, sub_cols, drop = FALSE]
  tmp_test_label <- big_test_label[sub_cols]
  
  # Store as a list (genes x 400 for data, length-400 for label)
  benchmark_subsets[[i]] <- list(data  = tmp_bench_data,
                                 label = tmp_bench_label)
  test_subsets[[i]]      <- list(data  = tmp_test_data,
                                 label = tmp_test_label)
}

save(benchmark_subsets, test_subsets, file = "./data/data_source/MSKpair_300_classification.RData")
