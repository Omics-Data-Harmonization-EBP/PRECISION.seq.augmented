library(limma)
set.seed(12345)

## 01. biological effects, and handling effects
biological.effects <- function(benchmark, group, c){
  log.benchmark <- log2(benchmark + 1)
  log.benchmark[log.benchmark > 25] <- 25
  group.level <- levels(factor(group))
  sample.g1 <- rowMeans(log.benchmark[,group == group.level[1]])
  sample.g2 <- rowMeans(log.benchmark[,group == group.level[2]])

  amplify_ind = (sample.g2 > sample.g1) + 1
  signal.amplify = abs(sample.g1 - sample.g2)*c

  log.benchmark[, group == group.level[1]] <- log.benchmark[, group == group.level[1]] + (signal.amplify * (amplify_ind == 1))
  log.benchmark[, group == group.level[2]] <- log.benchmark[, group == group.level[2]] + (signal.amplify * (amplify_ind == 2))
  log.benchmark[log.benchmark > 25] <- 25
  log.benchmark[log.benchmark < 0] <- 0

  amp.data <- round(2^log.benchmark - 1)

  return(list(data  = amp.data,
              group = group))
}

handling.effects <- function(clean.input, benchmark, test, group, d){
  log.benchmark <- log2(benchmark + 1)
  log.benchmark[log.benchmark > 25] <- 25
  log.benchmark[log.benchmark < 0] <- 0
  log.test <- log2(test + 1)
  log.test[log.test > 25] <- 25
  log.test[log.test < 0] <- 0

  effects <- log.test - log.benchmark
  log.benchmark.handled <- log2(clean.input + 1)
  log.benchmark.handled <- log.benchmark.handled + effects[,sample(ncol(effects))] * d
  log.benchmark.handled[log.benchmark.handled > 25] <- 25
  log.benchmark.handled[log.benchmark.handled < 0] <- 0

  simulated.benchmark.handled <- round(2^log.benchmark.handled - 1)

  return(list(data   = simulated.benchmark.handled,
              group  = group))
}


## 02. batch effect
# batch.effects <- function(data, groups, correlation_type = "independent", c = 1, partial_index = 0.6){
#   mu1 <- c
#   mu2 <- 1.8 * c
#   sigma <- 1
#   n <- ncol(data)
#   m <- nrow(data)
#
#   batch_assignments <- vector("character", n)
#   if (correlation_type == "independent") {
#     batch_assignments <- sample(c("Batch1", "Batch2"), n, replace = TRUE)
#   } else if (correlation_type == "partial") {
#     unique_groups <- unique(groups)
#     for (group in unique_groups) {
#       indices <- which(groups == group)
#       size <- length(indices)
#       num_batch1 <- round(size * partial_index)
#       num_batch2 <- size - num_batch1
#       if (group == unique_groups[1]) {
#         batch_assignments[indices] <- c(rep("Batch1", num_batch1), rep("Batch2", num_batch2))
#       } else if (group == unique_groups[2]) {
#         batch_assignments[indices] <- c(rep("Batch2", num_batch1), rep("Batch1", num_batch2))
#       }
#       batch_assignments[indices] <- sample(batch_assignments[indices])
#     }
#   } else if (correlation_type == "aligned") {
#     unique_groups <- unique(groups)
#     batch_assignments <- groups
#     batch_assignments[batch_assignments == unique_groups[1]] <- "Batch1"
#     batch_assignments[batch_assignments == unique_groups[2]] <- "Batch2"
#   }
#
#   batch_effects <- matrix(nrow = m, ncol = n)
#   for (j in 1:n) {
#     batch_effects[, j] <- ifelse(batch_assignments[j] == "Batch1",
#                                  rnorm(m, mean = mu1, sd = sigma),
#                                  rnorm(m, mean = mu2, sd = sigma))
#   }
#
#   data[data < 0] <- 0
#   log_data <- log2(data + 1)
#   adjusted_data <- log_data + batch_effects
#   adjusted_data[adjusted_data > 25] <- 25
#   adjusted_data[adjusted_data < 0] <- 0
#   batch_added_data <- round(2^adjusted_data - 1)
#
#   return(list(data = batch_added_data,
#               group = groups,
#               batch = batch_assignments))
# }
#

batch.effects <- function(data, groups, c = 1, confound_level = 0.6){

  data <- as.matrix(data)
  n <- ncol(data)
  m <- nrow(data)

  batch_means <- c(-1.5 * c, 1.5 * c, -1 * c, 1 * c,  0)
  batch_sds <- rep(1, 5)
  batch_assignments <- rep(0, n)

  group1_indices <- which(groups == 'MXF')
  group2_indices <- which(groups == 'PMFH')

  batch_size <- floor(n / 5)
  num_group1_batch12 <- 2 * round(batch_size * confound_level)
  num_group2_batch12 <- 2 * batch_size - num_group1_batch12
  num_group2_batch34 <- 2 * round(batch_size * confound_level)
  num_group1_batch34 <- 2 * batch_size - num_group2_batch34
  group1_batch12 <- sample(group1_indices, num_group1_batch12)
  group1_batch34 <- sample(setdiff(group1_indices, group1_batch12), num_group1_batch34)
  group2_batch12 <- sample(group2_indices, num_group2_batch12)
  group2_batch34 <- sample(setdiff(group2_indices, group2_batch12), num_group2_batch34)

  batch_assignments[group1_batch12] <- if(num_group1_batch12 %% 2 == 0) {sample(rep(c(1, 2), each = num_group1_batch12 / 2))} else {c(sample(rep(c(1, 2), each = num_group1_batch12 %/% 2)), sample(c(1,2), 1))}
  batch_assignments[group2_batch12] <- if(num_group2_batch12 %% 2 == 0) {sample(rep(c(1, 2), each = num_group2_batch12 / 2))} else {c(sample(rep(c(1, 2), each = num_group2_batch12 %/% 2)), sample(c(1,2), 1))}
  batch_assignments[group1_batch34] <- if(num_group1_batch34 %% 2 == 0) {sample(rep(c(3, 4), each = num_group1_batch34 / 2))} else {c(sample(rep(c(3, 4), each = num_group1_batch34 %/% 2)), sample(c(3,4), 1))}
  batch_assignments[group2_batch34] <- if(num_group2_batch34 %% 2 == 0) {sample(rep(c(3, 4), each = num_group2_batch34 / 2))} else {c(sample(rep(c(3, 4), each = num_group2_batch34 %/% 2)), sample(c(3,4), 1))}
  unassigned_indices <- which(batch_assignments == 0)
  batch_assignments[unassigned_indices] <- 5

  batch_effects <- matrix(nrow = m, ncol = n)
  for (j in 1:n) {
    batch_effects[, j] <- rnorm(m, mean = batch_means[batch_assignments[j]], sd = batch_sds[batch_assignments[j]])
  }

  log_data <- log2(data + 1)
  log_data[log_data < 0] <- 0
  adjusted_data <- log_data + batch_effects
  adjusted_data[adjusted_data > 25] <- 25
  adjusted_data[adjusted_data < 0] <- 0
  batch_added_data <- round(2^adjusted_data - 1)

  return(list(data = batch_added_data,
              group = groups,
              batch = batch_assignments))
}








