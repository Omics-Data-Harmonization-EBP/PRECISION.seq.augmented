library(limma)
set.seed(12345)

## 01. biological effects, and handling effects
biological.effects <- function(benchmark, group, c){
  log.benchmark <- log2(benchmark + 1)
  group.level <- levels(factor(group))
  sample.g1 <- rowMeans(log.benchmark[,group == group.level[1]])
  sample.g2 <- rowMeans(log.benchmark[,group == group.level[2]])

  amplify_ind = (sample.g2 > sample.g1) + 1
  signal.amplify = abs(sample.g1 - sample.g2)*c

  log.benchmark[, group == group.level[1]] <- log.benchmark[, group == group.level[1]] + (signal.amplify * (amplify_ind == 1))
  log.benchmark[, group == group.level[2]] <- log.benchmark[, group == group.level[2]] + (signal.amplify * (amplify_ind == 2))
  log.benchmark[log.benchmark > log2(25)] <- log2(25)

  amp.data <- round(2^log.benchmark - 1)
  amp.data[amp.data < 0] <- 0

  return(list(data  = amp.data,
              group = group))
}

handling.effects <- function(clean.input, benchmark, test, group, d){
  effects <- log2(test + 1) - log2(benchmark + 1)
  log.benchmark.handled <- log2(clean.input + 1)
  log.benchmark.handled <- log.benchmark.handled + effects * d
  log.benchmark.handled[log.benchmark.handled > log2(25)] <- log2(25)

  simulated.benchmark.handled <- round(2^log.benchmark.handled - 1)
  simulated.benchmark.handled[simulated.benchmark.handled < 0] <- 0

  return(list(data   = simulated.benchmark.handled,
              group  = group))
}




