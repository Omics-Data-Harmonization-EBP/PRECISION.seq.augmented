#' Add biological effects to benchmark data
#'
#' This function simulates biological effects on benchmark data by amplifying
#' the signal between \emph{two} groups based on a specified amplification factor.
#' It modifies the input data and returns the amplified data along with group information.
#'
#' @param benchmark Numeric matrix of expression values, (\code{p} by \code{n}), where
#' \code{p} is the number of features and \code{n} is the number of samples.
#' @param group Factor or vector of character strings defining sample groups
#'  for each sample. Must have exactly two levels.
#' @param c Numeric amplification factor
#' @return List containing modified data, group information, and amplification factor.
#' @export
biological.effects <- function(benchmark, group, c) {
  log.benchmark <- log2(benchmark + 1)
  log.benchmark[log.benchmark > 25] <- 25
  group.level <- levels(factor(group))
  sample.g1 <- rowMeans(log.benchmark[, group == group.level[1]])
  sample.g2 <- rowMeans(log.benchmark[, group == group.level[2]])

  amplify_ind <- (sample.g2 > sample.g1) + 1
  signal.amplify <- abs(sample.g1 - sample.g2) * c

  log.benchmark[, group == group.level[1]] <-
    log.benchmark[, group == group.level[1]] + (signal.amplify * (amplify_ind == 1))
  log.benchmark[, group == group.level[2]] <-
    log.benchmark[, group == group.level[2]] + (signal.amplify * (amplify_ind == 2))
  log.benchmark[log.benchmark > 25] <- 25
  log.benchmark[log.benchmark < 0] <- 0

  amp.data <- round(2^log.benchmark - 1)

  return(list(
    data = amp.data,
    group = group,
    amplification_factor = c
  ))
}

#' Add handling effects to benchmark data
#'
#' This function simulates handling effects by modifying the
#' clean input data based on the differences between the benchmark and test data.
#'
#' @param clean.input Numeric matrix of clean input data, (\code{p} by \code{n}), where
#' \code{p} is the number of features and \code{n} is the number of samples.
#' @param benchmark Numeric matrix of benchmark data, (\code{p} by \code{n}), where
#' \code{p} is the number of features and \code{n} is the number of samples.
#' @param test Numeric matrix of test data, (\code{p} by \code{n}), where
#' \code{p} is the number of features and \code{n} is the number of samples.
#' @param group Factor or vector of character strings defining sample groups
#'  for each sample. Must have exactly two levels.
#' @param d Numeric effect strength factor
#' @return List containing modified data, group information, and effect size
#' @export
handling.effects <- function(clean.input, benchmark, test, group, d) {
  log.benchmark <- log2(benchmark + 1)
  log.benchmark[log.benchmark > 25] <- 25
  log.benchmark[log.benchmark < 0] <- 0
  log.test <- log2(test + 1)
  log.test[log.test > 25] <- 25
  log.test[log.test < 0] <- 0

  effects <- log.test - log.benchmark
  log.benchmark.handled <- log2(clean.input + 1)
  log.benchmark.handled <- log.benchmark.handled + effects[, sample(ncol(effects))] * d
  log.benchmark.handled[log.benchmark.handled > 25] <- 25
  log.benchmark.handled[log.benchmark.handled < 0] <- 0

  simulated.benchmark.handled <- round(2^log.benchmark.handled - 1)

  return(list(
    data = simulated.benchmark.handled,
    group = group,
    effect_size = d
  ))
}
