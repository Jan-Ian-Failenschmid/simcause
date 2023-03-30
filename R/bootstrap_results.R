#' Extract Results from FCI Bootstrap Object
#'
#' Extracts proportion of time a graph feature has been identified by a causal inference algorithm
#' within the bootstrap samples created by fci_bootsrap.
#'
#' @param bootstrap_samples bootstrap_samples object created by fci_bootstrap.
#' @param feature Feature for which results should be extracted. Options are: "Edges", "Arrowheads", "Tails", and "Circles".
#' @param CI Logical whether CI for proportions should be calculated.
#' @param alpha Alpha level that should be used for the calculated confidence intervals.
#'
#' @returns A data frame holding the proportion of times each featre has been identified and a corresponding confidence interval.
#'
#' @examples
#' # Edges_results <- bootstrap_results(
#' #   bootstrap_samples = bootstrap_samples,
#' #   feature = "Edges"
#' # )
#' # Arrowheads_results <- bootstrap_results(
#' #   bootstrap_samples = bootstrap_samples,
#' #   feature = "Arrowheads"
#' # )
#' # Tails_results <- bootstrap_results(
#' #   bootstrap_samples = bootstrap_samples,
#' #   feature = "Tails"
#' # )
#' # Circles_results <- bootstrap_results(
#' #   bootstrap_samples = bootstrap_samples,
#' #   feature = "Circles"
#' # )
#'
#' @export

bootstrap_results <- function(bootstrap_samples, feature, CI = TRUE, alpha = 0.05) {

  # Extract label columns from the bootstrap_samples object
  result <- bootstrap_samples[, feature][[1]][, c(1, 2)]

  # Apply summary_stat mean to get Proportion of identified feature
  result$Proportion <- summary_stat(bootstrap_samples[, feature], mean)

  # Apply summary_stat sd to get standard deviation of identified feature
  result$SD <- summary_stat(bootstrap_samples[, feature], sd)

  # Calculate confidence interval
  if (CI) {
    n <- length(bootstrap_samples$Sample_numb)

    t <- stats::qt(p = 1 - (alpha / 2), df = n - 1)

    result$LB <- result$Proportion - t * result$SD / sqrt(n)

    result$UB <- result$Proportion + t * result$SD / sqrt(n)
  }

  return(result)
}

#' Helper function to extract summary statistic, to be used within bootstrap_results.
#'
#' @param data_frame_list Feature column of the bootrstrap_samples object created by fci_bootstrap
#' @param stat Statistics to be extracted, i.e. mean, sd
#'
#' @returns Results data frame
#'
#' @examples # summary_stat(bootstrap_samples[, feature], sd)

summary_stat <- function(data_frame_list, stat) {

  # Extract value column from data.frame list
  data_vector_list <- lapply(
    data_frame_list,
    function(x) {
      return(x[, 3])
    }
  )

  # Combine value vectors into single data.frame
  stat_matrix <- as.data.frame(do.call(cbind, data_vector_list))

  # Apply function to each row of data.frame
  result <- apply(stat_matrix, 1, stat)

  return(result)
}
