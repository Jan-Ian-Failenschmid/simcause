#' Bootstrap Sampling for the Fast Causal Inference Algorithm
#'
#' Bootstrap analysis of FCI results for a given data set.
#' Can be used to quantify the uncertainty that is present for any graph feature within the data.
#'
#' @param data Data frame of the variables to be analysed.
#' @param labels Character vector of variable names or node labels.
#' @param amount Amount of bootstrap samples to be drawn.
#' @param sample_proportion Relative sample size of the bootstrap samples to the observed data.
#' @param replace Logical indicating whether resampling should be done with or without replacement. Needs to be true for sample_proportion > 1.
#' @param alpha Alpha value for significance tests within FCI.
#'
#' @returns Returns data frame holding all bootstrap samples
#'
#' @examples
#'
#' # bootstrap_n <- 200
#' # sample_proportion <- 1
#' # threshold <- 1
#' # bootstrap_samples <- fci_bootstrap(
#' #   data = data, labels = colnames(data), amount = bootstrap_n,
#' #   sample_proportion = sample_proportion
#' # )
#'
#' @export

fci_bootstrap <- function(data,
                          labels,
                          amount = 100,
                          sample_proportion = 1,
                          replace = TRUE,
                          alpha = 0.01) {

  # If sample_proportion is larger than 1, replace has to be TRUE
  if (sample_proportion > 1 & replace == FALSE) replace <- TRUE

  # Precalculate sample size to safe operations later
  sample_size <- round(nrow(data) * sample_proportion, digits = 0)

  # Create Grid for the correct amount of samples
  bootstrap_samples <- data.frame(
    Sample_numb = c(1:amount)
  )

  # Draw bootstrap samples
  bootstrap_samples$Samples <- lapply(
    bootstrap_samples$Sample_numb,
    function(Sample_numb, x = x, size = size, replace = replace) {
      x[sample(x = nrow(x), size = size, replace = replace), ]
    },
    x = data,
    size = sample_size,
    replace = replace
  )

  # Specify suffStat for each sample
  bootstrap_samples$suffStat <- lapply(
    bootstrap_samples$Samples,
    function(Samples, n = n) {
      list(C = cor(Samples), n = n)
    },
    n = sample_size
  )

  # Run FCI over the grid
  bootstrap_samples$Fit_FCI <- lapply(
    bootstrap_samples$suffStat,
    pcalg::fci,
    indepTest = gaussCItest,
    alpha = alpha,
    labels = labels,
    selectionBias = FALSE
  )

  # Extract amat matrix from FCI output
  bootstrap_samples$amat <- lapply(
    bootstrap_samples$Fit_FCI,
    function(Fit_FCI) {
      return(Fit_FCI@amat)
    }
  )

  # Melt amat matrix into long format
  bootstrap_samples$melted <- lapply(
    bootstrap_samples$amat,
    function(amat) {
      melted <- reshape2::melt(amat)
      colnames(melted) <- c("From", "Into", "Mark")
      melted <- melted[melted$From != melted$Into, ]
      return(melted)
    }
  )

  # Index edges and marks for each sample

  bootstrap_samples$Edges <- lapply(
    bootstrap_samples$melted,
    function(melted) {
      Edges <- melted[!duplicated(t(apply(melted[, c("From", "Into")], 1, sort))), ]
      colnames(Edges) <- c("Vertex 1", "Vertex 2", "Edge")
      Edges$Edge[Edges$Edge > 0] <- 1
      return(Edges)
    }
  )

  bootstrap_samples$Arrowheads <- lapply(
    bootstrap_samples$melted,
    function(melted) {
      Arrowheads <- melted
      Arrowheads$Mark[Arrowheads$Mark != 2] <- 0
      Arrowheads$Mark[Arrowheads$Mark == 2] <- 1
      colnames(Arrowheads) <- c("From", "Into", "Arrowhead")
      return(Arrowheads)
    }
  )

  bootstrap_samples$Tails <- lapply(
    bootstrap_samples$melted,
    function(melted) {
      Tails <- melted
      Tails$Mark[Tails$Mark != 3] <- 0
      Tails$Mark[Tails$Mark == 3] <- 1
      colnames(Tails) <- c("From", "Into", "Tail")
      return(Tails)
    }
  )

  bootstrap_samples$Circles <- lapply(
    bootstrap_samples$melted,
    function(melted) {
      Circles <- melted
      Circles$Mark[Circles$Mark != 1] <- 0
      Circles$Mark[Circles$Mark == 1] <- 1
      colnames(Circles) <- c("From", "Into", "Circles")
      return(Circles)
    }
  )

  return(bootstrap_samples)
}
