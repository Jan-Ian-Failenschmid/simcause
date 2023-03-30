#' Simulate data from a potentially cyclical, directed, linear and Gaussian system at
#' different sample sizes. Apply the FCI algorithm to each sample and analyse
#' how often each feature of the graph has been identified correctly at each sample size.
#'
#' @param true_graph The causal structure of the system that is simulated from in form of a directed adjacency matrix.
#' @param lables Character vector of variable names or node labels.
#' @param beta Mean value for a normal distribution from which edge slopes are drawn.
#' @param target Target partial ancestral graph that optimally should be returned by FCI in the form of a pcalg pag amat matrix.
#' @param cycle_nodes Character vector of nodes that are part of a cycle.
#' @param alpha Alpha value for significance tests within FCI.
#' @param nmax Maximum sample size for the simulation.
#' @param nmin Minimum sample size for the simulation.
#' @param steps Steps of sample sizes to be simulated between the minimum and maximum.
#' @param plot_line Logical whether loess regression line should be added to result plot.
#' @param n_iter Amount of samples drawn at each sample size.
#' @param layout Layout for the plots illustrating the results of the simulation. Options are: "complete_m", "cycle_m", "dag_m", "complete_s", "cycle_s", "dag_s", "edges_s".
#' @param verbose Verbose argument for the FCI function
#'
#' @returns A list with a data frame holding the complete simulation and a data frame holding all results at each sample size. Additionally, plots the simulation results.
#'
#' @examples
#' # Nr of variables
#' p <- 5
#'
#'
#' # True SCM
#' true_graph <- matrix(c(
#'   0, 1, 0, 0, 0,
#'   0, 0, 1, 0, 0,
#'   0, 0, 0, 1, 0,
#'   0, 0, 0, 0, 1,
#'   0, 1, 0, 0, 0
#' ), ncol = p, byrow = TRUE)
#'
#' # Labels
#' labels <- c("A", "B", "C", "D", "E")
#'
#' # Define Betas to Go Through
#' beta_graph <- c(.5, 1.1, 2, 3, 5)
#'
#' # Define Sample Range
#' nmin <- 20
#' nmax <- 800
#' steps <- 20
#' n_iter <- 80
#'
#' # Define Cycle Nodes
#' cycle_nodes <- c("B", "C", "D", "E")
#'
#' # Define Ploting Parameters
#' plot_line <- TRUE
#' layout <- "cycle_m"
#'
#' target <- matrix(c(
#'   0, 2, 0, 0, 2,
#'   1, 0, 1, 0, 1,
#'   0, 1, 0, 1, 0,
#'   0, 0, 1, 0, 1,
#'   1, 1, 0, 1, 0
#' ), ncol = p, byrow = TRUE)
#'
#' # Run Simulation
#' fci_precision(
#'   true_graph = true_graph, labels = labels,
#'   beta = beta_graph, target = target, nmin = nmin, nmax = nmax,
#'   steps = steps, n_iter = n_iter, cycle_nodes = cycle_nodes,
#'   plot_line = plot_line, layout = layout
#' )

fci_precision <- function(true_graph,
                          labels,
                          beta,
                          target,
                          cycle_nodes = NULL,
                          alpha = 0.05,
                          nmax = 800,
                          nmin = 50,
                          steps = 20,
                          n_iter = 30,
                          plot_line = TRUE,
                          layout = "complete_m",
                          verbose = FALSE) {

  if
  (!require(pcalg, quietly = T, warn.conflicts = F)) {
    stop("pcalg is not installed")
  }

  if (!require(tidyverse, quietly = T, warn.conflicts = F)) {
    stop("tidyverse is not installed")
  }

  if (!require(patchwork, quietly = T, warn.conflicts = F)) {
    stop("patchwork is not installed")
  }

  if (!require(ggrepel, quietly = T, warn.conflicts = F)) {
    stop("ggrepel is not installed")
  }

  check_input(
    true_graph = true_graph,
    labels = labels,
    beta = beta,
    target = target,
    cycle_nodes = cycle_nodes,
    alpha = alpha,
    nmax = nmax,
    nmin = nmin,
    steps = steps,
    n_iter = n_iter,
    plot_line = plot_line,
    layout = layout,
    verbose = verbose
  )

  # Create simulation grid
  simulation_grid <- expand.grid(
    Beta = beta,
    N = round(seq(from = nmin, to = nmax, length.out = steps), digits = 0),
    Iteration = seq(1, n_iter)
  )

  # Draw samples markov to the true graphs accros the simulation grid
  simulation_grid$samples <- mapply(
    make_samples,
    simulation_grid$Beta,
    simulation_grid$N,
    MoreArgs = list(true_graph = true_graph),
    SIMPLIFY = F
  )

  # Specify suffStat for each sample
  simulation_grid$suffStat <- mapply(function(X, n) {
    list(C = cor(X), n = n)
  },
  X = simulation_grid$samples,
  n = simulation_grid$N,
  SIMPLIFY = F
  )

  # Define independence test
  indepTest <- gaussCItest

  # Run FCI over the grid
  simulation_grid$Fit_FCI <- lapply(
    simulation_grid$suffStat,
    fci,
    indepTest = indepTest,
    alpha = alpha,
    labels = labels,
    selectionBias = FALSE
  )

  # Extract graph comparisons
  simulation_grid$Results <- lapply(
    simulation_grid$Fit_FCI,
    compare_graphs,
    target = target,
    cycle_nodes = cycle_nodes,
    labels = labels
  )

  # Reformat results for output
  results_output <- format_results(simulation_grid = simulation_grid)

  # Plot results
  plot_results(results_output = results_output, plot_line = plot_line, layout = layout)

  # Return full simulation data frames on request
  if (verbose) {
    return(list(simulation_grid, results_output))
  }
}
