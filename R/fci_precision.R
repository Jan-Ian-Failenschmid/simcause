#' Sensitivity Simulations for Directed Linear and Gaussian Causal Systems
#'
#' Simulate data from a potentially cyclical, directed, linear and Gaussian system at
#' different sample sizes. Apply the FCI algorithm to each sample and analyse
#' how often each feature of the graph has been identified correctly at each sample size.
#'
#' @param true_graph The causal structure of the system that is simulated from in form of a directed adjacency matrix.
#' @param labels Character vector of variable names or node labels.
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
#' # # Nr of variables
#' # p <- 5
#' #
#' #
#' # # True SCM
#' # true_graph <- matrix(c(
#' #   0, 1, 0, 0, 0,
#' #   0, 0, 1, 0, 0,
#' #   0, 0, 0, 1, 0,
#' #   0, 0, 0, 0, 1,
#' #   0, 1, 0, 0, 0
#' # ), ncol = p, byrow = TRUE)
#' #
#' # # Labels
#' # labels <- c("A", "B", "C", "D", "E")
#' #
#' # # Define Betas to Go Through
#' # beta_graph <- c(.5, 1.1, 2, 3, 5)
#' #
#' # # Define Sample Range
#' # nmin <- 20
#' # nmax <- 800
#' # steps <- 20
#' # n_iter <- 80
#' #
#' # # Define Cycle Nodes
#' # cycle_nodes <- c("B", "C", "D", "E")
#' #
#' # # Define Ploting Parameters
#' # plot_line <- TRUE
#' # layout <- "cycle_m"
#' #
#' # target <- matrix(c(
#' #   0, 2, 0, 0, 2,
#' #   1, 0, 1, 0, 1,
#' #   0, 1, 0, 1, 0,
#' #   0, 0, 1, 0, 1,
#' #   1, 1, 0, 1, 0
#' # ), ncol = p, byrow = TRUE)
#' #
#' # # Run Simulation
#' # fci_precision(
#' #   true_graph = true_graph, labels = labels,
#' #   beta = beta_graph, target = target, nmin = nmin, nmax = nmax,
#' #   steps = steps, n_iter = n_iter, cycle_nodes = cycle_nodes,
#' #   plot_line = plot_line, layout = layout
#' # )
#'
#' @export

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


  # check_input(
  #   true_graph = true_graph,
  #   labels = labels,
  #   beta = beta,
  #   target = target,
  #   cycle_nodes = cycle_nodes,
  #   alpha = alpha,
  #   nmax = nmax,
  #   nmin = nmin,
  #   steps = steps,
  #   n_iter = n_iter,
  #   plot_line = plot_line,
  #   layout = layout,
  #   verbose = verbose
  # )

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
    list(C = stats::cor(X), n = n)
  },
  X = simulation_grid$samples,
  n = simulation_grid$N,
  SIMPLIFY = F
  )

  # Run FCI over the grid
  simulation_grid$Fit_FCI <- lapply(
    simulation_grid$suffStat,
    pcalg::fci,
    indepTest = gaussCItest,
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


#' Helper function to check the input of fci_precision and return appropriate error messages.
#' Does nothing else.

check_input <- function(true_graph,
                        labels,
                        beta,
                        target,
                        cycle_nodes,
                        alpha,
                        nmax,
                        nmin,
                        steps,
                        n_iter,
                        plot_line,
                        layout,
                        verbose) {
  if (!is.matrix(true_graph) | !is.numeric(true_graph)) {
    stop("Argument true_graph must be a numeric matrix")
  }

  if (!is.vector(labels) | !is.character(labels)) {
    stop("Argument labels must be a character vector")
  }

  if (!is.vector(beta) | !is.numeric(beta)) {
    stop("Argument beta must be a numeric vector")
  }

  if (!is.matrix(target) | !is.numeric(target)) {
    stop("Argument target must be a numeric matrix")
  }

  if (!is.null(cycle_nodes)) {
    if (!is.vector(cycle_nodes) | !is.character(cycle_nodes)) {
      stop("Argument labels must be a character vector or NULL")
    }

    if (!length(unique(cycle_nodes)) == length(cycle_nodes)) {
      stop("All elements of the argument cycle_nodes must be unique")
    }
  }

  if (!is.numeric(alpha)) {
    stop("Argument alpha must be numeric")
  }

  if (!is.numeric(nmax)) {
    stop("Argument nmax must be numeric")
  }

  if (!is.numeric(nmin)) {
    stop("Argument nmin must be numeric")
  }

  if (!is.numeric(steps)) {
    stop("Argument steps must be numeric")
  }

  if (!is.numeric(n_iter)) {
    stop("Argument n_iter must be numeric")
  }

  if (!is.logical(plot_line) | !length(plot_line) == 1) {
    stop("Argument plot_line must be either TRUE or FALSE")
  }

  if (!layout %in% c("complete_s", "dag_s", "cycle_s", "edges_s", "complete_m", "dag_m", "cycle_m") |
      !length(layout) == 1) {
    stop("Argument layout is not one of the available options, choose one of the list:
          complete_s, dag_s, cycle_s, edges_s, complete_m, dag_m, cycle_m")
  }

  if (!is.logical(verbose) | !length(verbose) == 1) {
    stop("Argument verbose must be either TRUE or FALSE")
  }

  if (!dim(true_graph)[1] == dim(true_graph)[2]) {
    stop("Argument true_graph must be a square matrix")
  }

  if (!dim(target)[1] == dim(target)[2]) {
    stop("Argument target must be a square matrix")
  }

  if (!dim(true_graph) == dim(target)) {
    stop("Argument true_graph and argument target must have the same dimension")
  }

  if (!length(unique(labels)) == length(labels)) {
    stop("All elements of the argument labels must be unique")
  }

  if (is.null(cycle_nodes) & layout %in% c("complete_s", "cycle_s", "complete_m", "cycle_m")) {
    stop("If no cycle_nodes are specified layout must be either edges or dag")
  }

  rownames(true_graph) <- colnames(true_graph) <- labels

  if (sum(true_graph[cycle_nodes, cycle_nodes]) == 0 & layout %in% c("complete_s", "dag_s", "complete_m", "dag_m")) {
    stop("Arguemnt true_graph must have at least one edge that is not in or connected
         to the cycle in order to choose the layout complete or dag")
  }
}

#' Helper function to sample from a potentially cyclic, linear, Gassian graph,
#' to be used within fci_precision.
#'
#' @param Beta Mean value for a normal distribution from which edge slopes are drawn.
#' @param N Sample size
#' @param true_graph The causal structure of the system that is simulated from in form of a directed adjacency matrix.
#'
#' @returns A data frame with a sample based on the true graph.

make_samples <- function(Beta,
                         N,
                         true_graph) {

  # Create beta matrix to mirror true graph
  beta_matrix <- true_graph

  # Sample beta candidates from the standard normal
  beta_candidates <- stats::rnorm(length(beta_matrix[true_graph == 1]))

  # Ensure that beta values have a mean of beta
  beta_candidates <- beta_candidates - mean(beta_candidates) + Beta

  # Assign beta candidates to the beta matrix
  beta_matrix[true_graph == 1] <- beta_candidates

  # Calculate the equilibrium Solution of Beta
  IminBinv <- solve(diag(1, ncol(true_graph)) - beta_matrix)

  # Sample data from the standard normal
  Z <- matrix(rnorm(N * ncol(true_graph), sd = 1),
              ncol = ncol(true_graph),
              nrow = N
  )

  # Transform data to be markov to the graph
  X <- Z %*% IminBinv

  return(X)
}

#' Helper function to compare a graph generated by FCI with an expected optimal graph,
#' to be used within fci_precision
#'
#' @param probe The output of fci.
#' @param target Target partial ancestral graph that optimally should be returned by FCI in the form of a pcalg pag amat matrix.
#' @param cycle_nodes Character vector of nodes that are part of a cycle.
#' @param lables Character vector of variable names or node labels.
#'
#' @returns Returns a list with feature sensitivities and orientation precision.

compare_graphs <- function(probe,
                           target,
                           cycle_nodes,
                           labels) {

  # Extract graph from fci output
  probe <- probe@amat

  # Make sure that graph labels are applied to target

  rownames(target) <- colnames(target) <- labels

  # Distinguish Cycle and Dag parts of the Graph
  cycle_probe <- probe[cycle_nodes, cycle_nodes]

  cycle_target <- target[cycle_nodes, cycle_nodes]

  dag_probe <- probe[!labels %in% cycle_nodes, !labels %in% cycle_nodes]

  dag_target <- target[!labels %in% cycle_nodes, !labels %in% cycle_nodes]

  target[labels %in% cycle_nodes, !labels %in% cycle_nodes]

  target[!labels %in% cycle_nodes, labels %in% cycle_nodes]

  probe[labels %in% cycle_nodes, !labels %in% cycle_nodes]

  probe[!labels %in% cycle_nodes, labels %in% cycle_nodes]

  # Count true positives and false negatives
  true_positive_edges <- sum(target > 0 & probe > 0) / 2

  false_negative_edges <- sum(target > 0 & probe == 0) / 2

  true_positive_circle_dag <- sum(dag_target == 1 & dag_probe == 1)

  false_negative_circle_dag <- sum(dag_target == 1 & dag_probe != 1)

  true_positive_arrowheads_dag <- sum(dag_target == 2 & dag_probe == 2)

  false_negative_arrowheads_dag <- sum(dag_target == 2 & dag_probe != 2)

  true_positive_tails_dag <- sum(dag_target == 3 & dag_probe == 3)

  false_negative_tails_dag <- sum(dag_target == 3 & dag_probe != 3)

  # Check cycle orientation
  cycle_no_orient <- all(cycle_probe <= 1)

  cycle_clockwise_orient <-
    sum(lower.tri(cycle_probe) & cycle_probe == 2) == (length(cycle_nodes) - 1) &
    sum(upper.tri(cycle_probe) & cycle_probe == 2) == 1

  cycle_counterclockwise_orient <-
    sum(upper.tri(cycle_probe) & cycle_probe == 2) == (length(cycle_nodes) - 1) &
    sum(lower.tri(cycle_probe) & cycle_probe == 2) == 1

  # Calculate sensitivity
  edge_sensitivity <- true_positive_edges / (true_positive_edges + false_negative_edges)

  circle_sensitivity <- true_positive_circle_dag / (true_positive_circle_dag + false_negative_circle_dag)

  arrowhead_sensitivity <- true_positive_arrowheads_dag / (true_positive_arrowheads_dag + false_negative_arrowheads_dag)

  tail_sensitivity <- true_positive_tails_dag / (true_positive_tails_dag + false_negative_tails_dag)

  return(list(
    edge_sensitivity = edge_sensitivity,
    circle_sensitivity = circle_sensitivity,
    arrowhead_sensitivity = arrowhead_sensitivity,
    tail_sensitivity = tail_sensitivity,
    cycle_no_orient = cycle_no_orient,
    cycle_clockwise_orient = cycle_clockwise_orient,
    cycle_counterclockwise_orient = cycle_counterclockwise_orient
  ))
}

#' Helper function to format results within fci_precision.
#'
#' @param simulation_grid Data frame holding the simulation results created within fci_precision.
#'
#' @returns A list of two data frames holding the simulation results in long format for plotting.

format_results <- function(simulation_grid) {
  # Prepare Data for plotting
  result_data <- data.frame(
    Beta = simulation_grid$Beta,
    N = simulation_grid$N,
    Edge_sensitivity = unlist(purrr::map(simulation_grid$Results, 1)),
    Circle_sensitivity = unlist(purrr::map(simulation_grid$Results, 2)),
    Arrowhead_sensitivity = unlist(purrr::map(simulation_grid$Results, 3)),
    Tail_sensitivity = unlist(purrr::map(simulation_grid$Results, 4)),
    Cycle_no_orient = unlist(purrr::map(simulation_grid$Results, 5)),
    Cycle_clockwise_orient = unlist(purrr::map(simulation_grid$Results, 6)),
    Cycle_counterclockwise_orient = unlist(purrr::map(simulation_grid$Results, 7))
  )

  # Transform plot_data to long format
  result_data_long <- dplyr::summarise(dplyr::group_by(tidyr::gather(data = result_data, key = "Result", value = "value", 3:9), Beta, N, Result),
    avg = mean(value),
    sd = sd(value)
  )


  result_data_long_complete <- tidyr::gather(data = result_data, key = "Result", value = "value", 3:9)

  return(list(result_data_long, result_data_long_complete))
}

#' Helper function to plot simulation results, to be used within fci_precision
#'
#' @param results_output List of data frames holding the formated simulation results, created by format_results
#' @param plot_line Logical whether loess regression line should be added to result plot.
#' @param layout Layout for the plots illustrating the results of the simulation. Options are: "complete_m", "cycle_m", "dag_m", "complete_s", "cycle_s", "dag_s", "edges_s".

plot_results <- function(results_output,
                         plot_line = TRUE,
                         layout = "complete_m") {

  requireNamespace("ggpubr", quietly = TRUE)

  requireNamespace("patchwork", quietly = TRUE)

  # Split up complete and summarized data sets
  plot_data_complete <- results_output[[2]] # complete data set will only be used for
  # the merged smooth stat

  plot_data <- results_output[[1]]

  # Split up layout argument for flow logic s = seperate, m = merged
  layout <- strsplit(layout, split = "_", fixed = TRUE)

  if (length(layout[[1]] == 2)) {
    plot_num <- layout[[1]][[2]]
  } else {
    plot_num <- "m"
  }

  layout <- layout[[1]][[1]]

  # Plot logic for separate plots for the different beta values
  if (plot_num == "s") {

    # Plot logic for the three different multi-aspect plot layouts
    if (layout %in% c("complete", "cycle", "dag")) {

      # Create Seperate plot for each Beta
      for (b in unique(plot_data$Beta)) {

        # Create subplots
        for (r in unique(plot_data$Result)) {

          # Subset data for each plot
          plot_data_subs <- plot_data[plot_data$Beta == b & plot_data$Result == r, ]

          # Create subplot
          assign(
            paste0("plot_", r),
            ggplot2::ggplot(
              data = plot_data_subs,
              mapping = aes(x = N, y = avg, ymin = avg - sd, ymax = avg + sd)
            ) +
              make_plot(r = r, plot_num = plot_num, layout = layout, plot_line = plot_line)
          )
        }

        # Create Beta label subplot
        plot_beta_label <- ggplot2::ggplot() +
          geom_text(aes(label = paste0("Beta = ", b), x = 1, y = 1), size = 24) +
          theme_void(
            base_size = 11,
            base_family = "",
            base_line_size = base_size / 22,
            base_rect_size = base_size / 22
          )

        if (layout == "complete") {
          # Compose complete plot
          layout_complete <- "
          ABB
          CDE
          FGH
        "

          # Arrange Subplots
          print(plot_beta_label + plot_Edge_sensitivity +
                  plot_Arrowhead_sensitivity + plot_Circle_sensitivity + plot_Tail_sensitivity +
                  plot_Cycle_no_orient + plot_Cycle_clockwise_orient + plot_Cycle_counterclockwise_orient +
                  patchwork::plot_layout(design = layout_complete))
        } else if (layout == "cycle") {
          # Compose cycle plot
          layout_cycle <- "
          ABB
          CDE
        "

          # Arrange Subplots
          print(plot_beta_label + plot_Edge_sensitivity +
                  plot_Cycle_no_orient + plot_Cycle_clockwise_orient + plot_Cycle_counterclockwise_orient +
                  patchwork::plot_layout(design = layout_cycle))
        } else if (layout == "dag") {

          # Compose cycle plot
          layout_dag <- "
          ABB
          CDE
        "
          # Arrange Subplots
          print(plot_beta_label + plot_Edge_sensitivity +
                  plot_Arrowhead_sensitivity + plot_Circle_sensitivity + plot_Tail_sensitivity +
                  patchwork::plot_layout(design = layout_dag))
        }
      }
    }

    if (layout == "edges") {
      for (b in unique(plot_data$Beta)) {

        # Subset data for each plot
        plot_data_subs <- plot_data[plot_data$Beta == b & plot_data$Result == "Edge_sensitivity", ]

        # Create subplot for the edges layout
        assign(
          paste0("plot_Edge_sensitivity_", b),
          ggplot2::ggplot(
            data = plot_data_subs,
            mapping = aes(x = N, y = avg, ymin = avg - sd, ymax = avg + sd)
          ) +
            make_plot(b = b, plot_num = plot_num, layout = layout, plot_line = plot_line)
        )
      }

      print(eval(parse(text = paste(paste0("plot_Edge_sensitivity_", unique(plot_data$Beta)), collapse = " / "))))
    }

    # Plot logic for the merged beta plots
  } else if (plot_num == "m") {

    # Plot logic for the three different multi-aspect plot layouts
    if (layout %in% c("complete", "cycle", "dag")) {
      for (r in unique(plot_data$Result)) {

        # Subset data for each plot
        plot_data_subs <- plot_data[plot_data$Result == r, ]
        plot_data_complete_subs <- plot_data_complete[plot_data_complete$Result == r, ]

        # Create subplot
        assign(
          paste0("plot_", r),
          ggplot2::ggplot(data = plot_data_complete_subs, mapping = aes(
            x = N, y = value, color = as.factor(Beta),
            fill = as.factor(Beta)
          )) +
            make_plot(
              r = r, plot_num = plot_num, layout = layout,
              plot_data_subs = plot_data_subs
            )
        )
      }

      if (layout == "complete") {
        # Compose complete plot
        layout_complete <- "
          AAA
          BCD
          EFG
        "
        # Arrange Subplots
        print(plot_Edge_sensitivity +
                plot_Arrowhead_sensitivity + plot_Circle_sensitivity + plot_Tail_sensitivity +
                plot_Cycle_no_orient + plot_Cycle_clockwise_orient + plot_Cycle_counterclockwise_orient +
                patchwork::plot_layout(design = layout_complete))
      } else if (layout == "cycle") {
        # Compose cycle plot
        layout_cycle <- "
          AB
          CD
        "
        # Arrange Subplots
        print(plot_Edge_sensitivity + plot_Cycle_no_orient +
                plot_Cycle_clockwise_orient + plot_Cycle_counterclockwise_orient +
                patchwork::plot_layout(design = layout_cycle))
      } else if (layout == "dag") {

        # Compose cycle plot
        layout_dag <- "
          AB
          CD
        "
        # Arrange Subplots
        print(plot_Edge_sensitivity + plot_Arrowhead_sensitivity +
                plot_Circle_sensitivity + plot_Tail_sensitivity +
                patchwork::plot_layout(design = layout_dag))
      }
    }
  }
}

#' Helper function to style plots, to be used within plot_results.
#'
#' @param r Type of plot to be styled
#' @param b Mean value for a normal distribution from which edge slopes are drawn in this plot.
#' @param plot_num Indicator whether the plots of different effect sizes should be "m"erged or "s"eperate.
#' @param plot_line Logical whether loess regression line should be added to result plot.
#' @param layout Layout for the plots illustrating the results of the simulation. Options are: "complete_m", "cycle_m", "dag_m", "complete_s", "cycle_s", "dag_s", "edges_s".
#' @param plot_data_subs Obtional argument for adding error bars to the graph. Currently, disabled as the grphs get too cluttered.

make_plot <- function(r = "Edge_sensitivity",
                      b = 0.5,
                      plot_num = "m",
                      plot_line = TRUE,
                      layout = "complete",
                      plot_data_subs = NULL) {
  if (plot_num == "s") {
    if (layout %in% c("complete", "cycle", "dag")) {

      # Plot structure for the seperated complete/cycle/dag layout plot
      list(
        geom_point(),
        geom_errorbar(),
        if (plot_line == TRUE) {
          geom_smooth(formula = y ~ sqrt(x), se = FALSE, method = "loess")
        },
        if (r == "Edge_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Edge Sensitivity",
            title = "Edge Sensitivity Estimates for the complete Graph"
          )
        } else if (r == "Circle_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Circle Sensitivity",
            title = "Circle Sensitivity Estimates for the DAG Part"
          )
        } else if (r == "Arrowhead_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Arrowhead Sensitivity",
            title = "Arrowhead Sensitivity Estimates for the DAG Part"
          )
        } else if (r == "Tail_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Tail Sensitivity",
            title = "Tail Sensitivity Estimates for the DAG Part"
          )
        } else if (r == "Cycle_no_orient") {
          labs(
            x = "Proportion",
            y = "Edge Sensitivity",
            title = "Proportion of Correctily Unoriented Cycles"
          )
        } else if (r == "Cycle_clockwise_orient") {
          labs(
            x = "Sample Size",
            y = "Proportion",
            title = "Proportion of Clockwise Oriented Cycles"
          )
        } else if (r == "Cycle_counterclockwise_orient") {
          labs(
            x = "Sample Size",
            y = "Proportion",
            title = "Proportion of Counterclockwise Oriented Cycles"
          )
        },
        theme_bw(),
        theme(
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10)
        )
      )
    } else if (layout == "edges") {

      # Plot structure for the separated edges layout
      list(
        geom_point(),
        geom_errorbar(),
        if (plot_line == TRUE) {
          geom_smooth(formula = y ~ sqrt(x), se = FALSE, method = "loess")
        },
        labs(
          x = "Sample Size",
          y = "Edge Sensitivity",
          title = paste0("Edge Sensitivity Estimates for the complete Graph, Beta = ", b)
        ),
        theme_bw(),
        theme(
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10)
        )
      )
    }
  } else if (plot_num == "m") {
    if (layout %in% c("cycle", "dag")) {
      # Plot structure for the merged cycle/dag layout plot
      list(
        geom_smooth(method = "loess", formula = y ~ x),
        geom_point(data = plot_data_subs, mapping = aes(
          x = N, y = avg, color = as.factor(Beta),
          fill = as.factor(Beta)
        )),
        if (r == "Edge_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Edge Sensitivity",
            title = "Edge Sensitivity Estimates for the complete Graph"
          )
        } else if (r == "Circle_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Circle Sensitivity",
            title = "Circle Sensitivity Estimates for the DAG Part"
          )
        } else if (r == "Arrowhead_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Arrowhead Sensitivity",
            title = "Arrowhead Sensitivity Estimates for the DAG Part"
          )
        } else if (r == "Tail_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Tail Sensitivity",
            title = "Tail Sensitivity Estimates for the DAG Part"
          )
        } else if (r == "Cycle_no_orient") {
          labs(
            x = "Proportion",
            y = "Edge Sensitivity",
            title = "Proportion of Correctily Unoriented Cycles"
          )
        } else if (r == "Cycle_clockwise_orient") {
          labs(
            x = "Sample Size",
            y = "Proportion",
            title = "Proportion of Clockwise Oriented Cycles"
          )
        } else if (r == "Cycle_counterclockwise_orient") {
          labs(
            x = "Sample Size",
            y = "Proportion",
            title = "Proportion of Counterclockwise Oriented Cycles"
          )
        },
        ggrepel::geom_text_repel(
          data = subset(plot_data_subs, N == max(N)),
          mapping = aes(
            label = paste0("Beta = ", as.factor(Beta)),
            y = avg, x = N
          ),
          size = 7,
          hjust = "left",
          direction = "y",
          nudge_x = 2,
          segment.color = NA,
          show.legend = FALSE
        ),
        # adds errorbars if desired, but clutters up the graph badly
        # geom_errorbar(data = plot_data_subs, mapping = aes(x = N, y = avg, ymin = avg - sd, ymax = avg + sd)),
        scale_x_continuous(
          expand = expansion(mult = c(.05, .16)),
          breaks = plot_data_subs$N
        ),
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)),
        theme_minimal(base_size = 24),
        theme(
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 20),
          legend.position = "none",
          axis.line.x.bottom = element_line(arrow = arrow(length = unit(0.25, "cm"))),
          axis.line.y.left = element_line(),
          axis.ticks.x.bottom = element_line(),
          axis.ticks.y.left = element_line(),
          axis.ticks.length = unit(-0.25, "cm")
        )
      )
    } else if (layout == "complete") {
      list(
        geom_smooth(method = "loess", formula = y ~ x),
        geom_point(data = plot_data_subs, mapping = aes(
          x = N, y = avg, color = as.factor(Beta),
          fill = as.factor(Beta)
        )),
        if (r == "Edge_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Edge Sensitivity",
            title = "Edge Sensitivity Estimates for the complete Graph",
            color = "Beta",
            fill = "Beta"
          )
        } else if (r == "Circle_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Circle Sensitivity",
            title = "Circle Sensitivity Estimates for the DAG Part"
          )
        } else if (r == "Arrowhead_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Arrowhead Sensitivity",
            title = "Arrowhead Sensitivity Estimates for the DAG Part"
          )
        } else if (r == "Tail_sensitivity") {
          labs(
            x = "Sample Size",
            y = "Tail Sensitivity",
            title = "Tail Sensitivity Estimates for the DAG Part"
          )
        } else if (r == "Cycle_no_orient") {
          labs(
            x = "Proportion",
            y = "Edge Sensitivity",
            title = "Proportion of Correctily Unoriented Cycles"
          )
        } else if (r == "Cycle_clockwise_orient") {
          labs(
            x = "Sample Size",
            y = "Proportion",
            title = "Proportion of Clockwise Oriented Cycles"
          )
        } else if (r == "Cycle_counterclockwise_orient") {
          labs(
            x = "Sample Size",
            y = "Proportion",
            title = "Proportion of Counterclockwise Oriented Cycles"
          )
        },
        # adds errorbars if desired, but clutters up the graph badly
        # geom_errorbar(data = plot_data_subs, mapping = aes(x = N, y = avg, ymin = avg - sd, ymax = avg + sd)),
        scale_x_continuous(
          expand = expansion(mult = c(.05, .16)),
          breaks = plot_data_subs$N
        ),
        scale_y_continuous(limits = c(0, 1), expand = c(0, 0)),
        theme_minimal(base_size = 24),
        theme(
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          legend.position = ifelse(r == "Edge_sensitivity", "right", "none"),
          axis.line.x.bottom = element_line(arrow = arrow(length = unit(0.25, "cm"))),
          axis.line.y.left = element_line(),
          axis.ticks.x.bottom = element_line(),
          axis.ticks.y.left = element_line(),
          axis.ticks.length = unit(-0.25, "cm")
        )
      )
    }
  }
}
