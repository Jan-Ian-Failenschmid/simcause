#' Plot results of the bootstrap graph feature analysis
#'
#' Plots bootstrap confidence intervals for each feature identified in the graph
#' and compares it against a predefined crtiterion.
#'
#' @param bootstrap_results_df Dataframe holding the results of the bootstrap analysis.
#' @param feature Graph feature that should be plotted
#' @param criterion Criterion proportion to compare the bootstrap proportion against.
#'
#' @examples
#'
#' # plot_bootstrap_results(Edges, "Edges")
#'
#' @export

plot_bootstrap_results <- function(bootstrap_results_df, feature, criterion = .8) {

  # Add Ylab column
  if (feature == "Edges") {
    bootstrap_results_df$Ylab <- paste0(bootstrap_results_df[, 1], " - ", bootstrap_results_df[, 2])
  } else if (feature == "Arrowheads") {
    bootstrap_results_df$Ylab <- paste0(bootstrap_results_df[, 1], " *-> ", bootstrap_results_df[, 2])
  } else if (feature == "Circles") {
    bootstrap_results_df$Ylab <- paste0(bootstrap_results_df[, 1], " *-o ", bootstrap_results_df[, 2])
  } else if (feature == "Tails") {
    bootstrap_results_df$Ylab <- paste0(bootstrap_results_df[, 1], " *-- ", bootstrap_results_df[, 2])
  }

  # Add significance indicator
  bootstrap_results_df$color[bootstrap_results_df$LB > 0] <- 1
  bootstrap_results_df$color[bootstrap_results_df$LB > criterion] <- 2
  bootstrap_results_df$color[is.na(bootstrap_results_df$color)] <- 0

  # Make plot
  print(ggplot2::ggplot(
    data = bootstrap_results_df,
    mapping = aes(
      x = Proportion,
      y = Ylab,
      xmin = LB,
      xmax = UB,
      color = as.factor(color)
    )
  ) +
    geom_point() +
    geom_errorbarh(height = 0) +
    geom_vline(xintercept = c(0, criterion)) +
    scale_color_manual(
      breaks = c(0, 1, 2),
      values = c("black", "blue", "red")
    ) +
    theme_minimal() +
    scale_x_continuous(breaks = seq(0, 1, .1)) +
    labs(title = paste0(feature, " Bootstrap Results")) +
    theme(
      axis.line.x.bottom = element_line(arrow = arrow(length = unit(0.25, "cm"))),
      legend.position = "none",
      axis.ticks.x.bottom = element_line(),
      axis.ticks.length = unit(-0.25, "cm"),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 4)
    ))
}
