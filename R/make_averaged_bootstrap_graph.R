#' Plot averaged partial ancestral graph based on bootstrap analysis
#'
#' Function to create a averaged partial ancestral graph based on the bootstrap analysis.
#' For this first all edges are included that occur in a larger proportion of
#' bootstrap samples than a preset criterion. Then all marks are applied to the
#' ends of these edges that occur in a larger proportion of
#' bootstrap samples than a preset criterion. Lastly, circles are assigned to the remaining edge
#' to indicate missing evidence.
#'
#' @param Arrowheads_results Dataframe holding the results of the bootstrap analysis for the arrowheads.
#' @param Tails_results Dataframe holding the results of the bootstrap analysis for the tails.
#' @param Circles_results Dataframe holding the results of the bootstrap analysis for the circles.
#' @param Edges_results Dataframe holding the results of the bootstrap analysis for the edges.
#' @param sig Boolean to indicate whether only significant edges should be plotted.
#' @param edge_threshold If sig=False removes all edges that occur in a smaller proportion of bootstrap samples than this.
#' @param mark_threshold If sig=False removes all marks that occur in a smaller proportion of bootstrap samples than this.
#'
#' @examples
#'
#' # make_averaged_bootstrap_graph(
#' #   Arrowheads_results = Arrowheads_results_without,
#' #   Tails_results = Tails_results_without,
#' #   Circles_results = Circles_results_without,
#' #   Edges_results = Edges_results_without,
#' #   edge_threshold = .8,
#' #   mark_threshold = .8
#' # )
#'
#' @export


make_averaged_bootstrap_graph <- function(Arrowheads_results, Tails_results, Circles_results,
                                          Edges_results, sig = FALSE, edge_threshold = 0, mark_threshold = 0) {

  if (sig & edge_threshold > 0) {
    message("Since, both significance and edge_threshold have been chosen, the more
            conservative edge_threshold will be used and the significance argument will be ignored.")
    sig <- FALSE
  }

  # Remove insignificant Edges
  if (sig) {
    Edges_results_insig <- Edges_results[Edges_results$LB <= 0, c(1, 2)]
    Edges_results_insig <- data.frame(
      `Vertex 1` = c(Edges_results_insig$`Vertex 1`, Edges_results_insig$`Vertex 2`),
      `Vertex 2` = c(Edges_results_insig$`Vertex 2`, Edges_results_insig$`Vertex 1`)
    )
    Arrowheads_results[do.call(paste0, Arrowheads_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig) |
                         do.call(paste0, Arrowheads_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig), c(3, 4, 5, 6)] <- rep(0, 4)
    Tails_results[do.call(paste0, Tails_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig) |
                    do.call(paste0, Tails_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig), c(3, 4, 5, 6)] <- rep(0, 4)
    Circles_results[do.call(paste0, Circles_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig) |
                      do.call(paste0, Circles_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig), c(3, 4, 5, 6)] <- rep(0, 4)
  } else if (edge_threshold > 0) { # Remove edges below the edge_threshold
    Edges_results_insig <- Edges_results[Edges_results$LB <= edge_threshold, c(1, 2)]
    Edges_results_insig <- data.frame(
      `Vertex 1` = c(Edges_results_insig$`Vertex 1`, Edges_results_insig$`Vertex 2`),
      `Vertex 2` = c(Edges_results_insig$`Vertex 2`, Edges_results_insig$`Vertex 1`)
    )
    Arrowheads_results[do.call(paste0, Arrowheads_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig) |
                         do.call(paste0, Arrowheads_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig), c(3, 4, 5, 6)] <- rep(0, 4)
    Tails_results[do.call(paste0, Tails_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig) |
                    do.call(paste0, Tails_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig), c(3, 4, 5, 6)] <- rep(0, 4)
    Circles_results[do.call(paste0, Circles_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig) |
                      do.call(paste0, Circles_results[, c(1, 2)]) %in% do.call(paste0, Edges_results_insig), c(3, 4, 5, 6)] <- rep(0, 4)
  }

  Arrowheads_results <- Arrowheads_results[, c(1, 2, 3, 5)]
  Tails_results <- Tails_results[, c(1, 2, 3, 5)]
  Circles_results <- Circles_results[, c(1, 2, 3, 5)]

  colnames(Arrowheads_results) <- c("From", "Into", "Proportion_Arrowheads", "LB_Arrowheads")
  colnames(Tails_results) <- c("From", "Into", "Proportion_Tails", "LB_Tails")
  colnames(Circles_results) <- c("From", "Into", "Proportion_Circles", "LB_Circles")

  overall_results <- dplyr::full_join(Arrowheads_results, Tails_results, by = c("From", "Into")) %>%
    dplyr::full_join(Circles_results, by = c("From", "Into"))

  if (sig) {
    overall_results$Mark <- mapply(
      function(Proportion_Arrowheads, Proportion_Tails, Proportion_Circles, LB_Arrowheads, LB_Tails, LB_Circles) {
        if (Proportion_Arrowheads == 0 & Proportion_Tails == 0 & Proportion_Arrowheads == 0) {
          return(0)
        } else if (Proportion_Arrowheads > Proportion_Tails & Proportion_Arrowheads > Proportion_Circles & LB_Arrowheads > 0) {
          return(2)
        } else if (Proportion_Tails > Proportion_Arrowheads & Proportion_Tails > Proportion_Circles & LB_Tails > 0) {
          return(3)
        } else if (Proportion_Circles >= Proportion_Arrowheads & Proportion_Circles >= Proportion_Tails & LB_Circles > 0) {
          return(1)
        } else if (LB_Arrowheads <= 0 & LB_Tails <= 0) {
          return(1)
        } else if (LB_Arrowheads > 0 & LB_Tails > 0 & Proportion_Arrowheads == Proportion_Tails) {
          return(1)
        }
      },
      Proportion_Arrowheads = overall_results$Proportion_Arrowheads,
      Proportion_Tails = overall_results$Proportion_Tails,
      Proportion_Circles = overall_results$Proportion_Circles,
      LB_Arrowheads = overall_results$LB_Arrowheads,
      LB_Tails = overall_results$LB_Tails,
      LB_Circles = overall_results$LB_Circles
    )
  } else if (mark_threshold > 0) {
    overall_results$Mark <- mapply(
      function(Proportion_Arrowheads, Proportion_Tails, Proportion_Circles, LB_Arrowheads, LB_Tails, LB_Circles, mark_threshold) {
        if (Proportion_Arrowheads == 0 & Proportion_Tails == 0 & Proportion_Arrowheads == 0) {
          return(0)
        } else if (Proportion_Arrowheads > Proportion_Tails & Proportion_Arrowheads > Proportion_Circles & LB_Arrowheads > mark_threshold) {
          return(2)
        } else if (Proportion_Tails > Proportion_Arrowheads & Proportion_Tails > Proportion_Circles & LB_Tails > mark_threshold) {
          return(3)
        } else if (Proportion_Circles >= Proportion_Arrowheads & Proportion_Circles >= Proportion_Tails & LB_Circles > mark_threshold) {
          return(1)
        } else if (LB_Arrowheads <= mark_threshold & LB_Tails <= mark_threshold) {
          return(1)
        } else if (LB_Arrowheads > mark_threshold & LB_Tails > mark_threshold & Proportion_Arrowheads == Proportion_Tails) {
          return(1)
        }
      },
      Proportion_Arrowheads = overall_results$Proportion_Arrowheads,
      Proportion_Tails = overall_results$Proportion_Tails,
      Proportion_Circles = overall_results$Proportion_Circles,
      LB_Arrowheads = overall_results$LB_Arrowheads,
      LB_Tails = overall_results$LB_Tails,
      LB_Circles = overall_results$LB_Circles,
      MoreArgs = list(mark_threshold = mark_threshold)
    )
  } else {
    overall_results$Mark <- mapply(
      function(Proportion_Arrowheads, Proportion_Tails, Proportion_Circles) {
        if (Proportion_Arrowheads == 0 & Proportion_Tails == 0 & Proportion_Arrowheads == 0) {
          return(0)
        } else if (Proportion_Arrowheads > Proportion_Tails & Proportion_Arrowheads > Proportion_Circles) {
          return(2)
        } else if (Proportion_Tails > Proportion_Arrowheads & Proportion_Tails > Proportion_Circles) {
          return(3)
        } else if (Proportion_Circles >= Proportion_Arrowheads & Proportion_Circles >= Proportion_Tails) {
          return(1)
        } else if (Proportion_Arrowheads == Proportion_Tails) {
          return(1)
        }
      },
      Proportion_Arrowheads = overall_results$Proportion_Arrowheads,
      Proportion_Tails = overall_results$Proportion_Tails,
      Proportion_Circles = overall_results$Proportion_Circles
    )
  }

  overall_results <- overall_results[, c(1, 2, 9)]

  overall_results <- as.data.frame(overall_results)

  overall_results[, c(1, 2)] <- sapply(overall_results[, c(1, 2)], as.character)

  overall_results[, c(3)] <- sapply(overall_results[, c(3)], as.numeric)

  overall_results <- reshape2::dcast(overall_results, From ~ Into, value.var = "Mark")

  overall_results <- overall_results[, !(names(overall_results) %in% "From")]

  row.names(overall_results) <- colnames(overall_results)

  overall_results[is.na(overall_results)] <- 0

  overall_results <- as.matrix(overall_results)

  overall_results[overall_results == 0 & t(overall_results) > 0] <- 1

  # Overall results was transposed, are from and to were the wrong way around?
  print(plot_expected_graph(t(overall_results), labels = colnames(overall_results)))
}
