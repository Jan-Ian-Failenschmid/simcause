#' Plots undirected bootsrap graph in which the proportion of times an eges has been
#' included corresponds to the saturation of the edge.
#'
#' @param bootstrap_edges Dataframe holding the proportion of time each edge has been included created by bootstrap_results(_, "Edges")
#' @param sig Logical if true only edges are included that are identified significantly more often than 0 times.
#' @param threshold Value between 0 and 1. If sig is false only edges are included with a proportion significantly greater than the threshold.
#'
#' @examples
#' undirected_bootstrap_graph(Edges_results, threshold = .8)

undirected_bootstrap_graph <- function(bootstrap_edges, sig = FALSE, threshold = 0) {
  if (!require(qgraph, quietly = T, warn.conflicts = F)) {
    stop("qgraph is not installed")
  }

  if (!require(reshape2, quietly = T, warn.conflicts = F)) {
    stop("reshape2 is not installed")
  }

  # Remove insignificant edges
  if (sig) {
    bootstrap_edges$Proportion[bootstrap_edges$LB <= 0] <- 0
  }

  # Remove edges under a threshold
  if (threshold > 0) {
    bootstrap_edges$Proportion[bootstrap_edges$LB <= threshold] <- 0
  }

  # Getting the data back into a qgraph compatible form is more work than anticipated
  # Subset relevant columns
  bootstrap_edges <- bootstrap_edges[, c(1, 2, 3)]
  # Reformat Node names to character
  bootstrap_edges[, c(1, 2)] <- lapply(bootstrap_edges[, c(1, 2)], as.character)
  # Add the nodes to each list that miss from it
  add <- unique(bootstrap_edges$`Vertex 2`[which(!bootstrap_edges$`Vertex 2` %in% bootstrap_edges$`Vertex 1`)])
  bootstrap_edges <- rbind(bootstrap_edges, c(as.character(add), as.character(add), 0))
  add <- unique(bootstrap_edges$`Vertex 1`[which(!bootstrap_edges$`Vertex 1` %in% bootstrap_edges$`Vertex 2`)])
  bootstrap_edges <- rbind(bootstrap_edges, c(as.character(add), as.character(add), 0))

  # A lot more work
  # Widen the data list to a matrix
  graph_matrix <- dcast(bootstrap_edges, `Vertex 1` ~ `Vertex 2`, value.var = "Proportion")
  # Change NA's to 0
  graph_matrix[is.na(graph_matrix)] <- 0
  # Remove Vertex 1 column
  graph_matrix <- graph_matrix[, !(names(graph_matrix) %in% "Vertex 1")]
  # Change Matrix to numeric
  graph_matrix <- sapply(graph_matrix, as.numeric)
  # Reformat as actual matrix
  graph_matrix <- as.matrix(graph_matrix)
  # Add transpose to make matrix symetrical
  graph_matrix <- graph_matrix + t(graph_matrix)
  # Add rownames
  row.names(graph_matrix) <- colnames(graph_matrix)

  qgraph(graph_matrix,
         layout = "spring", theme = "gray",
         nodeNames = row.names(graph_matrix), minimum = 0, maximum = 1
  )
}
