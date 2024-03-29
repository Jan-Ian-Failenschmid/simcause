#' Plot undirected graph based on the bootstrap graph feature analysis
#'
#' Create undirected graph in which only edges are included that have been assessed
#' against a predefined criterion in the bootstrap analysis
#'
#' @param bootstrap_edges Dataframe holding the results of the bootstrap analysis for the edges.
#' @param sig Boolean to indicate whether only significant edges should be plotted.
#' @param threshold If sig=False removes all edges that occur in a smaller proportion of bootstrap samples than this.
#'
#' @examples
#'
#' # make_undirected_bootstrap_graph(Edges_results)
#'
#' @export

make_undirected_bootstrap_graph <- function(bootstrap_edges, sig = FALSE, threshold = 0) {

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
  graph_matrix <- reshape2::dcast(bootstrap_edges, `Vertex 1` ~ `Vertex 2`, value.var = "Proportion")
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

  qgraph::qgraph(graph_matrix,
         layout = "spring", theme = "gray",
         nodeNames = row.names(graph_matrix), minimum = 0, maximum = 1
  )
}
