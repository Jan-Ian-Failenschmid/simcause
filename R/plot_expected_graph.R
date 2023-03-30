#' Plot Graph from Amat Matrix
#'
#' Function to create a target or expected pag based on a pcalg pag amat style adjacency matrix.
#' Wrapper for Rgraphviz
#'
#' @param acyc_graph_amat Adjacency matrix of the graph that should be plotted, in the style of a pcalg pag amat matrix.
#' @param labels Character vector of variable names or node labels.
#'
#' @export


plot_expected_graph <- function(acyc_graph_amat, labels) {
  rownames(acyc_graph_amat) <- colnames(acyc_graph_amat) <- labels # Add lables for cosmetic reasons
  g <- methods::as(acyc_graph_amat, "graphNEL")
  nn <- graph::nodes(g)
  p <- graph::numNodes(g)
  n.edges <- graph::numEdges(g) / 2
  ahs <- ats <- rep("none", n.edges)
  nms <- character(n.edges)
  cmat <- array(c(
    "0" = "none", "1" = "odot",
    "2" = "normal", "3" = "none"
  )[as.character(acyc_graph_amat)],
  dim = dim(acyc_graph_amat), dimnames = dimnames(acyc_graph_amat)
  )
  iE <- 0L
  for (i in seq_len(p - 1)) {
    x <- nn[i]
    for (j in (i + 1):p) {
      y <- nn[j]
      if (acyc_graph_amat[x, y] != 0) {
        iE <- iE + 1L
        ahs[[iE]] <- cmat[x, y]
        ats[[iE]] <- cmat[y, x]
        nms[[iE]] <- paste0(x, "~", y)
      }
    }
  }
  names(ahs) <- names(ats) <- nms
  graph::edgeRenderInfo(g) <- list(arrowhead = ahs, arrowtail = ats)
  Rgraphviz::renderGraph(Rgraphviz::layoutGraph(g))
}
