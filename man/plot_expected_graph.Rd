% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_expected_graph.R
\name{plot_expected_graph}
\alias{plot_expected_graph}
\title{Plot Graph from Amat Matrix}
\usage{
plot_expected_graph(acyc_graph_amat, labels)
}
\arguments{
\item{acyc_graph_amat}{Adjacency matrix of the graph that should be plotted, in the style of a pcalg pag amat matrix.}

\item{labels}{Character vector of variable names or node labels.}
}
\description{
Function to create a target or expected pag based on a pcalg pag amat style adjacency matrix.
Wrapper for Rgraphviz
}
