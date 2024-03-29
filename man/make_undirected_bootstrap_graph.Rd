% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_undirected_bootstrap_graph.R
\name{make_undirected_bootstrap_graph}
\alias{make_undirected_bootstrap_graph}
\title{Plot undirected graph based on the bootstrap graph feature analysis}
\usage{
make_undirected_bootstrap_graph(bootstrap_edges, sig = FALSE, threshold = 0)
}
\arguments{
\item{bootstrap_edges}{Dataframe holding the results of the bootstrap analysis for the edges.}

\item{sig}{Boolean to indicate whether only significant edges should be plotted.}

\item{threshold}{If sig=False removes all edges that occur in a smaller proportion of bootstrap samples than this.}
}
\description{
Create undirected graph in which only edges are included that have been assessed
against a predefined criterion in the bootstrap analysis
}
\examples{

# make_undirected_bootstrap_graph(Edges_results)

}
