% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fci_precision.R
\name{make_plot}
\alias{make_plot}
\title{Helper function to style plots, to be used within plot_results.}
\usage{
make_plot(
  r = "Edge_sensitivity",
  b = 0.5,
  plot_num = "m",
  plot_line = TRUE,
  layout = "complete",
  plot_data_subs = NULL
)
}
\arguments{
\item{r}{Type of plot to be styled}

\item{b}{Mean value for a normal distribution from which edge slopes are drawn in this plot.}

\item{plot_num}{Indicator whether the plots of different effect sizes should be "m"erged or "s"eperate.}

\item{plot_line}{Logical whether loess regression line should be added to result plot.}

\item{layout}{Layout for the plots illustrating the results of the simulation. Options are: "complete_m", "cycle_m", "dag_m", "complete_s", "cycle_s", "dag_s", "edges_s".}

\item{plot_data_subs}{Obtional argument for adding error bars to the graph. Currently, disabled as the grphs get too cluttered.}
}
\description{
Helper function to style plots, to be used within plot_results.
}
