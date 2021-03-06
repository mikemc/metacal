% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compositional-mean.R
\name{center}
\alias{center}
\title{Compute the center (compositional mean) of a set of compositions}
\usage{
center(
  .data,
  weights = rep(1, nrow(.data)),
  method = "proj",
  in_scale = "linear",
  out_scale = "linear",
  denom = NULL,
  enframe = FALSE,
  components = FALSE
)
}
\arguments{
\item{.data}{All-numeric data frame or matrix with taxa as columns.}

\item{weights}{Sample (row) weights.}

\item{method}{Method for computing the center: "proj", "gm", or "rss".}

\item{in_scale}{"linear" or "log".}

\item{out_scale}{"linear" or "log".}

\item{denom}{Taxa to use in the denominator; if NULL, use all taxa.}

\item{enframe}{Whether to return the bias estimate as a two-column tibble.}
}
\description{
Unobserved values should be marked as NaN.
}
