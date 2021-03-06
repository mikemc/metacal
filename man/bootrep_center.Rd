% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compositional-mean.R
\name{bootrep_center}
\alias{bootrep_center}
\title{Generate bootstrap replicates of the sample center}
\usage{
bootrep_center(
  .data,
  R = 4000,
  N = nrow(.data),
  method = "proj",
  dist = "dirichlet",
  in_scale = "linear",
  out_scale = "linear",
  denom = NULL,
  enframe = FALSE
)
}
\arguments{
\item{.data}{All-numeric data frame or matrix with taxa as columns.}

\item{R}{Number of bootstrap replicates.}

\item{N}{Number of trials for multinomial resampling.}

\item{method}{Method for computing the center: "proj", "gm", or "rss".}

\item{dist}{Distribution for drawing the bootstrap weights: "dirichlet" or
"multinomial".}

\item{in_scale}{"linear" (default) or "log".}

\item{out_scale}{"linear" (default) or "log".}

\item{denom}{Taxa to use in the denominator; if NULL, use all taxa.}

\item{enframe}{Whether to "enframe" the bootstrap estimates in a tibble.}
}
\description{
Generate bootstrap replicates of the sample center
}
