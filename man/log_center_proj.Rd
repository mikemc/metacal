% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compositional-mean.R
\name{log_center_proj}
\alias{log_center_proj}
\title{Compute the log center by the projection method}
\usage{
log_center_proj(mat, weights = rep(1, nrow(mat)))
}
\description{
Compute the log center by the projection method, which allows for missing
observations.
}
\references{
van den Boogaart KG, Tolosana-Delgado R, Bren M. 2006. Concepts
for handling of zeros and missing values in compositional data. Proc IAMG
6:1-4. http://www.stat.boogaart.de/Publications/iamg06_s07_01.pdf

Bren M, Tolosana-Delgado R, van den Boogaart KG. 2008. News from
"compositions", the R package.
https://core.ac.uk/download/pdf/132548286.pdf
}
