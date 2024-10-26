% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segments.R
\name{get_segments}
\alias{get_segments}
\title{Calculate phase segments.}
\usage{
get_segments(
  phases,
  ma.win = ceiling(nrow(phases)) * 0.02,
  spar = 0.001,
  cut = TRUE,
  plot = FALSE,
  verb = 0
)
}
\description{
Calculate phase segments.
}
