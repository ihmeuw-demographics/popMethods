% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popReconstruct_settings.R
\name{transform_dt}
\alias{transform_dt}
\title{Helper function to transform values in a data.table}
\usage{
transform_dt(dt, value_col, transformation, transformation_arguments)
}
\arguments{
\item{dt}{[\code{data.table(1)}]\cr
The data.table to apply the transformation functions to.}

\item{value_col}{[\code{character(1)}]\cr
The column in \code{dt} to transform.}

\item{transformation}{[\verb{function(1)}]\cr
NULL if no transformation. Otherwise transformation function to
apply to \code{value_col}.}

\item{transformation_arguments}{[\code{list()}]\cr
NULL if no transformation. Otherwise list of arguments to provide to the
corresponding \code{transformation} function(s).}
}
\value{
invisibly return \code{dt} with transformed values.
}
\description{
Helper function to transform values in a data.table
}
\examples{
popMethods:::transform_dt(
  dt = data.table::data.table(year = 1950, value = 2.5),
  value_col = "value",
  transformation = demUtils::logit,
  transformation_arguments = list(domain_lower = 0, domain_upper = 5)
)
}
