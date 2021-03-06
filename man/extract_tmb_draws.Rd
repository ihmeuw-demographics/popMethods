% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popReconstruct_draws.R
\name{extract_tmb_draws}
\alias{extract_tmb_draws}
\title{Extract draws from the popReconstruct model TMB fit}
\usage{
extract_tmb_draws(fit, inputs, value_col, settings, detailed_settings)
}
\arguments{
\item{fit}{[\code{stanfit(1)}] or [\code{sdreport(1)}]\cr
Model fit object returned by \code{\link[=popReconstruct_fit]{popReconstruct_fit()}}.}

\item{inputs}{[\code{list()}]\cr
[\code{data.table()}] for each ccmpp input. Requires 'srb', 'asfr', 'baseline';
mortality estimates provided as just 'survival' or two of 'mx', 'ax', and
'qx'; and migration estimates provided as just 'net_migration' or both
'immigration' and 'emigration'. See \strong{Section: Inputs} for more
information on each of the required inputs.}

\item{value_col}{[\code{character(1)}]\cr
Name of the column containing the value of interest in each of the
\code{inputs}. Default is 'value'.}

\item{settings}{[\code{list()}]\cr
Named list of settings for running the popReconstruct model with. The
required settings are the same as those required for \code{\link[demCore:ccmpp]{demCore::ccmpp()}},
see \strong{Section: Settings} for these. The optional settings specific to
popReconstruct are described in
\strong{Section: Optional popReconstruct Settings}}

\item{detailed_settings}{[\code{list()}]\cr
Detailed settings for each ccmpp input and 'population'.}
}
\description{
Extract draws from the popReconstruct model TMB fit
}
\seealso{
\code{\link[=popReconstruct_posterior_draws]{popReconstruct_posterior_draws()}}
}
