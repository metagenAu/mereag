% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta_anaylsis.R
\name{calc_meta_effect}
\alias{calc_meta_effect}
\title{Calculate Meta-Analysis Effect Size using Random Effects Model}
\usage{
calc_meta_effect(data, sm = "RR")
}
\arguments{
\item{data}{A data frame containing the meta-analysis data.}

\item{sm}{A character string specifying the summary measure; default is 'RR' (relative risk).}
}
\value{
A data frame with columns: RandomEffect, p_value, LCL (Lower Confidence Limit), UCL (Upper Confidence Limit), and n_studies (number of studies included).
}
\description{
This function performs a meta-analysis using a random effects model. It computes the random effect size, its p-value, and the lower and upper confidence limits.
}
\examples{
# Example usage
# calc_meta_effect(data = my_data, sm = 'OR')
}
