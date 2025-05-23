% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/immunity_functions.R
\name{acquire_immunity_r}
\alias{acquire_immunity_r}
\title{Immunity function}
\usage{
acquire_immunity_r(exposure, u, d)
}
\arguments{
\item{exposure}{Vector of immunity-boosting exposure over time. Can be age-dependent EIR (epsilon) or infection hazard.}

\item{u}{Period in which immunity cannot be boosted (days since last patent infection).}

\item{d}{Average duration of immunity (days).}
}
\value{
Immunity over time.
}
\description{
Acquisition of immunity against infection (IB) and clinical malaria (ICA).
}
