% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infection_functions.R
\name{get_epsilon}
\alias{get_epsilon}
\title{Calculate age-specific entomological inoculation rate}
\usage{
get_epsilon(age, eir, rho = 0.85, a0 = 2920, gamma_llin = 1)
}
\arguments{
\item{age}{Vector of simulation age/time (in days) over a sufficiently long period, e.g. 1:(365 * 10).}

\item{eir}{EIR per day.}

\item{rho}{Age dependent biting parameter. Default = 0.85.}

\item{a0}{Age dependent biting parameter. Default = 2920 days.}

\item{gamma_llin}{Adjustment to exposure for individuals sleeping under a bednet. Default = 1 (no adjustment).}
}
\description{
Calculate age-specific entomological inoculation rate
}
