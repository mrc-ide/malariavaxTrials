% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infection_functions.R
\name{get_severe_hazard}
\alias{get_severe_hazard}
\title{Simulate average severe hazard for a cohort}
\usage{
get_severe_hazard(
  eir,
  age,
  age_at_enrollment,
  gamma_llin,
  vx,
  season,
  ivm,
  n,
  zeta,
  weight,
  rho,
  a0,
  ub,
  db,
  b0,
  b1,
  ib0,
  kb,
  uv,
  dv,
  theta0,
  theta1,
  iv0,
  kv,
  fv0,
  av,
  gammav,
  cpp = TRUE
)
}
\arguments{
\item{eir}{EIR per day.}

\item{age}{Vector of simulation age/time (in days) over a sufficiently long period, e.g. 1:(365 * 10). Cannot be longer than vx.}

\item{age_at_enrollment}{Age at enrollment into trial (days).}

\item{gamma_llin}{Adjustment to exposure for individuals sleeping under a bednet. Set to 1 if not adjusting for this.}

\item{vx}{Vector of vaccine efficacy against infection over time (proportion).}

\item{season}{Seasonality profile (vector of seasonality per day). Can be generated using \code{get_season}.}

\item{ivm}{Maternal severe immunity level.}

\item{n}{Number of heterogeneity groups (in mosquito bite exposure).}

\item{zeta}{Heterogeneity parameters.}

\item{weight}{Heterogeneity weights.}

\item{rho}{Age dependent biting parameter (0.85).}

\item{a0}{Age dependent biting parameter (2920 days).}

\item{ub}{Period in which pre-erythrocytic immunity cannot be boosted (7.2 days).}

\item{db}{Average duration of pre-erythrocytic immunity (3650 days).}

\item{b0}{Maximum probability of infection due to no pre-erythrocytic immunity (0.59).}

\item{b1}{Maximum reduction in probability of infection due to pre-erythrocytic immunity (0.5).}

\item{ib0}{Scale parameter for pre-erythrocytic immunity (43.9).}

\item{kb}{Shape parameter for pre-erythrocytic immunity (2.16).}

\item{uv}{Period in which severe immunity cannot be boosted (11.43 days).}

\item{dv}{Average duration of severe immunity (10950 days).}

\item{theta0}{Maximum probability of severe malaria due to no severe immunity (0.075).}

\item{theta1}{Maximum reduction in probability of severe malaria due to severe immunity (0.0001)}

\item{iv0}{Scale parameter for severe immunity (1.096).}

\item{kv}{Shape parameter for severe immunity (2.00).}

\item{fv0}{Age-dependent modifier (0.141).}

\item{av}{Age-dependent modifier (2493.41).}

\item{gammav}{Age-dependent modifier (2.91).}

\item{cpp}{Use cpp functions. Default = TRUE.}
}
\description{
Simulates severe malaria incidence in a cohort of children over time since birth (starting at age 0).
Severe malaria incidence is returned from enrollment into the trial, defined by \code{age_at_enrollment}.
\code{age_at_enrollment} also determines the time at which the vaccine effect represented in \code{vx} starts.
The seasonality pattern in \code{season} needs to start at trial enrollment, and gets repeated
during the simulation if shorter than the simulation period defined by \code{age}.
}
