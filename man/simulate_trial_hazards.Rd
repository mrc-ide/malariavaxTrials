% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_trial_hazards.R
\name{simulate_trial_hazards}
\alias{simulate_trial_hazards}
\title{Simulate clinical or severe hazard in trial cohort}
\usage{
simulate_trial_hazards(
  outcome = "clinical",
  eir,
  age,
  age_at_enrollment,
  gamma_llin,
  vx,
  r_clin,
  n = 10,
  season = rep(1, 365),
  overrides = list(s2 = 1.16),
  cpp = TRUE
)
}
\arguments{
\item{outcome}{"clinical" for clinical hazard or "severe" for severe malaria (default = "clinical").}

\item{eir}{EIR per day.}

\item{age}{Vector of simulation age/time (in days) over a sufficiently long period, e.g. 1:(365 * 10). Cannot be longer than vectors in vx.}

\item{age_at_enrollment}{Age at enrollment into trial (days).}

\item{gamma_llin}{Adjustment to exposure for individuals sleeping under a bednet. Set to 1 if not adjusting for this.}

\item{vx}{Vaccine efficacy against infection over time (proportion). Must be a list of vectors. Can contain multiple vectors for different vaccine arms in the trial (e.g. different doses, booster schedules).}

\item{r_clin}{Adjustment for case definition.}

\item{n}{Number of heterogeneity groups (in mosquito bite exposure). Default = 10.}

\item{season}{Seasonality profile (vector of seasonality per day). Default = rep(1, 365) - no seasonality. Can also be generated using \code{get_season}.}

\item{overrides}{A named list of parameter values to use instead of defaults. Can be used to change malaria model parameters.}

\item{cpp}{Use cpp functions. Default = TRUE.}
}
\value{
List of simulated clinical hazards for each trial arm cohort
}
\description{
Simulate clinical/severe hazard in the vaccine and control arms of a trial.
Clinical/severe incidence is returned individually for the vaccine arms with and without bednets
and the control arm with and without bednets to allow fitting to invididual-level
data by ITN use. Option to adjust for ITN use in the trial using \code{gamma_llin}.
}
