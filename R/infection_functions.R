#' Calculate age-specific entomological inoculation rate

#' @param age Vector of simulation age/time (in days) over a sufficiently long period, e.g. 1:(365 * 10).
#' @param eir EIR per day.
#' @param rho Age dependent biting parameter. Default = 0.85.
#' @param a0 Age dependent biting parameter. Default = 2920 days.
#' @param gamma_llin Adjustment to exposure for individuals sleeping under a bednet. Default = 1 (no adjustment).
#'
#' @export
get_epsilon <- function(age, eir, rho = 0.85, a0 = 2920, gamma_llin = 1){
  gamma_llin * eir *  (1 - rho * exp(-age / a0))
}

#' Calculate the risk of infection
#'
#' @param epsilon Vector of age-specific EIR.
#' @param b Probability of infection upon mosquito bite.
#' @param seasonality Seasonality profile (vector of seasonality per day). Default = NULL - no seasonality. Can also be generated using `get_season`.
#' @param vaccine_efficacy Vaccine efficacy against infection over time (proportion). Default = NULL - no vaccine.
#'
#' @export
get_infection_hazard <- function(epsilon, b, seasonality = NULL, vaccine_efficacy = NULL){

  if(!is.null(seasonality) & length(seasonality) != length(epsilon)) {
    "seasonality and epsilon must have the same length"
  }

  if(!is.null(vaccine_efficacy) & length(vaccine_efficacy) != length(epsilon)) {
    "vaccine_efficacy and epsilon must have the same length"
  }

  if(any(seasonality < 0)) {
    stop("All values in season must be positive")
  }

  ih <- epsilon * b

  if(!is.null(seasonality)){
    ih <- ih * seasonality
  }

  if(!is.null(vaccine_efficacy)){
    ih <- ih * (1 - vaccine_efficacy)
  }

  return(ih)
}

#' Calculate the risk of clinical malaria in the absence of biting heterogeneity
#'
#' @param infection_hazard Infection hazard over time.
#' @param clinical_probability Probability of developing clinical symptoms.
#'
get_clinical_hazard_no_het <- function(infection_hazard, clinical_probability){
  infection_hazard * clinical_probability
}

#' Simulate average clinical hazard for a cohort
#'
#' Simulates clinical incidence in a cohort of children over time since birth (starting at age 0).
#' Clinical incidence is returned from enrollment into the trial, defined by `age_at_enrollment`.
#' `age_at_enrollment` also determines the time at which the vaccine effect represented in `vx` starts.
#' The seasonality pattern in `season` needs to start at trial enrollment, and gets repeated
#' during the simulation if shorter than the simulation period defined by `age`.
#'
#' @param eir EIR per day.
#' @param age Vector of simulation age/time (in days) over a sufficiently long period, e.g. 1:(365 * 10). Cannot be longer than vx.
#' @param age_at_enrollment Age at enrollment into trial (days).
#' @param gamma_llin Adjustment to exposure for individuals sleeping under a bednet. Set to 1 if not adjusting for this.
#' @param vx Vector of vaccine efficacy against infection over time (proportion).
#' @param season Seasonality profile (vector of seasonality per day). Can be generated using `get_season`.
#' @param icm Maternal immunity level.
#' @param n Number of heterogeneity groups (in mosquito bite exposure).
#' @param zeta Heterogeneity parameters.
#' @param weight Heterogeneity weights.
#' @param rho Age dependent biting parameter (0.85).
#' @param a0 Age dependent biting parameter (2920 days).
#' @param ub Period in which pre-erythrocytic immunity cannot be boosted (7.2 days).
#' @param db Average duration of pre-erythrocytic immunity (3650 days).
#' @param b0 Maximum probability of infection due to no pre-erythrocytic immunity (0.59).
#' @param b1 Maximum reduction in probability of infection due to pre-erythrocytic immunity (0.5).
#' @param ib0 Scale parameter for pre-erythrocytic immunity (43.9).
#' @param kb Shape parameter for pre-erythrocytic immunity (2.16).
#' @param uc Period in which clinical immunity cannot be boosted (6.06 days).
#' @param dc Average duration of clinical immunity (10950 days).
#' @param phi0 Maximum probability of clinical disease due to no clinical immunity (0.792).
#' @param phi1 Maximum reduction in probability of clinical disease due to clinical immunity (0.00074).
#' @param ic0 Scale parameter for clinical immunity (18).
#' @param kc Shape parameter for clinical immunity (2.37).
#' @param cpp Use cpp functions. Default = TRUE.

#' @export
get_clinical_hazard <- function(eir, age, age_at_enrollment, gamma_llin, vx,
                                season, icm, n, zeta, weight,
                                rho, a0, ub, db, b0, b1, ib0, kb, uc, dc,
                                phi0, phi1, ic0, kc, cpp = TRUE) {



  if(length(vx) < length(age)) {
    stop("vx must be the same length as the age vector or longer")
  }


  if(max(age) < age_at_enrollment) {
    stop("age_at_enrollment must be within age vector")
  }

  if(length(season) %% 365 != 0){
    stop("The length of the season vector must be a multiple of 365")
  }

  if(n < 1 | n != floor(n)) {
    stop("n must be an integer >= 1")
  }


  time <- round(age - age_at_enrollment)
  vx_shift <- c(rep(0, sum(time < 0)), vx)[age]
  seasonality <-  season[(time %% length(season)) + 1]

  clinical_hazard <- rep(0, length(which(time > 0)))

  for(h in 1:n){

    # Exposure to infectious bites
    epsilon <- get_epsilon(age = age, eir = eir * zeta[h], rho = rho, a0 = a0,
                           gamma_llin = gamma_llin)

    # Immunity against malaria infection
    ib <- acquire_immunity_r(exposure = epsilon, u = ub, d = db)

    # Probability infection
    if(cpp){
      b <- get_b_cpp(ib = ib, b0 = b0, b1 = b1, ib0 = ib0, kb = kb)
    } else {
      b <- get_b_r(ib = ib, b0 = b0, b1 = b1, ib0 = ib0, kb = kb)
    }

    # Hazard of infection
    infection_hazard <- get_infection_hazard(epsilon = epsilon, b = b,
                                             seasonality = seasonality,
                                             vaccine_efficacy = vx_shift)

    # Immunity against clinical malaria
    ica <- acquire_immunity_r(exposure = infection_hazard, u = uc, d = dc)

    # Probability an infection is clinical
    if(cpp){
      phi <- get_phi_cpp(ica = ica, icm = icm, phi0 = phi0, phi1 = phi1,
                       ic0 = ic0, kc = kc)
    } else {
      phi <- get_phi_r(ica = ica, icm = icm, phi0 = phi0, phi1 = phi1,
                       ic0 = ic0, kc = kc)
    }
    index_period <- which(time > 0)

    # Hazard of clinical infection
    clinical_hazard <- clinical_hazard +
      get_clinical_hazard_no_het(infection_hazard = infection_hazard,
                                 clinical_probability = phi)[index_period] * weight[h]

  }

  return(clinical_hazard)
}
