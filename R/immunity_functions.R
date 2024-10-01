#' Immunity function
#'
#' Acquisition of immunity against infection (IB) and clinical
#' malaria (ICA). White et al (S26).
#'
#' @param exposure Vector of immunity-boosting exposure over time. Can be age-dependent EIR (epsilon) or infection hazard.
#' @param u Period in which immunity cannot be boosted (days since last patent infection).
#' @param d Average duration of immunity (days).
#' @return Immunity over time.

#' @export
acquire_immunity_r <- function(exposure, u, d) {

  n <- length(exposure)
  immunity <- numeric(n)
  immunity[1] <- 0

  for (i in 2:n) {
    immunity[i] <- immunity[i - 1] + (exposure[i] / (exposure[i] * u + 1)) - (immunity[i - 1] / d)
  }

  return(immunity)
}

#' Maternal immunity function
#'
#' Maternal immunity is a function of the average immunity levels of a 20-year old woman.
#' @param age Vector of age/time (in days) over a sufficiently long time, e.g. 1:(365 * 12).
#' @param ica20 Average level of immunity in a 20 year old women.
#' @param pm Maternal immunity relative to mother's immunity level.
#' @param dm Average duration of maternal immunity.
#' @return Maternal immunity over time.

#' @export
get_icm <- function(age, ica20, pm, dm){
  pm * ica20 * exp(- age / dm)
}

#' Estimate the average level of immunity in a 20 year old women
#' to inform the maternal immunity
#'
#' @param eir EIR per day.
#' @param gamma_llin Adjustment to exposure for individuals sleeping under a bednet.
#' @param season Seasonality profile.
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
#' @param cpp Use cpp functions. Default = TRUE
#'
#' @export
get_maternal_start <- function(eir, gamma_llin,
                               season, rho, a0, ub, db, b0, b1, ib0, kb, uc, dc, cpp = TRUE){

  age <- 0:(365*20)
  vx <- rep(0, length(age))
  seasonality <-  season[(age %% length(season)) + 1]

  # Exposure to infectious bites
  epsilon <- get_epsilon(age = age, eir = eir, rho = rho, a0 = a0,
                         gamma_llin = gamma_llin)

  # Immunity against malaria infection
  if(cpp){
    ib <- acquire_immunity_cpp(exposure = epsilon, u = ub, d = db)
  } else {
    ib <- acquire_immunity_r(exposure = epsilon, u = ub, d = db)
  }

  # Probability infection
  if(cpp){
    b <- get_b_cpp(ib = ib, b0 = b0, b1 = b1, ib0 = ib0, kb = kb)
  } else {
    b <- get_b_r(ib = ib, b0 = b0, b1 = b1, ib0 = ib0, kb = kb)
  }

  # Hazard of infection
  infection_hazard <- get_infection_hazard(epsilon = epsilon, b = b,
                                           seasonality = seasonality,
                                           vaccine_efficacy = vx)

  # Immunity against clinical malaria
  if(cpp){
    ica <- acquire_immunity_cpp(exposure = infection_hazard, u = uc, d = dc)
  } else {
    ica <- acquire_immunity_r(exposure = infection_hazard, u = uc, d = dc)
  }


  return(ica[length(age)])
}

#' Calculate probability of an infection leading to a clinical malaria episode (phi)
#'
#' @param ica Clinical immunity level.
#' @param icm Maternal immunity level.
#' @param phi0 Maximum probability of clinical disease due to no clinical immunity (0.792).
#' @param phi1 Maximum reduction in probability of clinical disease due to clinical immunity (0.00074)
#' @param ic0 Scale parameter for clinical immunity (18).
#' @param kc Shape parameter for clinical immunity (2.37).
#'
#' @export
get_phi_r <- function(ica, icm, phi0, phi1, ic0, kc){
  phi0 * (phi1 + ((1 - phi1) / (1 + ((ica + icm) / ic0)^ kc)))
}

#' Calculate probability of infection upon receiving an infectious bite (b)
#'
#' @param ib Pre-erythrocytic immunity level.
#' @param b0 Maximum probability of infection due to no pre-erythrocytic immunity (0.59).
#' @param b1 Maximum reduction in probability of infection due to pre-erythrocytic immunity (0.5).
#' @param ib0 Scale parameter for pre-erythrocytic immunity (43.9).
#' @param kb Shape parameter for pre-erythrocytic immunity (2.16).

#' @export
get_b_r <- function(ib, b0, b1, ib0, kb){
  b0 * (b1 + ((1 - b1) / (1 + (ib / ib0)^ kb)))
}
