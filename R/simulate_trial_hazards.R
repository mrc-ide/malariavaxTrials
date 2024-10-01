#' Simulate clinical hazard in trial cohort
#'
#' Simulate clinical hazard in the vaccine and control arms of a trial.
#' Option to adjust for ITN use in the trial (@param gamma_llin).
#' Clinical incidence is returned individually for the vaccine arms with and without bednets
#' and the control arm with and without bednets to allow fitting to invididual-level
#' data by ITN use.
#'
#' @param eir EIR per day.
#' @param age_at_enrollment Age at enrollment into trial (days).
#' @param gamma_llin Adjustment to exposure for individuals sleeping under a bednet.
#' @param vx Vaccine efficacy over time (proportion). Must be a list of vectors. Can contain multiple vectors for different vaccine arms in the trial (e.g. different doses, booster schedules)
#' @param r_clin Adjustment for case definition.
#' @param age Vector of age/time (in days) over a sufficiently long time, e.g. 1:(365 * 12).
#' @param n Number of heterogeneity groups (in mosquito bite exposure)
#' @param season Seasonality profile. Default = rep(1, 365) - no seasonality.
#' @param overrides a named list of parameter values to use instead of defaults.
#' @param cpp Use cpp functions. Default = TRUE
#'
#' @return List of simulated clinical hazards for each trial arm cohort
#'
#' @export
simulate_trial_hazards <- function(eir, age_at_enrollment, gamma_llin, vx, r_clin,
                                   age, n = 10, season = rep(1, 365),
                                   overrides = list(s2 = 1.16), cpp = TRUE) {

  if(class(vx) != "list") {
    stop("vx must be a list and contain at least 1 vector")
  }

  if(any(sapply(vx, length)<length(age))) {
    stop("vx vectors must be the same length as the age vector or longer")
  }

  # Load fitted model parameters
    url <- "https://raw.github.com/mrc-ide/malariaEquilibrium/master/inst/extdata/Jamie_parameters.rds"
    p <- readRDS(gzcon(url(url)))

  # Override parameters with any client specified ones
  if (!is.list(overrides)) {
    stop('overrides must be a list')
  }

  for (name in names(overrides)) {
    if (!(name %in% names(p))) {
      stop(paste('unknown parameter', name, sep=' '))
    }
    p[[name]] <- overrides[[name]]
  }


  # Heterogeneity
  gh <- statmod::gauss.quad.prob(n = n, dist = "normal")
  zeta <- exp(-p$s2 * 0.5 + sqrt(p$s2) * gh$nodes)
  weight <- gh$weights

  # Maternally acquired immunity
  ica20_no_itn <- get_maternal_start(eir = eir, gamma_llin = 1, season = season,
                                     rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db,
                                     b0 = p$b0, b1 = p$b1, ib0 = p$IB0, kb = p$kb,
                                     uc = p$uc, dc = p$dc, cpp = cpp)

  ica20_itn <- get_maternal_start(eir = eir, gamma_llin = gamma_llin, season = season,
                                  rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db,
                                  b0 = p$b0, b1 = p$b1, ib0 = p$IB0, kb = p$kb,
                                  uc = p$uc, dc = p$dc, cpp = cpp)

  icm_no_itn <- get_icm(age = age, ica20 = ica20_no_itn, pm = p$PM, dm = p$dm)
  icm_itn <- get_icm(age = age, ica20 = ica20_itn, pm = p$PM, dm = p$dm)

  # Simulate clinical hazard

  n_vaccine_arms <- length(vx)
  vax_no_itn <- list()
  vax_itn <- list()


  for (i in 1:n_vaccine_arms) {

    # For vaccine cohort, without bednets
    vax_no_itn[[i]] <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                               gamma_llin = 1, vx = vx[[i]], season = season,
                                               icm = icm_no_itn, age = age, n = n, zeta = zeta, weight = weight,
                                               rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                                               b1 = p$b1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                                               dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                                               kc = p$kc, cpp = cpp)

    # For vaccine cohort, with bednets
    vax_itn[[i]] <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                            gamma_llin = gamma_llin, vx = vx[[i]], season = season,
                                            icm = icm_itn, age = age, n = n, zeta = zeta, weight = weight,
                                            rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                                            b1 = p$b1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                                            dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                                            kc = p$kc, cpp = cpp)
  }

  names(vax_no_itn) <- paste0(rep("vax_no_itn_arm", n_vaccine_arms), 1:n_vaccine_arms)
  names(vax_itn) <- paste0(rep("vax_itn_arm", n_vaccine_arms), 1:n_vaccine_arms)

  # For unvaccinated cohort, without bednets
  control_no_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                            gamma_llin = 1, vx = rep(0, length(age)),
                                            season = season, icm = icm_no_itn, age = age,
                                            n = n, zeta = zeta, weight = weight,
                                            rho = p$rho, a0 = p$a0, ub = p$ub,
                                            db = p$db, b0 = p$b0, b1 = p$b1,
                                            ib0 = p$IB0, kb = p$kb, uc = p$uc, dc = p$dc,
                                            phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0, kc = p$kc, cpp = cpp)

  control_no_itn <- list(control_no_itn = control_no_itn)

  # For unvaccinated cohort, with bednets
  control_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                         gamma_llin = gamma_llin, vx = rep(0, length(age)),
                                         season = season, icm = icm_itn, age = age,
                                         n = n, zeta = zeta, weight = weight,
                                         rho = p$rho, a0 = p$a0, ub = p$ub,
                                         db = p$db, b0 = p$b0, b1 = p$b1,
                                         ib0 = p$IB0, kb = p$kb, uc = p$uc, dc = p$dc,
                                         phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0, kc = p$kc, cpp = cpp)

  control_itn <- list(control_itn = control_itn)

  out <- append(append(append(vax_itn, vax_no_itn), control_itn), control_no_itn)

  return(out)

}




