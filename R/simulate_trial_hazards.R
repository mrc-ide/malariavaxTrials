#' Simulate clinical or severe hazard in trial cohort
#'
#' Simulate clinical/severe hazard in the vaccine and control arms of a trial.
#' Clinical/severe incidence is returned individually for the vaccine arms with and without bednets
#' and the control arm with and without bednets to allow fitting to invididual-level
#' data by ITN use. Option to adjust for ITN use in the trial using `gamma_llin`.
#'
#' @param outcome "clinical" for clinical hazard or "severe" for severe malaria (default = "clinical").
#' @param eir EIR per day.
#' @param age Vector of simulation age/time (in days) over a sufficiently long period, e.g. 1:(365 * 10). Cannot be longer than vectors in vx.
#' @param age_at_enrollment Age at enrollment into trial (days).
#' @param gamma_llin Adjustment to exposure for individuals sleeping under a bednet. Set to 1 if not adjusting for this.
#' @param vx Vaccine efficacy against infection over time (proportion). Must be a list of vectors. Can contain multiple vectors for different vaccine arms in the trial (e.g. different doses, booster schedules).
#' @param r_clin Adjustment for case definition.
#' @param n Number of heterogeneity groups (in mosquito bite exposure). Default = 10.
#' @param season Seasonality profile (vector of seasonality per day). Default = rep(1, 365) - no seasonality. Can also be generated using `get_season`.
#' @param overrides A named list of parameter values to use instead of defaults. Can be used to change malaria model parameters.
#' @param cpp Use cpp functions. Default = TRUE.
#'
#' @return List of simulated clinical hazards for each trial arm cohort
#'
#' @export
simulate_trial_hazards <- function(outcome = "clinical", eir, age, age_at_enrollment,
                                   gamma_llin, vx, r_clin,
                                   n = 10, season = rep(1, 365),
                                   overrides = list(s2 = 1.16), cpp = TRUE) {

  if(class(vx) != "list") {
    stop("vx must be a list and contain at least 1 vector")
  }

  if(any(sapply(vx, length)<length(age))) {
    stop("vx vectors must be the same length as the age vector or longer")
  }

  if(max(age) < age_at_enrollment) {
    stop("age_at_enrollment must be within age vector")
  }

  if(length(season) %% 365 != 0) {
    stop("The length of the season vector must be a multiple of 365")
  }

  if(any(season < 0)) {
    stop("All values in season must be positive")
  }


  if(n < 1 | n != floor(n)) {
    stop("n must be an integer >= 1")
  }

  # Load fitted model parameters
    #url <- "https://raw.github.com/mrc-ide/malariaEquilibrium/master/inst/extdata/Jamie_parameters.rds"
    name_full <- system.file("extdata/", "Jamie_parameters.rds", package = 'malariavaxTrials', mustWork = TRUE)
    p <- readRDS(name_full)

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

  if(outcome == "clinical") {

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
      if (gamma_llin != 1) {
        vax_itn[[i]] <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                                     gamma_llin = gamma_llin, vx = vx[[i]], season = season,
                                                     icm = icm_itn, age = age, n = n, zeta = zeta, weight = weight,
                                                     rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                                                     b1 = p$b1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                                                     dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                                                     kc = p$kc, cpp = cpp)
      }

    }

    names(vax_no_itn) <- paste0(rep("vax_no_itn_arm", n_vaccine_arms), 1:n_vaccine_arms)
    if (gamma_llin != 1) {
      names(vax_itn) <- paste0(rep("vax_itn_arm", n_vaccine_arms), 1:n_vaccine_arms)
    }

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
    if (gamma_llin != 1) {
      control_itn <- r_clin * get_clinical_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                                  gamma_llin = gamma_llin, vx = rep(0, length(age)),
                                                  season = season, icm = icm_itn, age = age,
                                                  n = n, zeta = zeta, weight = weight,
                                                  rho = p$rho, a0 = p$a0, ub = p$ub,
                                                  db = p$db, b0 = p$b0, b1 = p$b1,
                                                  ib0 = p$IB0, kb = p$kb, uc = p$uc, dc = p$dc,
                                                  phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0, kc = p$kc, cpp = cpp)

      control_itn <- list(control_itn = control_itn)
    }

    if(gamma_llin != 1) {
      out <- append(append(append(vax_itn, vax_no_itn), control_itn), control_no_itn)
    } else {
      out <- append(vax_no_itn, control_no_itn)
    }

  } else if(outcome == "severe") {

    # Maternally acquired immunity
    iva20_no_itn <- get_maternal_start(eir = eir, gamma_llin = 1, season = season,
                                       rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db,
                                       b0 = p$b0, b1 = p$b1, ib0 = p$IB0, kb = p$kb,
                                       uc = p$uv, dc = p$dv, cpp = cpp)

    iva20_itn <- get_maternal_start(eir = eir, gamma_llin = gamma_llin, season = season,
                                    rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db,
                                    b0 = p$b0, b1 = p$b1, ib0 = p$IB0, kb = p$kb,
                                    uc = p$uv, dc = p$dv, cpp = cpp)

    ivm_no_itn <- get_icm(age = age, ica20 = iva20_no_itn, pm = p$PVM, dm = p$dvm)
    ivm_itn <- get_icm(age = age, ica20 = iva20_itn, pm = p$PVM, dm = p$dvm)

    # Simulate clinical hazard

    n_vaccine_arms <- length(vx)
    vax_no_itn <- list()
    vax_itn <- list()

    for (i in 1:n_vaccine_arms) {

      # For vaccine cohort, without bednets
      vax_no_itn[[i]] <- r_clin * get_severe_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                                      gamma_llin = 1, vx = vx[[i]], season = season,
                                                      ivm = ivm_no_itn, age = age, n = n, zeta = zeta, weight = weight,
                                                      rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                                                      b1 = p$b1, ib0 = p$IB0, kb = p$kb, uv = p$uv,
                                                      dv = p$dv, theta0 = p$theta0, theta1 = p$theta1, iv0 = p$IV0,
                                                      kv = p$kv, fv0 = p$fv0, av = p$av, gammav = p$gammav, cpp = cpp)

      # For vaccine cohort, with bednets
      if (gamma_llin != 1) {
        vax_itn[[i]] <- r_clin * get_severe_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                                     gamma_llin = gamma_llin, vx = vx[[i]], season = season,
                                                     ivm = ivm_itn, age = age, n = n, zeta = zeta, weight = weight,
                                                     rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                                                     b1 = p$b1, ib0 = p$IB0, kb = p$kb, uv = p$uv,
                                                     dv = p$dv, theta0 = p$theta0, theta1 = p$theta1, iv0 = p$IV0,
                                                     kv = p$kv, fv0 = p$fv0, av = p$av, gammav = p$gammav, cpp = cpp)
      }

    }

    names(vax_no_itn) <- paste0(rep("vax_no_itn_arm", n_vaccine_arms), 1:n_vaccine_arms)
    if (gamma_llin != 1) {
      names(vax_itn) <- paste0(rep("vax_itn_arm", n_vaccine_arms), 1:n_vaccine_arms)
    }

    # For unvaccinated cohort, without bednets
    control_no_itn <- r_clin * get_severe_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                                 gamma_llin = 1, vx = rep(0, length(age)),
                                                 season = season, ivm = ivm_no_itn, age = age,
                                                 n = n, zeta = zeta, weight = weight,
                                                 rho = p$rho, a0 = p$a0, ub = p$ub,
                                                 db = p$db, b0 = p$b0, b1 = p$b1,
                                                 ib0 = p$IB0, kb = p$kb, uv = p$uv, dv = p$dv,
                                                 theta0 = p$theta0, theta1 = p$theta1, iv0 = p$IV0,
                                                 kv = p$kv, fv0 = p$fv0, av = p$av, gammav = p$gammav, cpp = cpp)

    control_no_itn <- list(control_no_itn = control_no_itn)

    # For unvaccinated cohort, with bednets
    if (gamma_llin != 1) {
      control_itn <- r_clin * get_severe_hazard(eir = eir, age_at_enrollment = age_at_enrollment,
                                                gamma_llin = gamma_llin, vx = rep(0, length(age)),
                                                season = season, ivm = ivm_itn, age = age,
                                                n = n, zeta = zeta, weight = weight,
                                                rho = p$rho, a0 = p$a0, ub = p$ub,
                                                db = p$db, b0 = p$b0, b1 = p$b1,
                                                ib0 = p$IB0, kb = p$kb, uv = p$uv, dv = p$dv,
                                                theta0 = p$theta0, theta1 = p$theta1, iv0 = p$IV0,
                                                kv = p$kv, fv0 = p$fv0, av = p$av, gammav = p$gammav, cpp = cpp)

      control_itn <- list(control_itn = control_itn)
    }

    if(gamma_llin != 1) {
      out <- append(append(append(vax_itn, vax_no_itn), control_itn), control_no_itn)
    } else {
      out <- append(vax_no_itn, control_no_itn)
    }

  }



  return(out)

}




