test_that("simulate_trial_hazards works", {
  library(boostr)

  # Time of doses (assuming first is dose 3 of primary series)
  td <- c(1, 365, 365 * 2)
  # Peaks following primary series or boost
  init_titres <- c(100, 80, 90)
  # Proportion short-lived
  prop_short <- c(1, 0.5, 0.6)
  # Decay short
  dur_short <- c(60, 60, 60)
  # Decay long
  dur_long <- c(300, 300, 300)

  # One vaccine cohort
  ab <- ab(timesteps=10*365, dose_timesteps = td[1], init_titres = init_titres[1],
           prop_short = prop_short[1], dur_short = dur_short[1], dur_long = dur_long[1])

  vx <- efficacy(titre=ab, max_efficacy = 0.9, alpha = 0.6, beta = 100)

  out_cpp <- simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(vx),
                                 n = 10, r_clin = 1.1, age = 1:(10*365), cpp = TRUE)

  out_r <- simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(vx),
                                    n = 10, r_clin = 1.1, age = 1:(10*365), cpp = FALSE)

  expect_identical(unlist(out_cpp), unlist(out_r))
  expect_true(all(sapply(out_cpp, function(x) all(x > 0)))==TRUE)
  expect_length(out_cpp, length(list(vx))*2+2)
  expect_true(all(1-out_cpp$vax_itn_arm1/out_cpp$control_itn <= 1))
  expect_true(all(1-out_cpp$vax_no_itn_arm1/out_cpp$control_no_itn <= 1))

  out_cpp <- simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 1, vx = list(vx),
                                    n = 10, r_clin = 1.1, age = 1:(10*365), cpp = TRUE)
  expect_length(out_cpp, length(list(vx))+1)

  # Two vaccine cohorts
  ab2 <- ab(timesteps=10*365, dose_timesteps = td, init_titres = init_titres,
           prop_short = prop_short, dur_short = dur_short, dur_long = dur_long)

  vx2 <- efficacy(titre=ab2, max_efficacy = 0.9, alpha = 0.6, beta = 100)

  out_cpp <- simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(vx, vx2),
                                    n = 10, r_clin = 1.1, age = 1:(10*365), cpp = TRUE)

  out_r <- simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(vx, vx2),
                                  n = 10, r_clin = 1.1, age = 1:(10*365), cpp = FALSE)

  expect_identical(unlist(out_cpp), unlist(out_r))
  expect_true(all(sapply(out_cpp, function(x) all(x > 0)))==TRUE)
  expect_length(out_cpp, length(list(vx, vx2))*2+2)
  expect_true(all(1-out_cpp$vax_itn_arm2/out_cpp$control_itn <= 1))
  expect_true(all(1-out_cpp$vax_no_itn_arm2/out_cpp$control_no_itn <= 1))

  # No vaccine dose in vaccine arm
  # now 1 and 3 should be identical
  out_cpp <- simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(rep(0,10*365)),
                                    n = 10, r_clin = 1.1, age = 1:(10*365), cpp = TRUE)

  out_r <- simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(rep(0,10*365)),
                                    n = 10, r_clin = 1.1, age = 1:(10*365), cpp = FALSE)

  expect_identical(unlist(out_cpp), unlist(out_r))
  expect_true(all(sapply(out_cpp, function(x) all(x > 0)))==TRUE)
  expect_true(all.equal(out_cpp$vax_itn_arm1, out_cpp$control_itn))
  expect_true(all.equal(out_cpp$vax_no_itn_arm1, out_cpp$control_no_itn))

  # Test error messages
  expect_error(
    simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = vx,
                           n = 10, r_clin = 1.1, age = 1:(10*365), cpp = TRUE),
    "vx must be a list and contain at least 1 vector"
  )

  expect_error(
    simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(vx[1:(5*365)]),
                           n = 10, r_clin = 1.1, age = 1:(10*365), cpp = TRUE),
    "vx vectors must be the same length as the age vector or longer"
  )


  expect_error(
    simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(vx),
                           n = 10, r_clin = 1.1, age = 1:300, cpp = TRUE),
    "age_at_enrollment must be within age vector"
  )

  expect_error(
    simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(vx),
                           n = 10, season = rep(1, 400), r_clin = 1.1, age = 1:(10*365), cpp = TRUE),
    "The length of the season vector must be a multiple of 365"
  )

  seas <- get_season(t= 1:365, g0 = 0.6774, g1 = -0.2493, g2 = -0.1334, g3 = 0.1590, h1 = 0.8928, h2 = 0.1165,
                     h3 = 0.0188, floor = -0.1)
  expect_error(
    simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(vx),
                           n = 10, season = seas, r_clin = 1.1, age = 1:(10*365), cpp = TRUE),
    "All values in season must be positive"
  )

  expect_error(
    simulate_trial_hazards(eir = 20/365, age_at_enrollment = 12*30, gamma_llin = 0.8, vx = list(vx),
                           n = 0, r_clin = 1.1, age = 1:(10*365), cpp = TRUE),
    "n must be an integer >= 1"
  )

})

test_that("get_clinical_hazard works", {
  library(boostr)

  # Time of doses (assuming first is dose 3 of primary series)
  td <- c(1, 365, 365 * 2)
  # Peaks following primary series or boost
  init_titres <- c(100, 80, 90)
  # Proportion short-lived
  prop_short <- c(1, 0.5, 0.6)
  # Decay short
  dur_short <- c(60, 60, 60)
  # Decay long
  dur_long <- c(300, 300, 300)

  # One vaccine cohort
  ab <- ab(timesteps=10*365, dose_timesteps = td[1], init_titres = init_titres[1],
           prop_short = prop_short[1], dur_short = dur_short[1], dur_long = dur_long[1])

  vx <- efficacy(titre=ab, max_efficacy = 0.9, alpha = 0.6, beta = 100)

  # Load fitted model parameters
  name_full <- system.file("extdata/", "Jamie_parameters.rds", package = 'malariavaxTrials', mustWork = TRUE)
  p <- readRDS(name_full)

  # Heterogeneity
  n <- 10
  gh <- statmod::gauss.quad.prob(n = n, dist = "normal")
  zeta <- exp(-p$s2 * 0.5 + sqrt(p$s2) * gh$nodes)
  weight <- gh$weight

  # Maternal immunity
  ica20 <- get_maternal_start(eir = 20/365, gamma_llin = 1, season = rep(1,365),
                              rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db,
                              b0 = p$b0, b1 = p$b1, ib0 = p$IB0, kb = p$kb,
                              uc = p$uc, dc = p$dc, cpp = TRUE)

  icm <- get_icm(age = 1:(10*365), ica20 = ica20, pm = p$PM, dm =p$dm)

  ch <- get_clinical_hazard(eir= 20/365, age_at_enrollment = 12*30, gamma_llin= 1, vx = vx,
                          season=rep(1,365), icm = icm, age = 1:(10*365), n = n, zeta = zeta, weight = weight,
                          rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                          b1 = p$b1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                          dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                          kc = p$kc, cpp = TRUE)

  expect_true(all(ch > 0))


  expect_error(
    get_clinical_hazard(eir= 20/365, age_at_enrollment = 12*30, gamma_llin= 1, vx = vx,
                        season=rep(1,365), icm = icm, age = 1:(12*365), n = n, zeta = zeta, weight = weight,
                        rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                        b1 = p$b1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                        dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                        kc = p$kc, cpp = TRUE),
    "vx must be the same length as the age vector or longer"
  )

  expect_error(
    get_clinical_hazard(eir= 20/365, age_at_enrollment = 12*30, gamma_llin= 1, vx = vx,
                        season=rep(1,365), icm = icm, age = 1:300, n = n, zeta = zeta, weight = weight,
                        rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                        b1 = p$b1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                        dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                        kc = p$kc, cpp = TRUE),
    "age_at_enrollment must be within age vector"
  )


  expect_error(
    get_clinical_hazard(eir= 20/365, age_at_enrollment = 12*30, gamma_llin= 1, vx = vx,
                        season=rep(1,300), icm = icm, age = 1:(10*365), n = n, zeta = zeta, weight = weight,
                        rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                        b1 = p$b1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                        dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                        kc = p$kc, cpp = TRUE),
    "The length of the season vector must be a multiple of 365"
  )

  expect_error(
    get_clinical_hazard(eir= 20/365, age_at_enrollment = 12*30, gamma_llin= 1, vx = vx,
                        season=rep(1,365), icm = icm, age = 1:(10*365), n = 0, zeta = zeta, weight = weight,
                        rho = p$rho, a0 = p$a0, ub = p$ub, db = p$db, b0 = p$b0,
                        b1 = p$b1, ib0 = p$IB0, kb = p$kb, uc = p$uc,
                        dc = p$dc, phi0 = p$phi0, phi1 = p$phi1, ic0 = p$IC0,
                        kc = p$kc, cpp = TRUE),
    "n must be an integer >= 1"
  )

})
