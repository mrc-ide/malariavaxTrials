#' Generate seasonality curve
#'
#' @param t Time (days). Default = 1:365.
#' @param g0 Fourier series parameter.
#' @param g1 Fourier series parameter.
#' @param g2 Fourier series parameter.
#' @param g3 Fourier series parameter.
#' @param h1 Fourier series parameter.
#' @param h2 Fourier series parameter.
#' @param h3 Fourier series parameter.
#' @param floor Minimum value in seasonality curve.
#' @export
get_season <- function(t = 1:365, g0, g1, g2, g3, h1, h2, h3, floor){

  t <- t / 365
  prediction <- data.frame(
    g0 = g0,
    g1 = g1 * cos(2 * pi * t * 1),
    g2 = g2 * cos(2 * pi * t * 2),
    g3 = g3 * cos(2 * pi * t * 3),
    h1 = h1 * sin(2 * pi * t * 1),
    h2 = h2 * sin(2 * pi * t * 2),
    h3 = h3 * sin(2 * pi * t * 3))
  prediction <- rowSums(prediction)
  prediction <- pmax(floor, prediction)

}
