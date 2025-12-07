#' Rejection Region for Interval Hypothesis Tests
#' 
#' `interval_alt_critvals()` computes the critival values for $H_0 : \theta \leq \theta_1 \text{ or } \theta \geq \theta_2$ vs $H_1 : \theta_1 < \theta < \theta_2$.
#' 
#' @param theta1 A positive numeric value
#' @param theta2 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A numeric vector
#' @examples
#' interval_alt_critvals(1, sqrt(3), -2, 25, -1, 0.01)
#' @export
interval_alt_critvals <- function(theta1, theta2, r, n, m, alpha) {
  if (theta1 <= 0) stop("Error: Parameter 'theta1' must be positive")
  if (theta2 <= 0) stop("Error: Parameter 'theta2' must be positive")
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (n <= 0) stop("Error: Parameter 'n' must be a positive integer")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (alpha < 0 || alpha > 1) stop("Error: Parameter 'alpha' must be between 0 and 1")

  nu <- 2 * n * m / r

  a_seq <- seq(1e-5, qchisq(1 - alpha - 1e-5, df = nu), length.out = 1000)

  solution_found <- FALSE

  for (a in a_seq) {
    b <- qchisq(pchisq(a, df = nu) + alpha, df = nu)
    equations <- function(vars, nu, theta1, theta2, r, alpha) {
      c1 <- vars[1]
      c2 <- vars[2]

      eq1 <- pchisq(c2, nu) - pchisq(c1, nu) - alpha
      eq2 <- pchisq(theta2 ^ r / theta1 ^ r * c2, nu) - pchisq(theta2 ^ r / theta1 ^ r * c1, nu) - alpha

      return(c(eq1, eq2))
    }

    solution <- nleqslv::nleqslv(
      x  = c(a, b),
      fn = equations,
      nu = nu,
      theta1 = theta1,
      theta2 = theta2,
      r = r,
      alpha = alpha
    )

    c1 <- solution$x[1]
    c2 <- solution$x[2]

    if (!is.null(solution) && solution$termcd %in% c(1, 2) && c1 < c2 && pchisq(c2, nu) - pchisq(c1, nu) - alpha >= 0 && pchisq(theta2 ^ r / theta1 ^ r * c2, nu) - pchisq(theta2 ^ r / theta1 ^ r * c1, nu) - alpha >= 0) {
      solution_found <- TRUE
      break
    }
  }

  if (!solution_found) {
    a <- qchisq(1 - alpha - 1e-5, df = nu)
    b <- qchisq(theta2 ^ r / theta1 ^ r * pchisq(a, df = nu) + alpha, df = nu) * theta1 ^ r / theta2 ^ r
    solution <- nleqslv::nleqslv(
      x  = c(a, b),
      fn = equations,
      nu = nu,
      theta1 = theta1,
      theta2 = theta2,
      r = r,
      alpha = alpha
    )

    c1 <- solution$x[1]
    c2 <- solution$x[2]
  }

  return(c(c1, c2))
}

#' Rejection Region for Interval Hypothesis Tests
#' 
#' `interval_alt_rr()` computes the rejection region for $H_0 : \theta \leq \theta_1 \text{ or } \theta \geq \theta_2$ vs $H_1 : \theta_1 < \theta < \theta_2$.
#' 
#' @param theta1 A positive numeric value
#' @param theta2 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A character vector
#' @examples
#' interval_alt_rr(1, sqrt(3), -2, 25, -1, 0.01)
#' @export
interval_alt_rr <- function(theta1, theta2, r, n, m, alpha) {
  c <- interval_alt_critvals(theta1, theta2, r, n, m, alpha)
  if (r > 0)
    sprintf("(%f, %4f) U (%4f, %f)", -Inf, c[1], c[2], Inf)
  else
    sprintf("(%4f, %4f)", c[1], c[2])
}

#' Power for Interval Hypothesis Tests
#' 
#' `interval_alt_beta()` computes statistical power at $\theta$ for interval hypothesis tests.
#' 
#' @param theta1 A positive numeric value
#' @param theta2 A positive numeric value
#' @param theta A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A numeric value between 0 and 1
#' @examples
#' interval_alt_beta(1, sqrt(3), sqrt(2), -2, 25, -1, 0.01)
#' @export
interval_alt_beta <- function(theta1, theta2, theta, r, n, m, alpha) {
  c <- interval_alt_critvals(theta1, theta2, r, n, m, alpha)
  if (r > 0)
    1 - pchisq(theta ^ r / theta1 ^ r * c[2], df = 2 * n * m / r) + pchisq(theta ^ r / theta1 ^ r * c[1], df = 2 * n * m / r)
  else
    pchisq(theta ^ r / theta1 ^ r * c[2], df = 2 * n * m / r) - pchisq(theta ^ r / theta1 ^ r * c[1], df = 2 * n * m / r)
}