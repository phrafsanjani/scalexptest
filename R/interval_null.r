#' Critical Values for Interval-null UMPU Test
#' 
#' `interval_null_critvals` computes the critical values for the UMPU test of
#' H₀: θ₁ ≤ θ ≤ θ₂ vs H₁: θ < θ₁ or θ > θ₂.
#' 
#' @param theta1 A positive numeric value representing θ₁
#' @param theta2 A positive numeric value representing θ₂ (must be greater than θ₁)
#' @param r A non-zero numeric parameter of the scale-exponential family
#' @param n An integer representing the sample size
#' @param m A non-zero numeric parameter of the scale-exponential family
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @return A numeric length-two vector of critical values
#' @examples
#' interval_null_critvals(1.75, 2.25, -2, 200, -1, 0.05)
#' @export
interval_null_critvals <- function(theta1, theta2, r, n, m, alpha) {
  if (theta1 <= 0) stop("Error: Parameter 'theta1' must be positive")
  if (theta2 <= 0) stop("Error: Parameter 'theta2' must be positive")
  if (theta2 <= theta1) stop("theta2 must be greater than theta1")
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (n <= 0) stop("Error: Parameter 'n' must be a positive integer")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (alpha < 0 || alpha > 1) stop("Error: Parameter 'alpha' must be between 0 and 1")

  nu <- 2 * n * m / r

  a_seq <- seq(1e-5, qchisq(alpha - 1e-5, df = nu), length.out = 1000)

  solution_found <- FALSE

  for (a in a_seq) {
    b <- qchisq(pchisq(a, df = nu) + 1 - alpha, df = nu)
    equations <- function(vars, nu, theta1, theta2, r, alpha) {
      c1 <- vars[1]
      c2 <- vars[2]

      eq1 <- pchisq(c2, nu) - pchisq(c1, nu) - (1 - alpha)
      eq2 <- pchisq(theta2 ^ r / theta1 ^ r * c2, nu) - pchisq(theta2 ^ r / theta1 ^ r * c1, nu) - (1 - alpha)

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

    if (!is.null(solution) && solution$termcd %in% c(1, 2) && c1 < c2 && pchisq(c2, nu) - pchisq(c1, nu) - (1 - alpha) >= 0 && pchisq(theta2 ^ r / theta1 ^ r * c2, nu) - pchisq(theta2 ^ r / theta1 ^ r * c1, nu) - (1 - alpha) >= 0) {
      solution_found <- TRUE
      break
    }
  }

  if (!solution_found) {
    a <- qchisq(alpha - 1e-5, df = nu)
    b <- qchisq(theta2 ^ r / theta1 ^ r * pchisq(a, df = nu) + 1 - alpha, df = nu) * theta1 ^ r / theta2 ^ r
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

#' Rejection Region for Interval-null UMPU Test
#' 
#' `interval_null_rr` computes the rejection region for the UMPU test of
#' H₀: θ₁ ≤ θ ≤ θ₂ vs H₁: θ < θ₁ or θ > θ₂.
#' 
#' @param theta1 A positive numeric value representing θ₁
#' @param theta2 A positive numeric value representing θ₂
#' @param r A non-zero numeric parameter of the scale-exponential family
#' @param n An integer representing the sample size
#' @param m A non-zero numeric parameter of the scale-exponential family
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @return A character string describing the rejection region in interval notation
#' @examples
#' interval_null_rr(1.75, 2.25, -2, 200, -1, 0.05)
#' @export
interval_null_rr <- function(theta1, theta2, r, n, m, alpha) {
  c <- interval_null_critvals(theta1, theta2, r, n, m, alpha)
  sprintf("(%f, %4f) U (%4f, %f)", -Inf, c[1], c[2], Inf)
}

#' Power for Interval-null UMPU Test
#' 
#' `interval_null_beta` computes statistical power for the UMPU test of
#' H₀: θ₁ ≤ θ ≤ θ₂ vs H₁: θ < θ₁ or θ > θ₂ at θ.
#' 
#' @param theta1 A positive numeric value representing θ₁
#' @param theta2 A positive numeric value representing θ₂
#' @param theta A numeric value representing θ
#' @param r A non-zero numeric parameter of the scale-exponential family
#' @param n An integer representing the sample size
#' @param m A non-zero numeric parameter of the scale-exponential family
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @return A numeric value between 0 and 1 representing the power at θ
#' @examples
#' interval_null_beta(1.75, 2.25, 2.29, -2, 200, -1, 0.05)
#' @export
interval_null_beta <- function(theta1, theta2, theta, r, n, m, alpha) {
  c <- interval_null_critvals(theta1, theta2, r, n, m, alpha)
  1 - pchisq(theta ^ r / theta1 ^ r * c[2], df = 2 * n * m / r) + pchisq(theta ^ r / theta1 ^ r * c[1], df = 2 * n * m / r)
}