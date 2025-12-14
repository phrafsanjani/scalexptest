#' Critical Value for Two-Sided UMPU Test
#' 
#' `two_sided_critvals` computes the critical values for the UMPU test of
#' H₀: θ = θ₀ vs H₁: θ ≠ θ₀.
#' 
#' @param theta0 A positive numeric value representing θ₀
#' @param r A non-zero numeric parameter of the scale-exponential family
#' @param n An integer representing the sample size
#' @param m A non-zero numeric parameter of the scale-exponential family
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A numeric length-two vector of critical values
#' @examples
#' two_sided_critvals(1, -1, 1, -1, 0.05)
#' two_sided_critvals(1, -2, 4, -1, 0.05)
#' two_sided_critvals(5, -1, 500, -1, 0.05)
#' @export
two_sided_critvals <- function(theta0, r, n, m, alpha) {
  if (theta0 <= 0) stop("Error: Parameter 'theta0' must be positive")
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (n <= 0) stop("Error: Parameter 'n' must be a positive integer")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (alpha < 0 || alpha > 1) stop("Error: Parameter 'alpha' must be between 0 and 1")

  nu <- 2 * n * m / r

  equations <- function(vars, nu, alpha) {
    c1 <- vars[1]
    c2 <- vars[2]

    eq1 <- pchisq(c2, nu) - pchisq(c1, nu) - (1 - alpha)
    eq2 <- pchisq(c2, nu + 2) - pchisq(c1, nu + 2) - (1 - alpha)

    return(c(eq1, eq2))
  }

  solution <- nleqslv::nleqslv(
    x  = c(qchisq(alpha / 2, df = nu), qchisq(1 - alpha/2, df = nu)),
    fn = equations,
    nu = nu,
    alpha = alpha
  )

  return(c(solution$x[1], solution$x[2]))
}

#' Rejection Region for Two-Sided UMPU Test
#' 
#' `two_sided_rr` computes the rejection region for the UMPU test of
#' H₀: θ = θ₀ vs H₁: θ ≠ θ₀.
#' 
#' @param theta0 A positive numeric value representing θ₀
#' @param r A non-zero numeric parameter of the scale-exponential family
#' @param n An integer representing the sample size
#' @param m A non-zero numeric parameter of the scale-exponential family
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A character string describing the rejection region in interval notation
#' @examples
#' two_sided_rr(1, -1, 1, -1, 0.05)
#' two_sided_rr(1, -2, 4, -1, 0.05)
#' two_sided_rr(5, -1, 500, -1, 0.05)
#' @export
two_sided_rr <- function(theta0, r, n, m, alpha) {
  c <- two_sided_critvals(theta0, r, n, m, alpha)
  sprintf("(%f, %4f) U (%4f, %f)", -Inf, c[1], c[2], Inf)
}

#' Power for Two-Sided UMPU Test
#' 
#' `two_sided_beta` computes statistical power for the UMPU test of
#' H₀: θ = θ₀ vs H₁: θ ≠ θ₀ at `theta`.
#' 
#' @param theta0 A positive numeric value representing θ₀
#' @param theta A positive numeric value
#' @param r A non-zero numeric parameter of the scale-exponential family
#' @param n An integer representing the sample size
#' @param m A non-zero numeric parameter of the scale-exponential family
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A numeric value between 0 and 1 representing power at `theta`
#' @examples
#' two_sided_beta(1, -1, 1, -1, 0.05)
#' two_sided_beta(1, -2, 4, -1, 0.05)
#' two_sided_beta(5, 5.75, -1, 500, -1, 0.05)
#' @export
two_sided_beta <- function(theta0, theta, r, n, m, alpha) {
  if (theta <= 0) stop("Error: Parameter 'theta' must be positive")
  c <- two_sided_critvals(theta0, r, n, m, alpha)
  1 - pchisq(theta ^ r / theta0 ^ r * c[2], df = 2 * n * m / r) + pchisq(theta ^ r / theta0 ^ r * c[1], df = 2 * n * m / r)
}