#' Critical Value for Simple MP and One-Sided UMP Tests
#' 
#' `simple_critval` computes the critical value for the MP test of
#' H₀: θ = θ₀ vs H₁: θ = θ₁, or UMP tests of
#' H₀: θ ≤ θ₀ vs H₁: θ > θ₀ and H₀: θ ≥ θ₀ vs H₁: θ < θ₀
#' 
#' @param theta0 A positive numeric value representing θ₀
#' @param theta1 A positive numeric value representing θ₁
#' @param r A non-zero numeric parameter of the scale-exponential family
#' @param n An integer representing the sample size
#' @param m A non-zero numeric parameter of the scale-exponential family
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A numeric critical value
#' @examples
#' simple_critval(2, 4, -2, 20, -1, 0.01)
#' simple_critval(1 / 8, 1 / 10, -1, 9, -1, 0.05)
#' simple_critval(2, 1, 1, 1, 1, 0.05)
#' simple_critval(4, 6, 1, 1, 1, 0.05)
#' simple_critval(4, 2, 1, 1, 1, 0.05)
#' simple_critval(1.5, 3, -1, 10, -1, 0.05)
#' simple_critval(10, 14, -1, 31, -1, 0.05)
#' simple_critval(2, 1.25, -2, 40, -1, 0.025)
#' simple_critval(4, 12, -2, 25, -1, 0.05)
#' @export
simple_critval <- function(theta0, theta1, r, n, m, alpha) {
  if (theta0 <= 0) stop("Error: Parameter 'theta0' must be positive")
  if (theta1 <= 0) stop("Error: Parameter 'theta1' must be positive")
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (n <= 0) stop("Error: Parameter 'n' must be a positive integer")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (alpha < 0 || alpha > 1) stop("Error: Parameter 'alpha' must be between 0 and 1")

  if (theta1 ^ r > theta0 ^ r)
    qchisq(p = alpha, df = 2 * n * m / r)
  else
    qchisq(p = 1 - alpha, df = 2 * n * m / r)
}

#' Rejection Region for Simple MP and One-Sided UMP Tests
#' 
#' `simple_rr` computes the rejection region for the MP test of
#' H₀: θ = θ₀ vs H₁: θ = θ₁, or UMP tests of
#' H₀: θ ≤ θ₀ vs H₁: θ > θ₀ and H₀: θ ≥ θ₀ vs H₁: θ < θ₀
#'  
#' @param theta0 A positive numeric value representing θ₀
#' @param theta1 A positive numeric value representing θ₁
#' @param r A non-zero numeric parameter of the scale-exponential family
#' @param n An integer representing the sample size
#' @param m A non-zero numeric parameter of the scale-exponential family
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A character string describing the rejection region in interval notation
#' @examples
#' simple_rr(2, 4, -2, 20, -1, 0.01)
#' simple_rr(1 / 8, 1 / 10, -1, 9, -1, 0.05)
#' simple_rr(2, 1, 1, 1, 1, 0.05)
#' simple_rr(4, 6, 1, 1, 1, 0.05)
#' simple_rr(4, 2, 1, 1, 1, 0.05)
#' simple_rr(1.5, 3, -1, 10, -1, 0.05)
#' simple_rr(10, 14, -1, 31, -1, 0.05)
#' simple_rr(2, 1.25, -2, 40, -1, 0.025)
#' simple_rr(4, 12, -2, 25, -1, 0.05)
#' @export
simple_rr <- function(theta0, theta1, r, n, m, alpha) {
  c <- simple_critval(theta0, theta1, r, n, m, alpha)
  if (theta1 ^ r > theta0 ^ r)
    sprintf("(%f, %4f)", -Inf, c)
  else
    sprintf("(%4f, %f)", c, Inf)
}

#' Power for Simple MP and One-Sided UMP Tests
#' 
#' `simple_beta` computes statistical power for the MP test of
#' H₀: θ = θ₀ vs H₁: θ = θ₁, or UMP tests of
#' H₀: θ ≤ θ₀ vs H₁: θ > θ₀ and H₀: θ ≥ θ₀ vs H₁: θ < θ₀,
#' at θ = θ₁.
#' 
#' @param theta0 A positive numeric value representing θ₀
#' @param theta1 A positive numeric value representing θ₁
#' @param r A non-zero numeric parameter of the scale-exponential family
#' @param n An integer representing the sample size
#' @param m A non-zero numeric parameter of the scale-exponential family
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A numeric value between 0 and 1 representing power at θ = θ₁
#' @examples
#' simple_beta(2, 4, -2, 20, -1, 0.01)
#' simple_beta(1 / 8, 1 / 10, -1, 9, -1, 0.05)
#' simple_beta(2, 1, 1, 1, 1, 0.05)
#' simple_beta(4, 6, 1, 1, 1, 0.05)
#' simple_beta(4, 2, 1, 1, 1, 0.05)
#' simple_beta(1.5, 3, -1, 10, -1, 0.05)
#' simple_beta(10, 14, -1, 31, -1, 0.05)
#' simple_beta(2, 1.25, -2, 40, -1, 0.025)
#' simple_beta(2, 12, -2, 25, -1, 0.05)
#' @export
simple_beta <- function(theta0, theta1, r, n, m, alpha) {
  c <- simple_critval(theta0, theta1, r, n, m, alpha)
  if (theta1 ^ r > theta0 ^ r)
    pchisq(q = theta1 ^ r / theta0 ^ r * c, df = 2 * n * m / r)
  else
    1 - pchisq(q = theta1 ^ r / theta0 ^ r * c, df = 2 * n * m / r)
}