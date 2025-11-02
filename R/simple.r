simple_critval <- function(theta0, theta1, r, n, m, alpha) {
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (theta0 <= 0) stop("Error: Parameter 'theta0' must be positive")
  if (theta1 <= 0) stop("Error: Parameter 'theta1' must be positive")

  if (theta1 ^ r > theta0 ^ r)
    qchisq(p = alpha, df = 2 * n * m / r) / (2 * theta0 ^ r)
  else
    qchisq(p = 1 - alpha, df = 2 * n * m / r) / (2 * theta0 ^ r)
}

#' Rejection Region for Simple and One-Sided Hypotheses
#' 
#' `simple_rr()` computes the rejection region for simple vs simple,
#' simple vs one-sided and one-sided vs one-sided hypotheses.
#' 
#' @param theta0 A positive numeric value
#' @param theta1 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A character vector
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

#' Power for Simple and One-Sided Hypotheses
#' 
#' `simple_beta()` computes statistical power at $\theta = \theta_1$ for simple vs simple,
#' simple vs one-sided and one-sided vs one-sided tests.
#' 
#' @param theta0 A positive numeric value
#' @param theta1 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A numeric value between 0 and 1
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
  if (theta1 ^ r > theta0 ^ r)
    pchisq(q = (theta1 / theta0) ^ r * qchisq(p = alpha, df = 2 * n * m / r), df = 2 * n * m / r)
  else
    1 - pchisq(q = (theta1 / theta0) ^ r * qchisq(p = 1 - alpha, df = 2 * n * m / r), df = 2 * n * m / r)
}