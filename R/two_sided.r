two_sided_critvals <- function(theta0, r, n, m, alpha, initial_guess) {
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (theta0 <= 0) stop("Error: Parameter 'theta0' must be positive")

  df <- n * m / r
  equations <- function(vars, nu, theta0, r, alpha) {
    c1 <- vars[1]
    c2 <- vars[2]

    eq1 <- expint::gammainc(df, c1 * theta0 ^ r) - expint::gammainc(df, c2 * theta0 ^ r) - (1 - alpha) * gamma(df)
    eq2 <- expint::gammainc(df + 1, c1 * theta0 ^ r) - expint::gammainc(df + 1, c2 * theta0 ^ r) - (1 - alpha) * df * gamma(df)

    return(c(eq1, eq2))
  }
  solution <- pracma::fsolve(equations, initial_guess, nu = nu, theta0 = theta0, r = r, alpha = alpha)

  c1 <- solution$x[1]
  c2 <- solution$x[2]

  return(c(c1, c2))
}

#' Rejection Region for Two-Sided Hypothesis Tests
#' 
#' `two_sided_rr()` computes the rejection region for simple vs two-sided hypotheses.
#' 
#' @param theta0 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @param initial_guess Numeric vector of length 2
#' @returns A character vector
#' @examples
#' two_sided_rr(1, -1, 1, -1, 0.05, c(0.05, 5))
#' two_sided_rr(1, -2, 4, -1, 0.05, c(0.3, 6))
#' @export
two_sided_rr <- function(theta0, r, n, m, alpha, initial_guess) {
  c <- two_sided_critvals(theta0, r, n, m, alpha, initial_guess)
  sprintf("(%f, %4f] U [%4f, %f)", -Inf, c[1], c[2], Inf)
}

#' Power for Two-Sided Hypothesis Tests
#' 
#' `two_sided_beta()` computes statistical power at $\theta = \theta_1$ for simple vs two-sided hypotheses.
#' 
#' @param theta0 A positive numeric value
#' @param theta1 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @param initial_guess Numeric vector of length 2
#' @returns A numeric value between 0 and 1
#' @examples
#' two_sided_beta(1, -1, 1, -1, 0.05, c(0.05, 5))
#' two_sided_beta(1, -2, 4, -1, 0.05, c(0.3, 6))
#' @export
two_sided_beta <- function(theta0, theta1, r, n, m, alpha, initial_guess) {
  c <- two_sided_critvals(theta0, r, n, m, alpha, initial_guess)
  if (theta1 <= 0)
    stop("Error: Parameter 'theta1' must be positive")
  1 - pchisq(q = 2 * c[2] * theta1 ^ r, df = 2 * n * m / r) + pchisq(q = 2 * c[1] * theta1 ^ r, df = 2 * n * m / r)
}