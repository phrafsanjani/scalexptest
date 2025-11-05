two_sided_critvals <- function(theta0, r, n, m, alpha, initial_guess) {
  if (theta0 <= 0) stop("Error: Parameter 'theta0' must be positive")
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (n <= 0) stop("Error: Parameter 'n' must be a positive integer")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (alpha < 0 || alpha > 1) stop("Error: Parameter 'alpha' must be between 0 and 1")
  if (!is.numeric(initial_guess) || length(initial_guess) != 2)
    stop("Error: Parameter 'initial_guess' must be a numeric two-element vector")
  if (initial_guess[1] >= initial_guess[2])
    stop("Error: initial_guess[1] must be less than initial_guess[2]")

  nu <- 2 * n * m / r
  equations <- function(vars, nu, theta0, r, alpha) {
    c1 <- vars[1]
    c2 <- vars[2]

    eq1 <- expint::gammainc(nu / 2, c1 / 2) - expint::gammainc(nu / 2, c2 / 2) - (1 - alpha) * gamma(nu / 2)
    eq2 <- expint::gammainc(nu / 2 + 1, c1 / 2) - expint::gammainc(nu / 2 + 1, c2 / 2) - (1 - alpha) * (nu / 2) * gamma(nu / 2)

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
  if (theta1 <= 0) stop("Error: Parameter 'theta1' must be positive")
  c <- two_sided_critvals(theta0, r, n, m, alpha, initial_guess)
  1 - pchisq(q = (theta1 / theta0) ^ r * c[2], df = 2 * n * m / r) + pchisq(q = (theta1 / theta0) ^ r * c[1], df = 2 * n * m / r)
}