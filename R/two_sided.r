two_sided_critvals <- function(theta0, r, n, m, alpha) {
  if (theta0 <= 0) stop("Error: Parameter 'theta0' must be positive")
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (n <= 0) stop("Error: Parameter 'n' must be a positive integer")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (alpha < 0 || alpha > 1) stop("Error: Parameter 'alpha' must be between 0 and 1")

  nu <- 2 * n * m / r

  equations <- function(vars, nu, theta0, r, alpha) {
    c1 <- vars[1]
    c2 <- vars[2]

    eq1 <- pchisq(2 * c2 * theta0 ^ r, nu) - pchisq(2 * c1 * theta0 ^ r, nu) - (1 - alpha)
    eq2 <- expint::gammainc(nu / 2 + 1, c1 * theta0 ^ r) - expint::gammainc(nu / 2 + 1, c2 * theta0 ^ r) - (1 - alpha) * (nu / 2) * gamma(nu / 2)

    return(c(eq1, eq2))
  }

  eps <- 0.01
  c1 <- 0
  max_iterations <- 1000  # Prevent infinite loops
  found_solution <- FALSE

  for (i in seq_len(max_iterations)) {
    c2 <- qchisq(1 - alpha + pchisq(2 * c1 * theta0 ^ r, nu), df = nu) / (2 * theta0 ^ r)
    
    # Try to find solution
    solution <- suppressWarnings(
      tryCatch({
        pracma::fsolve(equations, c(c1, c2), nu = nu, theta0 = theta0, r = r, alpha = alpha)
      }, error = function(e) NULL)
    )
    
    # Check if solution is valid
    if (!is.null(solution) && all(is.finite(solution$x))) {
      found_solution <- TRUE
      break
    }
    
    c1 <- c1 + eps
  }

  if (!found_solution) {
    stop("No solution found after ", max_iterations, " iterations. Consider adjusting eps or initial conditions.")
  }

  return(c(solution$x[1], solution$x[2]))
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
#' @returns A character vector
#' @examples
#' two_sided_rr(1, -1, 1, -1, 0.05)
#' two_sided_rr(1, -2, 4, -1, 0.05)
#' @export
two_sided_rr <- function(theta0, r, n, m, alpha) {
  c <- two_sided_critvals(theta0, r, n, m, alpha)
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
#' @returns A numeric value between 0 and 1
#' @examples
#' two_sided_beta(1, -1, 1, -1, 0.05)
#' two_sided_beta(1, -2, 4, -1, 0.05)
#' @export
two_sided_beta <- function(theta0, theta1, r, n, m, alpha) {
  if (theta1 <= 0) stop("Error: Parameter 'theta1' must be positive")
  c <- two_sided_critvals(theta0, r, n, m, alpha)
  1 - pchisq(q = 2 * c[2] * theta1 ^ r, df = 2 * n * m / r) + pchisq(q = 2 * c[1] * theta1 ^ r, df = 2 * n * m / r)
}