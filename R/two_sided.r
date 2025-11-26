two_sided_critvals <- function(theta0, r, n, m, alpha, a = 0, eps = 0.01, max_iterations = 5000) {
  if (theta0 <= 0) stop("Error: Parameter 'theta0' must be positive")
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (n <= 0) stop("Error: Parameter 'n' must be a positive integer")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (alpha < 0 || alpha > 1) stop("Error: Parameter 'alpha' must be between 0 and 1")

  nu <- 2 * n * m / r

  equations <- function(vars, nu, theta0, r, alpha) {
    c1 <- vars[1]
    c2 <- vars[2]

    eq1 <- pchisq(c2, nu) - pchisq(c1, nu) - (1 - alpha)
    eq2 <- pchisq(c2, nu + 2) - pchisq(c1, nu + 2) - (1 - alpha)

    return(c(eq1, eq2))
  }

  found_solution <- FALSE

  for (i in seq_len(max_iterations)) {
    b <- qchisq(1 - alpha + pchisq(a, nu), df = nu)
    
    solution <- suppressWarnings(
      tryCatch({
        pracma::fsolve(equations, c(a, b), nu = nu, theta0 = theta0, r = r, alpha = alpha)
      }, error = function(e) NULL)
    )
    
    if (!is.null(solution) && all(is.finite(solution$x)) && solution$x[1] < solution$x[2]) {
      found_solution <- TRUE
      break
    }
    
    a <- a + eps
  }

  if (!found_solution)
    stop("No solution found after ", max_iterations, " iterations. Consider adjusting 'eps' or 'max_iterations'.")

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
#' two_sided_rr(5, -1, 500, -1, 0.05, a = 2000, eps = 1)
#' @export
two_sided_rr <- function(theta0, r, n, m, alpha, a = 0, eps = 0.01, max_iterations = 5000) {
  c <- two_sided_critvals(theta0, r, n, m, alpha, a, eps, max_iterations)
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
#' two_sided_beta(5, 5.75, -1, 500, -1, 0.05, a = 2000, eps = 1)
#' @export
two_sided_beta <- function(theta0, theta, r, n, m, alpha, a = 0, eps = 0.01, max_iterations = 5000) {
  if (theta <= 0) stop("Error: Parameter 'theta' must be positive")
  c <- two_sided_critvals(theta0, r, n, m, alpha, a, eps, max_iterations)
  1 - pchisq(theta ^ r / theta0 ^ r * c[2], df = 2 * n * m / r) + pchisq(theta ^ r / theta0 ^ r * c[1], df = 2 * n * m / r)
}