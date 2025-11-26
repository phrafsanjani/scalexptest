#' Rejection Region for Interval Hypothesis Tests
#' 
#' `interval_null_critvals()` computes the critival values for $H_0 : \theta_1 < \theta < \theta_2$ vs $H_1 : \theta \leq \theta_1 \text{ or } \theta \geq \theta_2$.
#' 
#' @param theta1 A positive numeric value
#' @param theta2 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A numeric vector
#' @examples
#' interval_null_critvals(1.75, 2.25, -2, 200, -1, 0.05)
#' @export
interval_null_critvals <- function(theta1, theta2, r, n, m, alpha, a = 0, eps = 0.01, max_iterations = 5000) {
  if (theta1 <= 0) stop("Error: Parameter 'theta1' must be positive")
  if (theta2 <= 0) stop("Error: Parameter 'theta2' must be positive")
  if (r == 0) stop("Error: Parameter 'r' cannot be 0")
  if (n <= 0) stop("Error: Parameter 'n' must be a positive integer")
  if (m == 0) stop("Error: Parameter 'm' cannot be 0")
  if (alpha < 0 || alpha > 1) stop("Error: Parameter 'alpha' must be between 0 and 1")

  nu <- 2 * n * m / r

  equations <- function(vars, nu, theta1, theta2, r, alpha) {
    c1 <- vars[1]
    c2 <- vars[2]

    eq1 <- pchisq(c2, nu) - pchisq(c1, nu) - (1 - alpha)
    eq2 <- pchisq(theta2 ^ r / theta1 ^ r * c2, nu) - pchisq(theta2 ^ r / theta1 ^ r * c1, nu) - (1 - alpha)

    return(c(eq1, eq2))
  }

  found_solution <- FALSE

  for (i in seq_len(max_iterations)) {
    b <- qchisq(1 - alpha + pchisq(a, nu), nu)
    
    solution <- suppressWarnings(
      tryCatch({
        pracma::fsolve(equations, c(a, b), nu = nu, theta1 = theta1, theta2 = theta2, r = r, alpha = alpha)
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

#' Rejection Region for Interval Hypothesis Tests
#' 
#' `interval_null_rr()` computes the rejection region for $H_0 : \theta_1 < \theta < \theta_2$ vs $H_1 : \theta \leq \theta_1 \text{ or } \theta \geq \theta_2$.
#' 
#' @param theta1 A positive numeric value
#' @param theta2 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A character vector
#' @examples
#' interval_null_rr(1.75, 2.25, -2, 200, -1, 0.05)
#' @export
interval_null_rr <- function(theta1, theta2, r, n, m, alpha, a = 0, eps = 0.01, max_iterations = 5000) {
  c <- interval_null_critvals(theta1, theta2, r, n, m, alpha, a, eps, max_iterations)
  sprintf("(%f, %4f) U (%4f, %f)", -Inf, c[1], c[2], Inf)
}

#' Power for Interval Hypothesis Tests
#' 
#' `interval_null_beta()` computes statistical power at $\theta$ for interval hypothesis tests.
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
#' interval_null_beta(1.75, 2.25, 2.29, -2, 200, -1, 0.05)
#' @export
interval_null_beta <- function(theta1, theta2, theta, r, n, m, alpha, a = 0, eps = 0.01, max_iterations = 5000) {
  c <- interval_null_critvals(theta1, theta2, r, n, m, alpha, a, eps, max_iterations)
  1 - pchisq(theta ^ r / theta1 ^ r * c[2], df = 2 * n * m / r) + pchisq(theta ^ r / theta1 ^ r * c[1], df = 2 * n * m / r)
}