#' Rejection Region for Interval Hypothesis Tests
#' 
#' `interval_critvals()` computes the critival values for $H_0 : \theta \leq \theta_1 \text{ or } \theta \geq \theta_2$ vs $H_1 : \theta_1 < \theta < \theta_2$.
#' 
#' @param theta1 A positive numeric value
#' @param theta2 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A numeric vector
#' @examples
#' interval_critvals(1, sqrt(3), -2, 25, -1, 0.01)
#' @export
interval_critvals <- function(theta1, theta2, r, n, m, alpha) {
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

    if (r > 0) {
      eq1 <- pchisq(2 * c2 * theta1 ^ r, nu) - pchisq(2 * c1 * theta1 ^ r, nu) - (1 - alpha)
      eq2 <- pchisq(2 * c2 * theta2 ^ r, nu) - pchisq(2 * c1 * theta2 ^ r, nu) - (1 - alpha)
    } else {
      eq1 <- pchisq(2 * c2 * theta1 ^ r, nu) - pchisq(2 * c1 * theta1 ^ r, nu) - alpha
      eq2 <- pchisq(2 * c2 * theta2 ^ r, nu) - pchisq(2 * c1 * theta2 ^ r, nu) - alpha
    }

    return(c(eq1, eq2))
  }

  eps <- 0.01
  a <- 0
  max_iterations <- 5000
  found_solution <- FALSE

  for (i in seq_len(max_iterations)) {
    if (r > 0)
      b <- qchisq(1 - alpha + pchisq(2 * a * theta1 ^ r, nu), nu) / (2 * theta1 ^ r)
    else
      b <- qchisq(alpha + pchisq(2 * a * theta1 ^ r, nu), nu) / (2 * theta1 ^ r)
    
    solution <- suppressWarnings(
      tryCatch({
        pracma::fsolve(equations, c(a, b), nu = nu, theta1 = theta1, theta2 = theta2, r = r, alpha = alpha)
      }, error = function(e) NULL)
    )
    
    if (!is.null(solution) && all(is.finite(solution$x))) {
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
#' `interval_rr()` computes the rejection region for $H_0 : \theta \leq \theta_1 \text{ or } \theta \geq \theta_2$ vs $H_1 : \theta_1 < \theta < \theta_2$.
#' 
#' @param theta1 A positive numeric value
#' @param theta2 A positive numeric value
#' @param r A non-zero numeric value
#' @param n An integer represeting the sample size
#' @param m A non-zero numeric value
#' @param alpha Numeric value between 0 and 1 representing the significance level
#' @returns A character vector
#' @examples
#' interval_rr(1, sqrt(3), -2, 25, -1, 0.01)
#' @export
interval_rr <- function(theta1, theta2, r, n, m, alpha) {
  c <- interval_critvals(theta1, theta2, r, n, m, alpha)
  if (r > 0)
    sprintf("(%f, %4f) U (%4f, %f)", -Inf, c[1], c[2], Inf)
  else
    sprintf("(%4f, %4f)", c[1], c[2])
}

#' Power for Interval Hypothesis Tests
#' 
#' `interval_beta()` computes statistical power at $\theta$ for interval hypothesis tests.
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
#' interval_beta(1, sqrt(3), sqrt(2), -2, 25, -1, 0.01)
#' @export
interval_beta <- function(theta1, theta2, theta, r, n, m, alpha) {
  c <- interval_critvals(theta1, theta2, r, n, m, alpha)
  if (r > 0)
    1 - pchisq(2 * c[2] * theta ^ r, df = 2 * n * m / r) + pchisq(2 * c[1] * theta ^ r, df = 2 * n * m / r)
  else
    pchisq(2 * c[2] * theta ^ r, df = 2 * n * m / r) - pchisq(2 * c[1] * theta ^ r, df = 2 * n * m / r)
}