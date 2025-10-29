interval_critvals <- function(theta0, theta1, r, n, m, alpha, initial_guess) {
    s <- n * m / r
    equations <- function(vars, s, theta0, theta1, r, alpha) {
      c1 <- vars[1]
      c2 <- vars[2]
      
      if (r > 0) {
          eq1 <- expint::gammainc(s, c1 * theta0^r) - expint::gammainc(s, c2 * theta0^r) - alpha * gamma(s)
          eq2 <- expint::gammainc(s, c1 * theta1^r) - expint::gammainc(s, c2 * theta1^r) - alpha * gamma(s)
      } else {
          eq1 <- expint::gammainc(s, c1 * theta0^r) - expint::gammainc(s, c2 * theta0^r) - (1 - alpha) * gamma(s)
          eq2 <- expint::gammainc(s, c1 * theta1^r) - expint::gammainc(s, c2 * theta1^r) - (1 - alpha) * gamma(s)
      }
      
      return(c(eq1, eq2))
    }
    solution <- pracma::fsolve(equations, initial_guess, s = s, theta0 = theta0, theta1 = theta1, r = r, alpha = alpha)

    c1 <- solution$x[1]
    c2 <- solution$x[2]
  
    return(c(c1, c2))
}

#' Rejection Region for Interval Hypothesis Tests
#' 
#' `interval_rr()` computes the rejection region for $H_0 : \theta \leq \theta_1 \text{ or } \theta \geq \theta_2$ vs $H_1 : \theta_1 < \theta < \theta_2$.
#' 
#' @param theta0 $\theta_0$
#' @param r $r$
#' @param n $n$
#' @param m $m$
#' @param alpha Significance level $\alpha$
#' @returns A character vector
#' @export
interval_rr <- function(theta0, theta1, r, n, m, alpha, initial_guess) {
    c <- interval_critvals(theta0, theta1, r, n, m, alpha, initial_guess)
    if (r > 0)
        sprintf("(%4f, %4f)", c[1], c[2])
    else
        sprintf("(%f, %4f) U (%4f, %f)", -Inf, c[1], c[2], Inf)
}

#' Power for Interval Hypothesis Tests
#' 
#' `interval_beta()` computes statistical power at $\theta = \theta_1$ for interval hypothesis tests.
#' 
#' @param theta0 $\theta_0$
#' @param theta1 $\theta_1$
#' @param r $r$
#' @param n $n$
#' @param m $m$
#' @param alpha Significance level $\alpha$
#' @returns A probability value between 0 and 1
#' @export
inetrval_beta <- function(theta0, theta1, r, n, m, alpha, initial_guess) {
    c <- interval_critvals(theta0, theta1, r, n, m, alpha, initial_guess)
    if (r > 0)
        pgamma(c[2], shape = n * m / r, scale = 1 / theta1 ^ r) - pgamma(c[1], shape = n * m / r, scale = 1 / theta1 ^ r)
    else
        1 - pgamma(c[2], shape = n * m / r, scale = 1 / theta1 ^ r) + pgamma(c[1], shape = n * m / r, scale = 1 / theta1 ^ r)
}