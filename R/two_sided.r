two_sided_critvals <- function(theta0, r, n, m, alpha, initial_guess) {
    s <- n * m / r
    equations <- function(vars, s, theta0, r, alpha) {
      c1 <- vars[1]
      c2 <- vars[2]
      
      eq1 <- expint::gammainc(s, c1 * theta0 ^ r) - expint::gammainc(s, c2 * theta0 ^ r) - (1 - alpha) * gamma(s)
      eq2 <- expint::gammainc(1 + s, c1 * theta0 ^ r) - expint::gammainc(1 + s, c2 * theta0 ^ r) - (1 - alpha) * s * gamma(s)

      return(c(eq1, eq2))
    }
    solution <- pracma::fsolve(equations, initial_guess, s = s, theta0 = theta0, r = r, alpha = alpha)

    c1 <- solution$x[1]
    c2 <- solution$x[2]

    return(c(c1, c2))
}

#' Rejection Region for Two-Sided Hypothesis Tests
#' 
#' `two_sided_rr()` computes the rejection region for simple vs two-sided hypotheses.
#' 
#' @param theta0 $\theta_0$
#' @param r $r$
#' @param n $n$
#' @param m $m$
#' @param alpha Significance level $\alpha$
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
#' @param theta0 $\theta_0$
#' @param theta1 $\theta_1$
#' @param r $r$
#' @param n $n$
#' @param m $m$
#' @param alpha Significance level $\alpha$
#' @returns A probability value between 0 and 1
#' @examples
#' two_sided_beta(1, -1, 1, -1, 0.05, c(0.05, 5))
#' two_sided_beta(1, -2, 4, -1, 0.05, c(0.3, 6))
#' @export
two_sided_beta <- function(theta0, theta1, r, n, m, alpha, initial_guess) {
  c <- two_sided_critvals(theta0, r, n, m, alpha, initial_guess)
  1 - pgamma(c[2], shape = n * m / r, scale = 1 / theta1 ^ r) + pgamma(c[1], shape = n * m / r, scale = 1 / theta1 ^ r)
}