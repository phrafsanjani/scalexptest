simple_critval <- function(theta0, theta1, r, n, m, alpha) {
    if (theta1 ^ r > theta0 ^ r) {
      qgamma(alpha, shape = n * m / r, scale = 1 / theta0 ^ r)
    } else {
      qgamma(1 - alpha, shape = n * m / r, scale = 1 / theta0 ^ r)
    }
}

#' Rejection Region for Simple and One-Sided Hypotheses
#' 
#' `simple_rr()` computes the rejection region for simple vs simple,
#' simple vs one-sided and one-sided vs one-sided hypotheses.
#' 
#' @param theta0 $\theta_0$
#' @param theta1 $\theta_1$
#' @param r $r$
#' @param n $n$
#' @param m $m$
#' @param alpha Significance level $\alpha$
#' @returns A character vector
#' @examples
#' simple_rr(3, 5, -2, 25, -1, 0.05)
#' simple_rr(3, 2, -2, 25, -1, 0.05)
#' simple_rr(4, 12, -2, 25, -1, 0.05)
#' simple_rr(2, 12, -2, 28, -1, 0.05)
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
#' @param theta0 $\theta_0$
#' @param theta1 $\theta_1$
#' @param r $r$
#' @param n $n$
#' @param m $m$
#' @param alpha Significance level $\alpha$
#' @returns A probability value between 0 and 1
#' @examples
#' simple_beta(3, 5, -2, 25, -1, 0.05)
#' simple_beta(3, 2, -2, 25, -1, 0.05)
#' simple_beta(4, 12, -2, 25, -1, 0.05)
#' simple_beta(2, 12, -2, 28, -1, 0.05)
#' @export
simple_beta <- function(theta0, theta1, r, n, m, alpha) {
    c <- simple_critval(theta0, theta1, r, n, m, alpha)
    if (theta1 ^ r > theta0 ^ r)
        pgamma(c, shape = n * m / r, scale = 1 / theta1 ^ r)
    else
        1 - pgamma(c, shape = n * m / r, scale = 1 / theta1 ^ r)
}