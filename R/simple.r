simple_critval() <- function(theta0, theta1, r, n, m, alpha) {
    if (theta1 ^ r > theta0 ^ r) {
      qgamma(alpha, shape = n * m / r, scale = 1 / theta0 ^ r)
    } else {
      qgamma(1 - alpha, shape = n * m / r, scale = 1 / theta0 ^ r)
    }
}

#' @export
simple_rr() <- function(theta0, theta1, r, n, m, alpha) {
    c <- simple_critval(theta0, theta1, r, n, m, alpha)
    if (theta1 ^ r > theta0 ^ r)
        sprintf("(%f, %4f)", -Inf, c)
    else
        sprintf("(%4f, %f)", c, Inf)
}

#' @export
simple_beta() <- function(theta0, theta1, r, n, m, alpha) {
    c <- simple_critval(theta0, theta1, r, n, m, alpha)
    if (theta1 ^ r > theta0 ^ r)
        pgamma(c, shape = n * m / r, scale = 1 / theta1 ^ r)
    else
        1 - pgamma(c, shape = n * m / r, scale = 1 / theta1 ^ r)
}