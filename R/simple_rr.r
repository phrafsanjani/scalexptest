simple_rr() <- function(theta0, theta1, r, n, m, alpha) {
    if (theta1 ^ r > theta0 ^ r) {
        c <- qgamma(alpha, shape = n * m / r, scale = 1 / theta0 ^ r)
        sprintf("(%f, %4f)", -Inf, c)
    } else {
        c <- qgamma(1 - alpha, shape = n * m / r, scale = 1 / theta0 ^ r)
        sprintf("(%4f, %f)", c, Inf)
    }
}