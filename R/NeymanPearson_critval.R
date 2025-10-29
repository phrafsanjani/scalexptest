NeymanPearson_critval() <- function(theta0, theta1, r, n, m, alpha) {
    if (theta1 ^ r > theta0 ^ r) {
      qgamma(alpha, shape = n * m / r, scale = 1 / theta0 ^ r)
    } else {
      qgamma(1 - alpha, shape = n * m / r, scale = 1 / theta0 ^ r)
    }
}