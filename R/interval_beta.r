inetrval_beta() <- function(theta0, theta1, r, n, m, alpha, initial_guess) {
    c <- interval_critvals(theta0, theta1, r, n, m, alpha, initial_guess)
    if (r > 0)
        pgamma(c[2], shape = n * m / r, scale = 1 / theta1 ^ r) - pgamma(c[1], shape = n * m / r, scale = 1 / theta1 ^ r)
    else
        1 - pgamma(c[2], shape = n * m / r, scale = 1 / theta1 ^ r) + pgamma(c[1], shape = n * m / r, scale = 1 / theta1 ^ r)
}