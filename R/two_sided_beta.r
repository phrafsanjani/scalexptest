two_sided_beta() <- function(theta0, theta1, r, n, m, alpha, initial_guess) {
  c <- two_sided_critvals(theta0, r, n, m, alpha, initial_guess)
  1 - pgamma(c[2], shape = n * m / r, scale = 1 / theta1 ^ r) + pgamma(c[1], shape = n * m / r, scale = 1 / theta1 ^ r)
}