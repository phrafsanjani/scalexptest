library(pracma)
library(expint)

two_sided_critvals() <- function(theta0, r, n, m, alpha, initial_guess) {
    s <- n * m / r
    equations <- function(vars, s, theta0, r, alpha) {
      c1 <- vars[1]
      c2 <- vars[2]
      
      eq1 <- gammainc(s, c1 * theta0 ^ r) - gammainc(s, c2 * theta0 ^ r) - (1 - alpha) * gamma(s)
      eq2 <- gammainc(1 + s, c1 * theta0 ^ r) - gammainc(1 + s, c2 * theta0 ^ r) - (1 - alpha) * s * gamma(s)

      return(c(eq1, eq2))
    }
    solution <- fsolve(equations, initial_guess, s = s, theta0 = theta0, r = r, alpha = alpha)

    c1 <- solution$x[1]
    c2 <- solution$x[2]

    return(c(c1, c2))
}