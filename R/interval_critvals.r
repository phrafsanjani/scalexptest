library(pracma)
library(expint)

interval_critvals() <- function(theta0, theta1, r, n, m, alpha, initial_guess) {
    s <- n * m / r
    equations <- function(vars, s, theta0, theta1, r, alpha) {
      c1 <- vars[1]
      c2 <- vars[2]
      
      if (r > 0) {
          eq1 <- gammainc(s, c1 * theta0^r) - gammainc(s, c2 * theta0^r) - alpha * gamma(s)
          eq2 <- gammainc(s, c1 * theta1^r) - gammainc(s, c2 * theta1^r) - alpha * gamma(s)
      } else {
          eq1 <- gammainc(s, c1 * theta0^r) - gammainc(s, c2 * theta0^r) - (1 - alpha) * gamma(s)
          eq2 <- gammainc(s, c1 * theta1^r) - gammainc(s, c2 * theta1^r) - (1 - alpha) * gamma(s)
      }
      
      return(c(eq1, eq2))
    }
    solution <- fsolve(equations, initial_guess, s = s, theta0 = theta0, theta1 = theta1, r = r, alpha = alpha)

    c1 <- solution$x[1]
    c2 <- solution$x[2]
  
    return(c(c1, c2))
}