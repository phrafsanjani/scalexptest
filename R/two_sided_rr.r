two_sided_rr() <- function(theta0, r, n, m, alpha, initial_guess) {
    c <- two_sided_critvals(theta0, r, n, m, alpha, initial_guess)
    sprintf("(%f, %4f] U [%4f, %f)", -Inf, c[1], c[2], Inf)
}