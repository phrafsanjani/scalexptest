interval_rr() <- function(theta0, theta1, r, n, m, alpha, initial_guess) {
    c <- interval_critvals(theta0, theta1, r, n, m, alpha, initial_guess)
    if (r > 0)
        sprintf("(%4f, %4f)", c[1], c[2])
    else
        sprintf("(%f, %4f) U (%4f, %f)", -Inf, c[1], c[2], Inf)
}