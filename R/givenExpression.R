givenExpression <- function(cD) {
    nelik = exp(rowSums(log(1 - exp(cD@locLikelihoods)), na.rm = TRUE))
    cD@posteriors <- (cD@posteriors) + log(1 - nelik)
    cD@nullPosts <- matrix(log(nelik), ncol = 1)
    cD
}
