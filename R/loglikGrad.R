loglikGrad <- function(paramInit,matY,matX) {
  N <- dim(matY)[1]
  c <- dim(matY)[2]
  qbTot <- dim(matX)[2]
  paramIter <- c(-1e+10, paramInit[1:(c - 1)], 1e+10, paramInit[c : (c - 1 + qbTot)])

  z <- as.matrix(rep(1, N)) %*% paramIter[1 : (c + 1)] + matX %*% paramIter[(c + 2) : (length(paramIter))] %*% t(as.matrix(rep(1, c + 1)))
  z_j <- z[, 2 : (c + 1)]
  z_jm1 <- z[, 1 : c]
  CF_j <- ifelse(is.finite(exp(z_j) / (1 + exp(z_j))), exp(z_j) / (1 + exp(z_j)), 1)
  CF_jm1 <- exp(z_jm1) / (1 + exp(z_jm1))
  DF_j <- ifelse(is.finite(exp(z_j) / ((1 + exp(z_j)) ^ 2)), exp(z_j) / ((1 + exp(z_j)) ^ 2), 0)
  DF_jm1 <- exp(z_jm1) / ((1 + exp(z_jm1)) ^ 2)
  zOwn_j <- ifelse(is.finite((exp(z_j) - 1) / ((1 + exp(z_j)))), (exp(z_j) - 1) / ((1 + exp(z_j))), 0)
  zOwn_jm1 <- (exp(z_jm1) - 1)/((1 + exp(z_jm1)))

  gradient <- c(rep(0, c - 1 + qbTot))

  k <- 1
  for (k in 1 : (c - 1)){
    kron_jk <- matrix(0, ncol = c, nrow = N)
    kron_jm1k <- matrix(0, ncol = c, nrow = N)
    kron_jk[, k] <- 1
    kron_jm1k[, k + 1] <- 1
    der_a_k <- matY * (kron_jk * DF_j - kron_jm1k * DF_jm1) / (CF_j - CF_jm1)
    gradient[k] <- sum(sum(der_a_k))
  }

  k <- 1
  for (k in 1 : qbTot){
    der_b_k <- matY * (matX[,k] %*% t(as.matrix(rep(1, c)))) * (DF_j - DF_jm1) / (CF_j - CF_jm1)
    gradient[k + c - 1] <- sum(sum(der_b_k))
  }

  return(gradient)
}
