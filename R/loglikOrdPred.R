loglikOrdPred <- function(paramInit,matY,matX) {
  N <- dim(matY)[1]
  c <- dim(matY)[2]
  qbTot <- dim(matX)[2]
  paramIter <- c(-1e+10, paramInit[1 : (c - 1)], 1e+10, paramInit[c : (c - 1 + qbTot)])

  z <- as.matrix(rep(1, N)) %*% paramIter[1 : (c + 1)] + matX %*% paramIter[(c + 2) : (length(paramIter))] %*% t(as.matrix(rep(1, c + 1)))
  z_j <- z[, 2 : (c + 1)]
  z_jm1 <- z[, 1 : c]
  CF_j <- ifelse(is.finite(exp(z_j) / (1 + exp(z_j))), exp(z_j) / (1 + exp(z_j)), 1)
  CF_jm1 <- exp(z_jm1) / (1 + exp(z_jm1))

  return(sum(matY * ifelse(is.finite(log(CF_j - CF_jm1)), log(CF_j - CF_jm1), -1e+5)))
}

