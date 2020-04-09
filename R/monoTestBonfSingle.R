monoTestBonfSingle <- function(simultAlpha=0.05, OP_UMLE, OP_SE){
  q_cat <- length(OP_UMLE)

  OP_LB <- OP_UMLE - qnorm((1 - (simultAlpha / q_cat) / 2), 0, 1) * OP_SE
  OP_UB <- OP_UMLE + qnorm((1 - (simultAlpha / q_cat) / 2), 0, 1) * OP_SE

  OP_LB_1Lag <- c(0, OP_LB[-q_cat])
  OP_UB_1Lag <- c(0, OP_UB[-q_cat])

  OP_Labels_One_Mat <- (matrix(rep(OP_LB, q_cat), nrow = q_cat) >=
                          matrix(rep(OP_UB_1Lag, q_cat), nrow = q_cat, byrow = TRUE)) * 1

  OP_Labels_One_Mat[upper.tri(OP_Labels_One_Mat, diag = FALSE)] <- 0
  rownames(OP_Labels_One_Mat) <- 2:(q_cat + 1)

  OP_Labels_MinusOne_Mat <- (matrix(rep(OP_UB,q_cat), nrow = q_cat) <=
                               matrix(rep(OP_LB_1Lag, q_cat), nrow = q_cat, byrow = TRUE)) * -1
  OP_Labels_MinusOne_Mat[upper.tri(OP_Labels_MinusOne_Mat, diag = FALSE)] <- 0
  rownames(OP_Labels_MinusOne_Mat) <- 2:(q_cat + 1)

  directions <- unique(c(as.vector(OP_Labels_MinusOne_Mat), as.vector(OP_Labels_One_Mat)))
  testRes <- if (sum(c(-1,1)%in%directions) == 2) {"Reject H_0"} else {"Not Reject H_0"}

  list(testRes = testRes, simultAlpha=simultAlpha, indivAlpha = simultAlpha / q_cat)
}
