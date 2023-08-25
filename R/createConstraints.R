# @title Constraints creator
#
# @description Creates the contraints matrix and contrainst vector to be used in the optimisation problem.
# @param matY Matrix with binary values of the ordinal response variable.
# @param matX Matrix with binary values of the ordinal predictors variables and non-ordinal predictors.
# @param paramInit Vector of parameter values for intitialisation.
# @param q_cat_OrdPred Vector with the number of categories of the ordinal predictors.
# @param increasing Vector logical values. TRUE for isotonic associations and FALSE for the antitonic ones.
# @return A list of elements.
createConstraints<- function(matY, matX, paramInit, q_cat_OrdPred, increasing) {
  qxOrd <- length(q_cat_OrdPred)
  increasingPN <- increasing+(increasing - 1)
  xOrdD <- matX[, 1 : (sum(q_cat_OrdPred - 1))]
  qOrdCatPred <- rep(1 : qxOrd,q_cat_OrdPred - 1)
  indOrdCatPred <- rep(1, sum(q_cat_OrdPred - 1))
  indOrdCatPred[c(1, cumsum(q_cat_OrdPred - 1)[-qxOrd] + 1)] <- 0
  increasingOrdCat <- rep(increasing, q_cat_OrdPred - 1)
  constrOrdMat <- lapply(1 : qxOrd, function(x) {diag((q_cat_OrdPred - 1)[x]) * increasingPN[x]})
  for (nmat in 1 : qxOrd) {
    constrOrdMat[[nmat]][lower.tri(constrOrdMat[[nmat]])] <- increasingPN[nmat] * -1
    if (dim(constrOrdMat[[nmat]])[1] == 2 || dim(constrOrdMat[[nmat]])[1] == 3) {
      constrOrdMat[[nmat]][3] <- 0
    } else {
      constrOrdMat[[nmat]][which(lower.tri(constrOrdMat[[nmat]]) != 0, arr.ind = TRUE)[
        -cumsum(c(1, rev(seq(2, (q_cat_OrdPred - 1)[nmat] - 1, 1)))),]] <- 0
    }
  }

  qOrdCatPredUpDown <- qOrdCatPred * (increasingOrdCat * 2 - 1)

  N <- dim(matY)[1]
  c <- dim(matY)[2]
  qbOrd <- sum(q_cat_OrdPred - 1)
  qbTot <- dim(matX)[2]
  qbNonOrd <- qbTot-qbOrd

  A <- matrix(0, c - 1 + qbTot, c - 1 + qbTot, byrow = TRUE)
  A[c : (c - 1 + unname(q_cat_OrdPred - 1)[1]), c : (c - 1 + unname(q_cat_OrdPred - 1)[1])] <- constrOrdMat[[1]]
  if (qxOrd > 1){
    constMat <- 2
    for (constMat in 2 : qxOrd){
      A[(c + cumsum(unname(q_cat_OrdPred - 1))[constMat - 1]) : (c - 1 + cumsum(unname(q_cat_OrdPred - 1))[constMat]),
        (c + cumsum(unname(q_cat_OrdPred - 1))[constMat - 1]) : (c - 1 + cumsum(unname(q_cat_OrdPred - 1))[constMat])] <-
        constrOrdMat[[constMat]]
    }
  }
  B <- c(rep(1, c - 1), rep(0, qbOrd), rep(1, qbNonOrd))

  param_UMLE_OP <- paramInit[c : (c + qbOrd - 1)]

  id_from_OP <- c(1, (cumsum(unname(q_cat_OrdPred - 1)) + 1)[-length(q_cat_OrdPred - 1)])
  id_to_OP <- cumsum(unname(q_cat_OrdPred - 1))
  param_UDBlocks <- param_UMLE_OP
  for (OP in 1 : qxOrd) {
    param.aux <- (if(increasing[OP] == TRUE) {1} else {-1}) * c(0, param_UDBlocks[id_from_OP[OP] : id_to_OP[OP]])
    param.aux[which(param.aux < 0)] <- 0
    while (sum(param.aux[-1] >= param.aux[-length(param.aux)]) < (length(param.aux) - 1)) {
      startingBlockID <- which(cumsum(cumsum(param.aux[-1] < param.aux[-length(param.aux)])) == 1) ### 2, unique solution
      param.aux[startingBlockID : (startingBlockID + sum(cumsum(1 : (length(param.aux) - startingBlockID)) ==
                                                           cumsum(cumsum(param.aux[startingBlockID : length(param.aux)][-1] <=
                                                                           param.aux[startingBlockID : (length(param.aux) - 1)]))))]<-
        mean(param.aux[startingBlockID : (startingBlockID + sum(cumsum(1 : (length(param.aux) - startingBlockID)) ==
                                                                  cumsum(cumsum(param.aux[startingBlockID : length(param.aux)][-1] <=
                                                                                  param.aux[startingBlockID : (length(param.aux) - 1)]))))])
    }
    param.aux <- (if(increasing[OP] == TRUE) {1} else {-1}) * param.aux
    param_UDBlocks[id_from_OP[OP] : id_to_OP[OP]] <- param.aux[-1]
  }

  ConstrParamInit <- paramInit
  ConstrParamInit[c : (c + qbOrd - 1)] <- param_UDBlocks

  UpDown <- c(rep(0, c - 1), qOrdCatPredUpDown, rep(0, qbNonOrd))
  for (adjustParam in 2 : (c - 1 + qbTot)) {
    if (UpDown[adjustParam] > 0 &
        UpDown[adjustParam] == UpDown[adjustParam - 1] &
        ConstrParamInit[adjustParam] <= ConstrParamInit[adjustParam - 1]) {
      ConstrParamInit[adjustParam] <- ConstrParamInit[adjustParam - 1] + 1e-9
    }
    if (UpDown[adjustParam] > 0 &
        ConstrParamInit[adjustParam] <= 0) {
      ConstrParamInit[adjustParam] <-  1e-9
    }
    if (UpDown[adjustParam] < 0 &
        UpDown[adjustParam] == UpDown[adjustParam - 1] &
        ConstrParamInit[adjustParam] >= ConstrParamInit[adjustParam - 1]) {
      ConstrParamInit[adjustParam] <- ConstrParamInit[adjustParam - 1] - 1e-9
    }
    if (UpDown[adjustParam] < 0&
        ConstrParamInit[adjustParam] >= 0) {
      ConstrParamInit[adjustParam] <-  -1e-9
    }
  }
  ConstrParamInit
  list(A = A, B = B, ConstrParamInit = ConstrParamInit)
}
