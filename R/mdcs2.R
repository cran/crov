mdcs2 <- function(formula, data = NULL, monoDir = NULL, CLS1 = 0.95,
                  TLBS2 = 0.85, TLNS2 = 0.999, StepSizeCLS2 = 0.0001) {

  matmf.aux <- model.frame(formula, drop.unused.levels = TRUE, data = data)

  if (attr(attr(matmf.aux,"terms"),"dataClasses")[1]!="ordered") stop("No ordinal response, use ordered factors.")
  OP.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "ordered")[-1]

  if ( !is.null(monoDir) ) {
    names(monoDir) <- names(OP.ID.aux)
    monoDirText <- factor(levels = c("Isotonic", "Antitonic", "Estimate"))
    monoDirText[monoDir==-1] <- "Antitonic"
    monoDirText[monoDir==0]  <- "Estimate"
    monoDirText[monoDir==1]  <- "Isotonic"
    names(monoDirText) <- names(OP.ID.aux)
  }

  if (length(OP.ID.aux) == 0) stop("No ordinal predictors")
  NonOP.Numer.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "numeric")
  NonOP.Nomin.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "factor")

  for(col in c(1, OP.ID.aux)) {matmf.aux[, col] <- factor(matmf.aux[, col], ordered = FALSE)}

  matX.aux <- model.matrix(formula, matmf.aux)[, -1]
  originalLocationX <- colnames(matX.aux)

  if ( !is.null(monoDir) ) {
    OP.ID.aux.monoDir <- OP.ID.aux[monoDir!=0]
    OP.ID.aux.monoDir0 <- OP.ID.aux[monoDir==0]
    matmf <- matmf.aux[, c(1, OP.ID.aux.monoDir, OP.ID.aux.monoDir0, NonOP.Numer.ID.aux, NonOP.Nomin.ID.aux)]
  } else {
    matmf <- matmf.aux[, c(1, OP.ID.aux, NonOP.Numer.ID.aux, NonOP.Nomin.ID.aux)]
  }

  OP.ID.matmf <- 2 : (1 + length(OP.ID.aux))
  names(OP.ID.matmf) <- names(OP.ID.aux)
  if (length(NonOP.Numer.ID.aux) > 0) {
    NonOP.Numer.ID.matmf <- (max(OP.ID.matmf) + 1) : (max(OP.ID.matmf) + length(NonOP.Numer.ID.aux))
    names(NonOP.Numer.ID.matmf) <- names(NonOP.Numer.ID.aux)
  } else {NonOP.Numer.ID.matmf <- NULL}

  if (length(NonOP.Nomin.ID.aux) > 0) {
    NonOP.Nomin.ID.matmf <- (max(NonOP.Numer.ID.matmf,OP.ID.matmf) + 1) : (max(NonOP.Numer.ID.matmf,OP.ID.matmf) + length(NonOP.Nomin.ID.aux))
    names(NonOP.Nomin.ID.matmf) <- names(NonOP.Nomin.ID.aux)
  } else {NonOP.Nomin.ID.matmf <- NULL}


  OrdResp <- factor(matmf[, 1], ordered = TRUE)
  OrdPredMatNumeric <- sapply(matmf[, OP.ID.matmf], as.numeric)

  if (length(NonOP.Numer.ID.matmf) == 0 & length(NonOP.Nomin.ID.matmf) == 0) {
    NonOrdPredMatExtended <- NULL
  }

  if (length(NonOP.Numer.ID.matmf) == 0 & length(NonOP.Nomin.ID.matmf) > 0) {
    NonOrdPredMatExtended <- matX.aux[, c(unlist(apply(sapply(names(NonOP.Nomin.ID.matmf), function (x) {
      grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )))]
  }

  if (length(NonOP.Numer.ID.matmf) > 0 & length(NonOP.Nomin.ID.matmf) == 0) {
    NonOrdPredMatExtended <- matX.aux[, c(unlist(apply(sapply(names(NonOP.Numer.ID.matmf), function (x) {
      grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )))]
  }

  if (length(NonOP.Numer.ID.matmf) > 0 & length(NonOP.Nomin.ID.matmf) > 0) {
    NonOrdPredMatExtended <- matX.aux[, c(unlist(apply(sapply(names(NonOP.Numer.ID.matmf), function (x) {
      grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )),
      unlist(apply(sapply(names(NonOP.Nomin.ID.matmf), function (x) {
        grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )))]
  }

  if (length(NonOP.Numer.ID.matmf)+length(NonOP.Nomin.ID.matmf) == 0) {
    matX <- matX.aux[, c(unlist(apply(sapply(names(OP.ID.matmf), function (x) {
      grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )))]
  } else {
    if (length(NonOP.Numer.ID.matmf) == 0) {
      matX <- matX.aux[, c(unlist(apply(sapply(names(OP.ID.matmf), function (x) {
        grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )),
        unlist(apply(sapply(names(NonOP.Nomin.ID.matmf), function (x) {
          grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )))]
    } else {
      if (length(NonOP.Nomin.ID.matmf) == 0) {
        matX <- matX.aux[, c(unlist(apply(sapply(names(OP.ID.matmf), function (x) {
          grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )),
          unlist(apply(sapply(names(NonOP.Numer.ID.matmf), function (x) {
            grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )))]
      } else {
        matX <- matX.aux[, c(unlist(apply(sapply(names(OP.ID.matmf), function (x) {
          grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )),
          unlist(apply(sapply(names(NonOP.Numer.ID.matmf), function (x) {
            grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )),
          unlist(apply(sapply(names(NonOP.Nomin.ID.matmf), function (x) {
            grepl(x, colnames(matX.aux), fixed = TRUE)} ), 2, which )))]
      }
    }
  }

  matY <- model.matrix( ~ as.factor(matmf[,1]) - 1)

  q_cat_OrdPred <- apply(as.matrix(matmf[, OP.ID.matmf]), 2, function(x) length(unique(x)))

  vglm.result <- vglm(OrdResp~. , family = cumulative(parallel = TRUE), data = matmf[, -1])
  param_Names <- names(coefficients(vglm.result))
  param_UMLE <- unname(coefficients(vglm.result))
  se_UMLE <- unname(coef(summary(vglm.result))[,2])

  param_Names_OP <- param_Names[length(unique(matmf[,1])) : (length(unique(matmf[,1])) + sum(q_cat_OrdPred) - length(q_cat_OrdPred) - 1)]
  param_UMLE_OP  <-  param_UMLE[length(unique(matmf[,1])) : (length(unique(matmf[,1])) + sum(q_cat_OrdPred) - length(q_cat_OrdPred) - 1)]
  se_UMLE_OP     <-     se_UMLE[length(unique(matmf[,1])) : (length(unique(matmf[,1])) + sum(q_cat_OrdPred) - length(q_cat_OrdPred) - 1)]

  id_from_OP <- c(1,(cumsum(unname(q_cat_OrdPred - 1)) + 1)[-length(q_cat_OrdPred)])
  id_to_OP <- cumsum(unname(q_cat_OrdPred - 1))


  MDC_IndexS1 <- factor(levels = c("Isotonic", "Antitonic", "Both", "None"))
  if (!is.null(monoDir)) {
    start <- length(q_cat_OrdPred)-sum(monoDir==0)+1
  } else {
    start <- 1
  }

  for (OP in start : length(q_cat_OrdPred)) {
    OP_UMLE <- matrix(rep(t(param_UMLE_OP[id_from_OP[OP] : id_to_OP[OP]]), 1), byrow = TRUE, nrow = 1, ncol = (q_cat_OrdPred[OP] - 1))
    OP_SE <- matrix(rep(t(se_UMLE_OP[id_from_OP[OP] : id_to_OP[OP]]), 1), byrow = TRUE, nrow = 1, ncol = (q_cat_OrdPred[OP] - 1))

    OP_LB <- OP_UMLE - qnorm(CLS1+(1-CLS1)/2, 0, 1) * OP_SE
    OP_UB <- OP_UMLE + qnorm(CLS1+(1-CLS1)/2, 0, 1) * OP_SE

    OP_LB_1Lag <- c(0, OP_LB[-(q_cat_OrdPred[OP] - 1)])
    OP_UB_1Lag <- c(0, OP_UB[-(q_cat_OrdPred[OP] - 1)])

    OP_Labels_One_Mat <- (matrix(rep(OP_LB, q_cat_OrdPred[OP] - 1),
                                 nrow = q_cat_OrdPred[OP] - 1) >=
                            matrix(rep(OP_UB_1Lag, q_cat_OrdPred[OP] - 1),
                                   nrow = q_cat_OrdPred[OP] - 1, byrow = TRUE)) * 1
    OP_Labels_One_Mat[upper.tri(OP_Labels_One_Mat,diag = FALSE)] <- 0
    rownames(OP_Labels_One_Mat) <- 2 : q_cat_OrdPred[OP]

    OP_Labels_MinusOne_Mat <- (matrix(rep(OP_UB, q_cat_OrdPred[OP] - 1),
                                      nrow = q_cat_OrdPred[OP] - 1) <=
                                 matrix(rep(OP_LB_1Lag, q_cat_OrdPred[OP] - 1),
                                        nrow = q_cat_OrdPred[OP] - 1, byrow = TRUE)) * -1
    OP_Labels_MinusOne_Mat[upper.tri(OP_Labels_MinusOne_Mat, diag = FALSE)] <- 0
    rownames(OP_Labels_MinusOne_Mat) <- 2 : q_cat_OrdPred[OP]

    directions <- unique(c(as.vector(OP_Labels_MinusOne_Mat), as.vector(OP_Labels_One_Mat)))

    MDC_IndexS1[OP] <- if (sum(directions == 1) > 0 & sum(directions == -1) == 0) {"Isotonic"} else {
      if (sum(directions == 1) == 0 & sum(directions == -1) > 0) {"Antitonic"} else {
        if (sum(directions == 1) == 0 & sum(directions == -1) == 0) {"Both"} else {"None"}}}
  }

  names(MDC_IndexS1) <- names(matmf)[OP.ID.matmf]
  MDC <- data.frame(MDC_S1 = MDC_IndexS1, CL_S1 = rep(CLS1, length(q_cat_OrdPred)))

  MDC_IndexS2Both <- factor(levels = c("Isotonic","Antitonic","Both","None"))
  MDC_IndexS2Both <- MDC_IndexS1
  MDC <- cbind.data.frame(MDC,MDC_S2 = MDC$MDC_S1, CL_S2 = MDC$CL_S1)
  for (OP in which(MDC_IndexS1 == "Both")) {
    CL <- CLS1
    while (MDC_IndexS2Both[OP] == "Both" & CL > TLBS2) {
      CL <- round(CL - StepSizeCLS2, 9)
      OP_UMLE <- matrix(rep(t(param_UMLE_OP[id_from_OP[OP] : id_to_OP[OP]]), 1), byrow = TRUE, nrow = 1, ncol = (q_cat_OrdPred[OP] - 1))
      OP_SE <- matrix(rep(t(se_UMLE_OP[id_from_OP[OP] : id_to_OP[OP]]), 1), byrow = TRUE, nrow = 1, ncol = (q_cat_OrdPred[OP] - 1))

      OP_LB <- OP_UMLE - qnorm(CL+(1-CL)/2, 0, 1) * OP_SE
      OP_UB <- OP_UMLE + qnorm(CL+(1-CL)/2, 0, 1) * OP_SE

      OP_LB_1Lag <- c(0, OP_LB[-(q_cat_OrdPred[OP] - 1)])
      OP_UB_1Lag <- c(0, OP_UB[-(q_cat_OrdPred[OP] - 1)])

      OP_Labels_One_Mat <- (matrix(rep(OP_LB, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1) >=
                              matrix(rep(OP_UB_1Lag, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1, byrow = TRUE)) * 1
      OP_Labels_One_Mat[upper.tri(OP_Labels_One_Mat, diag = FALSE)] <- 0
      rownames(OP_Labels_One_Mat)<- 2 : q_cat_OrdPred[OP]

      OP_Labels_MinusOne_Mat <- (matrix(rep(OP_UB, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1) <=
                                   matrix(rep(OP_LB_1Lag, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1, byrow = TRUE)) * -1
      OP_Labels_MinusOne_Mat[upper.tri(OP_Labels_MinusOne_Mat, diag = FALSE)] <- 0
      rownames(OP_Labels_MinusOne_Mat) <- 2 : q_cat_OrdPred[OP]

      directions <- unique(c(as.vector(OP_Labels_MinusOne_Mat), as.vector(OP_Labels_One_Mat)))

      MDC_IndexS2Both[OP] <- if (sum(directions == 1) > 0 & sum(directions == -1) == 0) {"Isotonic"} else {
        if (sum(directions == 1) == 0 & sum(directions == -1) > 0) {"Antitonic"} else {
          if (sum(directions == 1) == 0 & sum(directions == -1) == 0) {"Both"} else {"None"}}}
      MDC$MDC_S2[OP] <- MDC_IndexS2Both[OP]
      MDC$CL_S2[OP] <- CL
    }
  }

  MDC_IndexS2None <- factor(levels=c("Isotonic","Antitonic","Both","None"))
  MDC_IndexS2None <- MDC_IndexS2Both
  for (OP in which(MDC_IndexS2Both == "None")) {
    CL <- CLS1
    while (MDC_IndexS2None[OP] == "None" & CL < TLNS2) {
      CL <- round(CL + StepSizeCLS2, 9)
      OP_UMLE <- matrix(rep(t(param_UMLE_OP[id_from_OP[OP] : id_to_OP[OP]]), 1),byrow = TRUE, nrow = 1, ncol = (q_cat_OrdPred[OP] - 1))
      OP_SE <- matrix(rep(t(se_UMLE_OP[id_from_OP[OP] : id_to_OP[OP]]), 1), byrow = TRUE, nrow = 1, ncol = (q_cat_OrdPred[OP] - 1))

      OP_LB <- OP_UMLE - qnorm(CL+(1-CL)/2, 0, 1) * OP_SE
      OP_UB <- OP_UMLE + qnorm(CL+(1-CL)/2, 0, 1) * OP_SE

      OP_LB_1Lag <- c(0, OP_LB[-(q_cat_OrdPred[OP] - 1)])
      OP_UB_1Lag <- c(0, OP_UB[-(q_cat_OrdPred[OP] - 1)])

      OP_Labels_One_Mat <- (matrix(rep(OP_LB, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1) >=
                              matrix(rep(OP_UB_1Lag, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1, byrow = T)) * 1
      OP_Labels_One_Mat[upper.tri(OP_Labels_One_Mat, diag = FALSE)] <- 0
      rownames(OP_Labels_One_Mat) <- 2 : q_cat_OrdPred[OP]

      OP_Labels_MinusOne_Mat <- (matrix(rep(OP_UB, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1) <=
                                   matrix(rep(OP_LB_1Lag, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1, byrow = TRUE)) * -1
      OP_Labels_MinusOne_Mat[upper.tri(OP_Labels_MinusOne_Mat, diag = FALSE)] <- 0
      rownames(OP_Labels_MinusOne_Mat) <- 2 : q_cat_OrdPred[OP]

      directions <- unique(c(as.vector(OP_Labels_MinusOne_Mat), as.vector(OP_Labels_One_Mat)))

      MDC_IndexS2None[OP] <- if (sum(directions == 1) > 0 & sum(directions == -1) == 0) {"Isotonic"} else {
        if (sum(directions == 1) == 0 & sum(directions == -1) > 0) {"Antitonic"} else {
          if (sum(directions == 1) == 0 & sum(directions == -1) == 0) {"Both"} else {"None"}}}
      MDC$MDC_S2[OP] <- MDC_IndexS2None[OP]
      MDC$CL_S2[OP] <- CL
    }
  }

  MDC_IndexS2 <- MDC$MDC_S2
  names(MDC_IndexS2) <- rownames(MDC)

  if (!is.null(monoDir)) {
    MDC_IndexS2[is.na(MDC_IndexS2)] <- monoDirText[monoDir!=0]
  }

  namesOPs <- names(MDC_IndexS2)[which(MDC_IndexS2%in%c("Antitonic","Isotonic"))]

  newData <- data

  if(length(namesOPs)==0) {toBeUnordered <- OP.ID.aux} else {
    toBeUnordered <- OP.ID.aux[-which(names(OP.ID.aux)%in%namesOPs)]}

  for (i in toBeUnordered) {newData[,i] <- factor(newData[,i],ordered=FALSE)}

  return(list(mdcs2results=MDC_IndexS2,
              namesOPs=namesOPs,
              newData=newData))

}
