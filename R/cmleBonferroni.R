cmleBonferroni <- function(formula, data = NULL, monoDir = NULL,
                           signLevelBonferroni = 0.95 ) {

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
# colnames(matmf)
  vglm.result <- vglm(OrdResp~. , family = cumulative(parallel = TRUE), data = matmf[, -1])
  param_Names <- names(coefficients(vglm.result))
  param_UMLE <- unname(coefficients(vglm.result))
  se_UMLE <- unname(coef(summary(vglm.result))[,2])

  param_Names_OP <- param_Names[length(unique(matmf[,1])) : (length(unique(matmf[,1])) + sum(q_cat_OrdPred) - length(q_cat_OrdPred) - 1)]
  param_UMLE_OP  <-  param_UMLE[length(unique(matmf[,1])) : (length(unique(matmf[,1])) + sum(q_cat_OrdPred) - length(q_cat_OrdPred) - 1)]
  se_UMLE_OP     <-     se_UMLE[length(unique(matmf[,1])) : (length(unique(matmf[,1])) + sum(q_cat_OrdPred) - length(q_cat_OrdPred) - 1)]

  OP_ID_from <- cumsum(c(1,unname(q_cat_OrdPred)[-length(q_cat_OrdPred)]-1))
  OP_ID_to <- cumsum(q_cat_OrdPred-1)

  resMonoTest.aux <- data.frame("OPname"=NA,
                                "NumParam"=rep(as.integer(NA),length(OP_ID_to)),
                                "indivAplha"=rep(as.double(NA),length(OP_ID_to)),
                                "simultAlpha"=rep(as.double(NA),length(OP_ID_to)),
                                "simultPvalue"=rep(as.double(NA),length(OP_ID_to)),
                                "RejectMonotonicity"=rep(as.logical(NA),length(OP_ID_to)),
                                "testRes"=NA)
  if (!is.null(monoDir)) {
    start <- length(q_cat_OrdPred)-sum(monoDir==0)+1
  } else {
    start <- 1
  }

  for (i in start : length(q_cat_OrdPred)) {
    resMonoTest.single <- monoTestBonf(simultAlpha = signLevelBonferroni,
                                       OP_UMLE=param_UMLE_OP[OP_ID_from[i]:OP_ID_to[i]],
                                       OP_SE=se_UMLE_OP[OP_ID_from[i]:OP_ID_to[i]])
    resMonoTest.aux$OPname[i] <- names(OP_ID_to)[i]
    resMonoTest.aux$NumParam[i] <- q_cat_OrdPred[i]-1
    resMonoTest.aux$indivAplha[i] <- resMonoTest.single$indivAlpha
    resMonoTest.aux$simultAlpha[i] <- resMonoTest.single$simultAlpha
    resMonoTest.aux$simultPvalue[i] <- resMonoTest.single$simultPvalue
    resMonoTest.aux$RejectMonotonicity[i] <- resMonoTest.single$simultPvalue <
      resMonoTest.single$simultAlpha
    resMonoTest.aux$testRes[i] <- resMonoTest.single$testRes
  }

  resMonoTest.aux

  namesOPs <-
    resMonoTest.aux[which(resMonoTest.aux[,"RejectMonotonicity"]==FALSE),"OPname"]

  if ( sum(monoDir!=0)>0 ) {namesOPs <- c(names(monoDir[monoDir!=0]),namesOPs)}

  newData <- data

  if(length(namesOPs)==0) {toBeUnordered <- OP.ID.aux} else {
    toBeUnordered <- resMonoTest.aux[which(resMonoTest.aux[,"RejectMonotonicity"]==TRUE),"OPname"]
  }

  for (i in toBeUnordered) {newData[,i] <- factor(newData[,i],ordered=FALSE)}

  return(list(cmleBonferroniresults=resMonoTest.aux,
              namesOPs=namesOPs,
              newData=newData))

}
