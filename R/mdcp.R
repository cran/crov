#' @title Monotonicity Direction Classification (MDC) procedure
#'
#' @description Fits a constrained regression model for an ordinal response with ordinal predictors
#' and possibly others, Espinosa and Hennig (2019) <https://doi.org/10.1007/s11222-018-9842-2>. The parameter estimates associated with an ordinal
#' predictor are constrained to be monotonic. If a monotonicity direction (isotonic or antitonic) is
#' not specified for an ordinal predictor by the user, then the monotonicity direction classification
#' procedure establishes it.
#' @param formula A \code{formula} to be fitted with ordinal response, one or more ordinal predictors, and possibly one or more other predictors.
#' For ordinal response and ordinal predictors use ordered factors.
#' @param data A data.frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data.frame), containing the
#' variables in \code{formula}. Neither a matrix nor an array will be accepted.
#' @param tryAllMonoDir A logical value that indicates whether one model should be fitted for each one of the possible combinations of monotonicity
#' directions. Use \code{TRUE} if none monotonicity direction is pre-specified using \code{monoDir} and the MDC procedure is not used.
#' @param monoDir Vector with monotonicity directions for the ordinal predictors to be used as constraints. Possible values for \code{monoDir} are
#' \code{TRUE} and \code{FALSE}. Use \code{TRUE} for "isotonic" and \code{FALSE} for "antitonic". The order of the elements in \code{monoDir} must be
#' the same as the order of the ordinal predictors in the object \code{formula}, i.e., the j-th element of \code{monoDir} must correspond to
#' the monotonicity direction of the j-th ordinal predictor in \code{formula}. If \code{tryAllMonoDir} and \code{monoDir} are not used (default option),
#' the monotonicity direction classification prodecure is executed to find the monotonicity directions associated to the model with the maximum log-likelihood.
#' @param CLS1 Numerical value for the confidence level to be used in the first step of the MDC procedure. This parameter is active if
#' \code{tryAllMonoDir} and \code{monoDir} are not used.
#' @param TLBS2 Numerical value for the tolerance level to be used in the second step of the MDC procedure over those ordinal
#' predictors classified as "Both" in the first step. This parameter is active if \code{tryAllMonoDir} and \code{monoDir} are not used.
#' @param TLNS2 Numerical value for the tolerance level to be used in the second step of the MDC procedure over those ordinal
#' predictors classified as "None" in the first step. This parameter is active if \code{tryAllMonoDir} and \code{monoDir} are not used.
#' @param StepSizeCLS2 Numerical value for the magnitude in which the confidence levels will be increased or decreased during the second step of
#' the MDC procedure. This parameter is active if \code{tryAllMonoDir} and \code{monoDir} are not used.
#' @param method The type of constrained method to be used among "MDCS1", "MDCS2", "MDCS3", "CMLEbonferroni", "CMLEconfReg", and "CMLEfiltered". Default value corresponds to "MDCS3".
#' @param monoTestSignLevel Significance level used when method is "CMLEbonferroni" or "CMLEconfReg". Default value 0.05.
#' @return \code{MDCproc}: Data frame with the monotonicity direction classification (Isotonic, Antitonic, Both, or None) used for each
#' ordinal predictor in each one of the steps of the MDC procedure (S1, S2 and S3), together with their individual confidence levels (CL). If
#' \code{monoDir} is used, \code{MDCproc} shows the monotonicity directions in \code{monoDir}.
#' @return \code{estimates}: Vector of parameter estimates of the model.
#' @return \code{estimates_se}: Vector of standard errors of the parameter estimates of the model.
#' @return \code{log.lik}: Value of the log-likelihood of the model.
#' @return \code{allModels}: Data frame with monotonicity directions, log-likelihood and parameter estimates of all models involved in the third step of the MDC
#' procedure. If parameter \code{monoDir} is used, \code{allModels} shows these results from the model with monotonicity directions
#' used in \code{monoDir} only. If parameter \code{tryAllMonoDir} is used, \code{allModels} shows these results from all the models according to all possible
#' combinations of monotonicity directions.
#' @return \code{constrOptimRes}: List with the outcomes provided by the function \code{\link[stats:constrOptim]{constrOptim}}.
#' @return \code{UMLE}: Vector with the parameter estimates of the unconstrained version of the model.
#' @return \code{UMLE_SE}: Vector with the standard errors of the unconstrained version of the model.
#' @import VGAM
#' @importFrom gtools permutations
#' @importFrom stats constrOptim qnorm
#' @examples # Ordinal predictors: EduLevel, IncQuint, Health,
#' # Overcrowd, and NumChildren
#' # mdcpExample <- mdcp(QoL ~ EduLevel + Age + IncQuint + Gender + Health +
#' # Overcrowd + Activity + NumChildren, data = crovData,
#' # CLS1 = 0.95, TLBS2 = 0.90, TLNS2 = 0.99, StepSizeCLS2 = 0.0002)
#' # mdcpExample$MDCproc
#' # cbind("CMLE"=mdcpExample$estimates,"UMLE"=mdcpExample$UMLE)
#' # mdcpExample$UMLE_SE
#' # mdcpExample$log.lik
#' # mdcpExample$allModels[1:6]
#' @seealso \code{\link[crov:monoTestBonf]{monoTestBonf}}, \code{\link[crov:monoTestConfReg]{monoTestConfReg}}, \code{\link[stats:constrOptim]{constrOptim}}.
#' @references Espinosa, J., Hennig, C. A constrained regression model for an ordinal response
#' with ordinal predictors. Stat Comput 29, 869-890 (2019). https://doi.org/10.1007/s11222-018-9842-2.
#' @export


mdcp<-function(formula, data = NULL, tryAllMonoDir = FALSE, monoDir = NULL,
               CLS1 = 0.95, TLBS2 = 0.85, TLNS2 = 0.999, StepSizeCLS2 = 0.0001,
               method = NULL, monoTestSignLevel = 0.05){

  matmf.aux <- model.frame(formula, drop.unused.levels = TRUE, data = data)
  if (attr(attr(matmf.aux,"terms"),"dataClasses")[1]!="ordered") stop("No ordinal response, use ordered factors.")
  OP.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "ordered")[-1]
  if (length(OP.ID.aux) == 0) warning("No ordinal predictors in the original data. Unconstrained model was fitted (see vglm() from package VGAM).")
  if (!is.null(monoDir) & length(OP.ID.aux)!=length(monoDir)) stop("The number of elements in monoDir must be the same as the
                                                                   number of ordinal predictors.")

  if (is.null(method)){method <- "MDCS3"}

  if (sum(method==c("MDCS1","MDCS2","MDCS3","CMLEbonferroni","CMLEconfReg","CMLEfiltered"))==0) stop("Method does not exist.")
  if (method!="MDCS3" & ( tryAllMonoDir==TRUE | !is.null(monoDir)) ) stop("This method is not compatible with imposing monotonicity to all OPs (tryAllMonoDir=TRUE or is.null(monoDir)=FALSE).")
  if (tryAllMonoDir == TRUE) {
    CLS1 <- 1
    TLBS2 <- 1
    StepSizeCLS2 <- 0 }
  if (tryAllMonoDir == TRUE & !is.null(monoDir)) stop("tryAllMonoDir and monoDir cannot be specified at the same time, one or none of them should be used.")
  if (!is.null(monoDir) & class(monoDir) != "logical") stop("monoDir must contain logical values only: TRUE or FALSE. Use TRUE for 'Isotonic' and FALSE for 'Antitonic'")


  if (method=="MDCS1") {
    resMDCS1      <- mdcs1(formula=formula, data = data, CLS1 = CLS1)
    mdcs1results  <- resMDCS1$mdcs1results
    newData <- resMDCS1$newData
    finalOPsNames <- resMDCS1$namesOPs
  }

  if (method=="MDCS2") {
    resMDCS2      <- mdcs2(formula, data = data, CLS1 = CLS1,
                           TLBS2 = TLBS2, TLNS2 = TLNS2, StepSizeCLS2 = StepSizeCLS2)
    mdcs2results  <- resMDCS2$mdcs2results
    newData <- resMDCS2$newData
    finalOPsNames <- resMDCS2$namesOPs
  }

  if ( method=="MDCS3" | is.null(method)) {
    newData       <- data
    finalOPsNames <- names(OP.ID.aux)
  }

  if (method=="CMLEbonferroni") {
    resCMLEbonferroni      <- cmleBonferroni(formula, data = data,
                                             signLevelBonferroni = monoTestSignLevel)
    CMLEbonferroniresults  <- resCMLEbonferroni$cmleBonferroniresults
    newData             <- resCMLEbonferroni$newData
    finalOPsNames       <- resCMLEbonferroni$namesOPs
  }

  if (method=="CMLEconfReg") {
    resCMLEconfReg      <- cmleConfReg(formula, data = data,
                                       signLevelConfReg = monoTestSignLevel)
    CMLEconfRegresults  <- resCMLEconfReg$cmleConfRegresults
    newData             <- resCMLEconfReg$newData
    finalOPsNames       <- resCMLEconfReg$namesOPs
  }

  if (method=="CMLEfiltered") {
    resCMLEfiltered      <- cmleFiltered(formula=formula, data = data, CLS1 = CLS1)
    CMLEfilteredresults  <- resCMLEfiltered$cmleFilteredresults
    newData <- resCMLEfiltered$newData
    finalOPsNames <- resCMLEfiltered$namesOPs
  }

  matmf.aux <- model.frame(formula, drop.unused.levels = TRUE, data = newData)
  OP.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "ordered")[-1]
  if (length(OP.ID.aux) == 0) warning("No ordinal predictors after dropping constraints. Unconstrained model was fitted (see vglm() from package VGAM).")
  if (length(OP.ID.aux) == 0) {

    vglm.result <- vglm(formula=formula, family = cumulative(parallel = TRUE),
                        data = newData)

    param_Names <- names(coefficients(vglm.result))
    param_UMLE <- unname(coefficients(vglm.result))
    se_UMLE <- unname(coef(summary(vglm.result))[,2])
    UMLE <- as.vector(t(param_UMLE))
    UMLE_SE <- as.vector(t(se_UMLE))

    list(MDCproc = NA,
         OPsNames = NA,
         estimates = UMLE,
         estimates_se = UMLE_SE,
         log.lik = (summary(vglm.result))@criterion$loglikelihood,
         allModels = NA,
         constrOptimRes = NA,
         UMLE = UMLE,
         UMLE_SE = UMLE_SE,
         UMLE_logLik = (summary(vglm.result))@criterion$loglikelihood
    )


  } else {

    NonOP.Numer.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "numeric")
    NonOP.Nomin.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "factor")

    for(col in c(1, OP.ID.aux)) {matmf.aux[, col] <- factor(matmf.aux[, col], ordered = FALSE)}

    matX.aux <- model.matrix(formula, matmf.aux)[, -1]
    originalLocationX <- colnames(matX.aux)

    matmf <- matmf.aux[, c(1, OP.ID.aux, NonOP.Numer.ID.aux, NonOP.Nomin.ID.aux)]

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

    if (is.null(monoDir)) {
      id_from_OP <- c(1,(cumsum(unname(q_cat_OrdPred - 1)) + 1)[-length(q_cat_OrdPred)])
      id_to_OP <- cumsum(unname(q_cat_OrdPred - 1))

      MDC_IndexS1 <- factor(levels = c("Isotonic", "Antitonic", "Both", "None"))
      for (OP in 1 : length(q_cat_OrdPred)) {
        OP_UMLE <- matrix(rep(t(param_UMLE_OP[id_from_OP[OP] : id_to_OP[OP]]), 1), byrow = TRUE, nrow = 1, ncol = (q_cat_OrdPred[OP] - 1))
        OP_SE <- matrix(rep(t(se_UMLE_OP[id_from_OP[OP] : id_to_OP[OP]]), 1), byrow = TRUE, nrow = 1, ncol = (q_cat_OrdPred[OP] - 1))

        OP_LB <- OP_UMLE - qnorm(CLS1+(1-CLS1)/2, 0, 1) * OP_SE
        OP_UB <- OP_UMLE + qnorm(CLS1+(1-CLS1)/2, 0, 1) * OP_SE

        OP_LB_1Lag <- c(0, OP_LB[-(q_cat_OrdPred[OP] - 1)])
        OP_UB_1Lag <- c(0, OP_UB[-(q_cat_OrdPred[OP] - 1)])

        OP_Labels_One_Mat <- (matrix(rep(OP_LB, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1) >=
                                matrix(rep(OP_UB_1Lag, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1, byrow = TRUE)) * 1
        OP_Labels_One_Mat[upper.tri(OP_Labels_One_Mat,diag = FALSE)] <- 0
        rownames(OP_Labels_One_Mat) <- 2 : q_cat_OrdPred[OP]

        OP_Labels_MinusOne_Mat <- (matrix(rep(OP_UB, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1) <=
                                     matrix(rep(OP_LB_1Lag, q_cat_OrdPred[OP] - 1), nrow = q_cat_OrdPred[OP] - 1, byrow = TRUE)) * -1
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

      increasing <- MDC$MDC_S2 == "Isotonic"
      increasing[which(MDC$MDC_S2 %in% c("Both","None"))] <- NA

      if (length(increasing[is.na(increasing) == TRUE]) >= 1) {
        increasing.mat <- matrix(rep(increasing, 2 ^ length(increasing[is.na(increasing) == TRUE])),
                                 nrow = 2 ^ length(increasing[is.na(increasing) == TRUE]),
                                 ncol = length(increasing), byrow = TRUE)
        if (length(increasing[is.na(increasing) == TRUE]) == 1) {
          increasing.mat[, which(is.na(increasing))] <- c(TRUE, FALSE)
        } else{
          increasing.mat[, which(is.na(increasing))] <- matrix(
            as.logical(permutations( n = 2, r = length(increasing[is.na(increasing) == TRUE]),
                                     v = 1 : length(increasing[is.na(increasing) == TRUE]), repeats.allowed = TRUE) - 1),
            nrow = 2 ^ length(increasing[is.na(increasing) == TRUE]),
            ncol = length(increasing[is.na(increasing) == TRUE]))
        }
      } else {
        increasing.mat <- matrix(increasing, nrow = 1, ncol = length(increasing))
      }

      num_models <- dim(increasing.mat)[1]
      ModelsRes <- as.data.frame(increasing.mat)
      colnames(ModelsRes) <- paste(rep("OP", length(q_cat_OrdPred)), 1 : length(q_cat_OrdPred), sep="")
      rownames(ModelsRes) <- paste(rep("Model", num_models), 1 : num_models, sep="")
      ModelsRes[ModelsRes == TRUE] <- "Isotonic"
      ModelsRes[ModelsRes == FALSE] <- "Antitonic"

      param_cmle <- matrix(NA, num_models, length(param_Names))
      colnames(param_cmle) <- param_Names
      ModelsRes <- cbind.data.frame(ModelsRes, log.likelihood.cmle = NA, param_cmle)
      MDC <- cbind.data.frame(MDC,MDC_S3 = MDC$MDC_S2)

    } else {

      num_models <- 1
      increasing.mat <- matrix(monoDir, nrow = 1, ncol = length(monoDir))

      ModelsRes <- as.data.frame(increasing.mat)
      colnames(ModelsRes) <- paste(rep("OP", length(q_cat_OrdPred)), 1 : length(q_cat_OrdPred), sep="")
      rownames(ModelsRes) <- paste(rep("Model", num_models), 1 : num_models, sep="")
      ModelsRes[ModelsRes == TRUE] <- "Isotonic"
      ModelsRes[ModelsRes == FALSE] <- "Antitonic"

      param_cmle <- matrix(NA, num_models, length(param_Names))
      colnames(param_cmle) <- param_Names
      ModelsRes <- cbind.data.frame(ModelsRes, log.likelihood.cmle = NA, param_cmle)

    }

    res.aux<-list()
    cmle.aux<-list()
    for (i in 1 : num_models) {
      res.aux[[i]] <- createConstraints(matY = matY, matX = matX, paramInit = param_UMLE,
                                        q_cat_OrdPred = q_cat_OrdPred, increasing = increasing.mat[i,])
      cmle.aux[[i]] <- constrOptim(res.aux[[i]]$ConstrParamInit, f = loglikOrdPred, grad = loglikGrad,
                                   ui = res.aux[[i]]$A, ci=-res.aux[[i]]$B,
                                   control = list(reltol = 1e-05, fnscale = -1),
                                   method = "BFGS", matY = matY, matX = matX,
                                   outer.iterations = 100, mu = 1e-04, outer.eps = 1e-05)
      ModelsRes$log.likelihood.cmle[i] <- cmle.aux[[i]]$value
      ModelsRes[i, -1 : -(length(q_cat_OrdPred) + 1)] <- cmle.aux[[i]]$par
    }

    best_model <- which(ModelsRes$log.likelihood.cmle == max(ModelsRes$log.likelihood.cmle))

    if (is.null(monoDir)) {
      MDC$MDC_S3 <- unname(t(ModelsRes[best_model, 1 : length(q_cat_OrdPred)]))
    } else {
      MDC <- data.frame("monoDir"=factor(rep("Isotonic",length(monoDir)), levels = c("Isotonic", "Antitonic")))
      MDC[which(monoDir == TRUE), 1] <- "Isotonic"
      MDC[which(monoDir == FALSE), 1] <- "Antitonic"
      rownames(MDC) <- names(matmf)[OP.ID.matmf]
    }

    reorderX <- c(1:(dim(matY)[2] - 1),
                  unlist(sapply(originalLocationX,
                                function(x) {which(rownames(t(ModelsRes[best_model, -1 : -(length(q_cat_OrdPred) + 1)]))==x)})))

    estimates <- as.vector(t(ModelsRes[best_model, -1 : -(length(q_cat_OrdPred) + 1)])[reorderX])
    names(estimates) <- c(rownames(t(ModelsRes[best_model, -1 : -(length(q_cat_OrdPred) + 1)]))[1:(dim(matY)[2] - 1)],
                          originalLocationX)

    hessMat <- loglikHess(matY = matY, matX = matX,
                          paramInit = as.vector(t(ModelsRes[best_model, -1 : -(length(q_cat_OrdPred) + 1)])))
    fisher_info<-solve(-hessMat)
    estimates_se <- sqrt(diag(fisher_info))
    estimates_se <- estimates_se[reorderX]
    names(estimates_se) <- names(estimates)

    UMLE <- as.vector(t(param_UMLE[reorderX]))
    UMLE_SE <- as.vector(t(se_UMLE[reorderX]))
    names(UMLE) <- names(estimates)
    names(UMLE_SE) <- names(estimates)

    constrOptimRes <- cmle.aux[[best_model]]
    constrOptimRes$par <- constrOptimRes$par[reorderX]

    list(matY = matY,
         matX = matX,
         MDCproc = MDC,
         OPsNames = finalOPsNames,
         estimates = estimates,
         estimates_se = estimates_se,
         log.lik = ModelsRes[best_model, (length(q_cat_OrdPred) + 1)],
         allModels = ModelsRes[,c(1:(dim(MDC)[1]+1),reorderX+(dim(MDC)[1]+1))],
         constrOptimRes = constrOptimRes,
         UMLE = UMLE,
         UMLE_SE = UMLE_SE,
         UMLE_logLik = (summary(vglm.result))@criterion$loglikelihood
    )

  }

}
