#' @title Monotonicity Direction Classification (MDC) procedure
#'
#' @description Fits a constrained regression model for an ordinal response with ordinal predictors
#' and possibly others, Espinosa and Hennig (2019) <DOI:10.1007/s11222-018-9842-2>. The parameter estimates associated with an ordinal
#' predictor are constrained to be monotonic. If a monotonicity direction (isotonic or antitonic) is
#' not specified for an ordinal predictor (OP) by the user, then a constrained method to be indicated in the option \code{method} establishes it or the
#' approach that tries all possible combinations of monotonicity directions an the chooses the one with maximum likelihood.
#' @param formula A \code{formula} to be fitted with ordinal response, one or more ordinal predictors, and possibly one or more other predictors.
#' For ordinal response and ordinal predictors use ordered factors.
#' @param data A data.frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data.frame), containing the
#' variables in \code{formula}. Neither a matrix nor an array will be accepted.
#' @param tryAllMonoDir A logical value that indicates whether one model should be fitted for each one of the possible combinations of monotonicity
#' directions of the ordinal predictors' effects. Use \code{TRUE} if no constrained method is used in \code{method}.
#' @param monoDir Vector with monotonicity directions for the ordinal predictors to be used as constraints. Possible values for \code{monoDir} are
#' \code{1}, \code{0} and \code{-1}. Use \code{1} for "isotonic" and \code{-1} for "antitonic". If the monotonicity direction of the prameters of an ordinal
#' predictor has to be estimated, then use \code{0}, which also allows to drop the monotonicity
#' assumption when \code{method} is not "MDCS3".
#' The order of the elements in \code{monoDir} must be the same as the order of the ordinal predictors in the object \code{formula}, i.e., the j-th element of \code{monoDir} must correspond to
#' the monotonicity direction of the j-th ordinal predictor in \code{formula}. For example,
#' \code{monoDir=c(0,-1,-1,1,0)} indicates that the monotonicity direction of the effects of
#' the first OP will be estimated; the monotonicity direction of the effects of the second and third OPs are antitonic;
#' the monotonicity direction of the effects of the fourth OP isotonic; and the monotonicity direction
#' of the effects of the fifth OP will also be estimated. If \code{tryAllMonoDir} and \code{monoDir} are not used (default option),
#' the monotonicity direction classification prodecure is executed on all of the ordinal predictors to find the monotonicity directions associated to the
#' model with the maximum log-likelihood.
#' @param CLS1 Numerical value for the confidence level to be used in the first step of the MDC procedure. This parameter is active if
#' \code{tryAllMonoDir} is not used.
#' @param TLBS2 Numerical value for the tolerance level to be used in the second step of the MDC procedure over those ordinal
#' predictors classified as "Both" in the first step. This parameter is active if \code{tryAllMonoDir} and \code{monoDir} are not used.
#' @param TLNS2 Numerical value for the tolerance level to be used in the second step of the MDC procedure over those ordinal
#' predictors classified as "None" in the first step. This parameter is active if \code{tryAllMonoDir} and \code{monoDir} are not used.
#' @param StepSizeCLS2 Numerical value for the magnitude in which the confidence levels will be increased or decreased during the second step of
#' the MDC procedure. This parameter is active if \code{tryAllMonoDir} and \code{monoDir} are not used.
#' @param method The type of constrained method to be used among \code{"MDCS1"}, \code{"MDCS2"}, \code{"MDCS3"},
#' \code{"CMLEbonferroni"}, \code{"CMLEconfReg"}, and \code{"CMLEfiltered"}. Default value is \code{"MDCS3"}.
#' \code{"MDCS1"} uses the first step of the MDC procedure described in Espinosa, J. and Christian H. (2019)
#' to drop the monotonicity constraints on those ordinal predictors (OPs) classified as "both" or "none".
#' \code{"MDCS2"} uses the second step of the MDC procedure described in Espinosa, J. and Christian H. (2019)
#' to drop the monotonicity constraints on those OPs classified as "both" or "none".
#' \code{"MDCS3"} performs the three steps of the MDC procedure described in Espinosa, J. and Christian H. (2019)
#' and does not drop monotonicity constraints on the OPs, being the most restrictive method.
#' \code{"CMLEbonferroni"} tests the null hypothesis of monotonic effects for an OP as described in Espinosa, J. and Christian H. (2019)
#' and drops the monotonicity constraints on those OPs for which the monotonicity test rejects the null hypothesis.
#' \code{"CMLEconfReg"} tests the null hypothesis of monotonic effects for an OP based on the analysis of confidence regions
#' and drops the monotonicity constraints on those OPs for which the monotonicity test rejects the null hypothesis.
#' \code{"CMLEfiltered"} uses the first step of the MDC procedure described in Espinosa, J. and Christian H. (2019)
#' to drop the monotonicity constraints on those ordinal predictors (OPs) classified as "none".
#' @param monoTestSignLevel Significance level used when \code{method} is \code{"CMLEbonferroni"} or \code{"CMLEconfReg"}. Default value 0.05.
#' @param reltol Passed to \code{constrOptim()}.
#' @param mu Passed to \code{constrOptim()}.
#' @param outer.eps Passed to \code{constrOptim()}.
#' @return \code{MDCproc}: Data frame with the monotonicity direction classification (Isotonic, Antitonic, Both, or None) used for each
#' ordinal predictor in each one of the steps of the MDC procedure (S1, S2 and S3), together with their individual confidence levels (CL). If
#' \code{monoDir} is used, \code{MDCproc} shows the monotonicity directions in \code{monoDir}.
#' @return \code{estimates}: Vector of parameter estimates of the model.
#' @return \code{log.lik}: Value of the log-likelihood of the model.
#' @return \code{allModels}: Data frame with monotonicity directions, log-likelihood and parameter estimates of all models involved in the third step of the MDC
#' procedure. If parameter \code{monoDir} is used, \code{allModels} shows these results from the model with monotonicity directions
#' used in \code{monoDir} only. If parameter \code{tryAllMonoDir} is used, \code{allModels} shows these results from all the models according to all possible
#' combinations of monotonicity directions.
#' @return \code{constrOptimRes}: List with the outcomes provided by the function \code{\link[stats:constrOptim]{constrOptim}}.
#' @return \code{UMLE}: Vector with the parameter estimates of the unconstrained version of the model.
#' @return \code{UMLE_SE}: Vector with the standard errors of the unconstrained version of the model.
#' @return \code{q_cat_OrdPred}: Vector with the number of ordinal categories for each ordinal predictor. Values are displayed in the same order as the ordinal predictors are incorporated in \code{formula}.
#' @import VGAM
#' @importFrom gtools permutations
#' @importFrom stats constrOptim qnorm
#' @examples # Ordinal predictors: EduLevel, IncQuint and Health,
#' mdcpExample <- mdcp(QoL ~ EduLevel + Age + IncQuint + Gender + Health, data = crovData,
#' CLS1 = 0.95, TLBS2 = 0.90, TLNS2 = 0.99, StepSizeCLS2 = 0.0002)
#' mdcpExample$MDCproc
#' cbind("CMLE"=mdcpExample$estimates,"UMLE"=mdcpExample$UMLE)
#' mdcpExample$UMLE_SE
#' mdcpExample$log.lik
#' mdcpExample$allModels[1:6]
#' @seealso \code{\link[crov:monoTestBonf]{monoTestBonf}},
#' \code{\link[crov:monoTestConfReg]{monoTestConfReg}},
#' \code{\link[crov:confRegUCRandUCCR]{confRegUCRandUCCR}},
#' \code{\link[crov:confRegCCR]{confRegCCR}},
#' \code{\link[crov:plotCMLE]{plotCMLE}},
#' \code{\link[stats:constrOptim]{constrOptim}}.
#' @references Espinosa, J., and Hennig, C. "A constrained regression model for an ordinal response
#' with ordinal predictors." Statistics and Computing 29.5 (2019): 869-890.
#' https://doi.org/10.1007/s11222-018-9842-2.
#' @export


mdcp<-function(formula, data = NULL, tryAllMonoDir = FALSE, monoDir = NULL,
               CLS1 = 0.95, TLBS2 = 0.85, TLNS2 = 0.999, StepSizeCLS2 = 0.0001,
               method = NULL, monoTestSignLevel = 0.05,
               reltol = 1e-05, mu = 1e-04, outer.eps = 1e-05){

  matmf.aux <- model.frame(formula, drop.unused.levels = TRUE, data = data)
  if ( attr(attr(matmf.aux,"terms"),"dataClasses")[1]!="ordered") stop("No ordinal response, use ordered factors.")
  OP.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "ordered")[-1]

  if ( !is.null(monoDir) ) {
    names(monoDir) <- names(OP.ID.aux)
    monoDirText <- factor(levels = c("Isotonic", "Antitonic", "Estimate"))
    monoDirText[monoDir==-1] <- "Antitonic"
    monoDirText[monoDir==0]  <- "Estimate"
    monoDirText[monoDir==1]  <- "Isotonic"
    names(monoDirText) <- names(OP.ID.aux)
  }

  if ( length(OP.ID.aux) == 0) warning("No ordinal predictors in the original data. Unconstrained model was fitted (see vglm() from package VGAM).")
  if ( !is.null(monoDir) & length(OP.ID.aux)!=length(monoDir)) stop("The number of elements in monoDir must be the same as the
                                                                    number of ordinal predictors.")

  if ( tryAllMonoDir == FALSE & is.null(method) ){method <- "MDCS3"}

  if (!is.null(monoDir)) {
    if ( tryAllMonoDir == TRUE & sum(monoDir!=0)==length(monoDir) ) stop("Option tryAllMonoDir=TRUE is not compatible with establishing all of the monotonicity directions in option monoDir")
  }

  if (!is.null(method)&!is.null(monoDir)) {
    if ( method!="MDCS3" & sum(monoDir!=0)==length(monoDir) ) stop("This method is not compatible with establishing all of the monotonicity directions in option monoDir")
  }


  if ( sum(method == c("MDCS1","MDCS2","MDCS3","CMLEbonferroni","CMLEconfReg",
                     "CMLEfiltered"))==0 & !is.null(method)) stop("Method does not exist.")
  if ( tryAllMonoDir == TRUE & !is.null(method) ) stop("This method is not compatible with trying all possible combinations of monotonicity directions (tryAllMonoDir=TRUE).")
  if ( tryAllMonoDir == TRUE & is.null(method)) {
    CLS1 <- 1
    TLBS2 <- 1
    StepSizeCLS2 <- 0
    method <- "MDCS3"}
  if ( !is.null(monoDir) & sum(monoDir%in%c(-1,0,1))!=length(monoDir) ) stop("monoDir must contain values -1, 0 or 1 only. Use -1 for 'Antitonic', 1 for 'Isotonic',
                                                                             and 0 for a monotonicity directions that is needed to be estimated.")

  if ( method=="MDCS1" ) {
    resMDCS1      <- mdcs1(formula=formula, data = data, monoDir = monoDir, CLS1 = CLS1)
    mdcs1results  <- resMDCS1$mdcs1results[names(monoDir)]
    newData <- resMDCS1$newData
    finalOPsNames <- colnames(newData)[colnames(newData)%in%resMDCS1$namesOPs]
  }

  if ( method=="MDCS2" ) {
    resMDCS2      <- mdcs2(formula, data = data, monoDir = monoDir, CLS1 = CLS1,
                           TLBS2 = TLBS2, TLNS2 = TLNS2, StepSizeCLS2 = StepSizeCLS2)
    mdcs2results  <- resMDCS2$mdcs2results[names(monoDir)]
    newData <- resMDCS2$newData
    finalOPsNames <- colnames(newData)[colnames(newData)%in%resMDCS2$namesOPs]
  }

  if ( method=="MDCS3") {
    newData       <- data
    finalOPsNames <- names(OP.ID.aux)
  }

  if (method=="CMLEbonferroni") {
    resCMLEbonferroni      <- cmleBonferroni(formula, data = data, monoDir = monoDir,
                                             signLevelBonferroni = monoTestSignLevel)
    CMLEbonferroniresults  <- resCMLEbonferroni$cmleBonferroniresults
    CMLEbonferroniresults[is.na(CMLEbonferroniresults[,1]),1]  <- names(monoDir[monoDir!=0])
    CMLEbonferroniresults  <- CMLEbonferroniresults[
      match(names(monoDir),CMLEbonferroniresults[,1]),]
    newData             <- resCMLEbonferroni$newData
    finalOPsNames       <- colnames(newData)[colnames(newData)%in%resCMLEbonferroni$namesOPs]
  }

  if (method=="CMLEconfReg") {
    resCMLEconfReg      <- cmleConfReg(formula, data = data, monoDir = monoDir,
                                       signLevelConfReg = monoTestSignLevel)
    CMLEconfRegresults  <- resCMLEconfReg$cmleConfRegresults[names(monoDir)]
    newData             <- resCMLEconfReg$newData
    finalOPsNames       <- colnames(newData)[colnames(newData)%in%resCMLEconfReg$namesOPs]
  }

  if (method=="CMLEfiltered") {
    resCMLEfiltered      <- cmleFiltered(formula=formula, data = data, monoDir = monoDir,
                                         CLS1 = CLS1)
    CMLEfilteredresults  <- resCMLEfiltered$cmleFilteredresults[names(monoDir)]
    newData <- resCMLEfiltered$newData
    finalOPsNames <- colnames(newData)[colnames(newData)%in%resCMLEfiltered$namesOPs]
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

    for( col in c(1, OP.ID.aux) ) {matmf.aux[, col] <- factor(matmf.aux[, col], ordered = FALSE)}

    # colnames(matmf.aux)
    matX.aux <- model.matrix(formula, matmf.aux)[, -1]
    originalLocationX <- colnames(matX.aux)

    if ( !is.null(monoDir) ) {
      monoDir <- monoDir[names(monoDir)%in%names(OP.ID.aux)]
      monoDirText <- monoDirText[names(monoDirText)%in%names(monoDir)]
      OP.ID.aux.monoDir <- OP.ID.aux[monoDir!=0]
      OP.ID.aux.monoDir0 <- OP.ID.aux[monoDir==0]
      matmf <- matmf.aux[, c(1, OP.ID.aux.monoDir, OP.ID.aux.monoDir0, NonOP.Numer.ID.aux, NonOP.Nomin.ID.aux)]
    } else {
      matmf <- matmf.aux[, c(1, OP.ID.aux, NonOP.Numer.ID.aux, NonOP.Nomin.ID.aux)]
    }

    # colnames(matmf)
    OP.ID.matmf <- 2 : (1 + length(OP.ID.aux))
    names(OP.ID.matmf) <- names(matmf)[2 : (1 + length(OP.ID.aux))]
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

    matmfNames<-as.data.frame(matmf[,-1], col.names=colnames(matmf)[-1])
    colnames(matmfNames) <- colnames(matmf)[-1]

    vglm.result <- vglm(OrdResp~. , family = cumulative(parallel = TRUE), data = matmfNames)
    param_Names <- names(coefficients(vglm.result))
    param_UMLE <- unname(coefficients(vglm.result))
    se_UMLE <- unname(coef(summary(vglm.result))[,2])

    param_Names_OP <- param_Names[length(unique(matmf[,1])) : (length(unique(matmf[,1])) + sum(q_cat_OrdPred) - length(q_cat_OrdPred) - 1)]
    param_UMLE_OP  <-  param_UMLE[length(unique(matmf[,1])) : (length(unique(matmf[,1])) + sum(q_cat_OrdPred) - length(q_cat_OrdPred) - 1)]
    se_UMLE_OP     <-     se_UMLE[length(unique(matmf[,1])) : (length(unique(matmf[,1])) + sum(q_cat_OrdPred) - length(q_cat_OrdPred) - 1)]

    # is.null(monoDir) | sum(monoDir==0)>0
    if ( is.null(monoDir) | sum(monoDir==0)>0 ) {
      id_from_OP <- c(1,(cumsum(unname(q_cat_OrdPred - 1)) + 1)[-length(q_cat_OrdPred)])
      id_to_OP <- cumsum(unname(q_cat_OrdPred - 1))

      MDC_IndexS1 <- factor(levels = c("Isotonic", "Antitonic", "Both", "None"))
      if ( sum(monoDir==0)>0 ) {
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

      if (!is.null(monoDir)) {
        MDC_IndexS1[is.na(MDC_IndexS1)] <- monoDirText[monoDir!=0]
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
      increasing.mat <- matrix(monoDir==1, nrow = 1, ncol = length(monoDir))

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
                                   control = list(reltol = reltol, fnscale = -1),
                                   method = "BFGS", matY = matY, matX = matX,
                                   outer.iterations = 100, mu = mu, outer.eps = outer.eps)
      ModelsRes$log.likelihood.cmle[i] <- cmle.aux[[i]]$value
      ModelsRes[i, -1 : -(length(q_cat_OrdPred) + 1)] <- cmle.aux[[i]]$par
    }

    best_model <- which(ModelsRes$log.likelihood.cmle == max(ModelsRes$log.likelihood.cmle))

    # is.null(monoDir) | (sum(monoDir==0)==length(monoDir))
    if ( is.null(monoDir) | (sum(monoDir==0)==length(monoDir)) | sum(monoDir==0)>0 ) {
      MDC$MDC_S3 <- unname(t(ModelsRes[best_model, 1 : length(q_cat_OrdPred)]))
      MDC <- as.data.frame(MDC[match(names(OP.ID.aux),rownames(MDC)),])
      if ( sum(monoDir%in%c(-1,1)) >0) {MDC[monoDir%in%c(-1,1),c(2,4)] <- c(NA,NA)}
    } else {
      MDC <- data.frame("monoDir"=factor(rep("Isotonic",length(monoDir)), levels = c("Isotonic", "Antitonic")))
      MDC$monoDir <- unname(t(ModelsRes[best_model, 1 : length(q_cat_OrdPred)]))
      MDC[which(monoDir == 1), 1] <- "Isotonic"
      MDC[which(monoDir == -1), 1] <- "Antitonic"
      rownames(MDC) <- names(q_cat_OrdPred)
    }

    if (tryAllMonoDir==TRUE) {
      MDC <- as.data.frame(MDC[,5], row.names=rownames(MDC))
      colnames(MDC) <- "tryAllResults"
    }

    reorderX <- c(1:(dim(matY)[2] - 1),
                  unlist(sapply(originalLocationX,
                                function(x) {
                                  which(rownames(t(ModelsRes[best_model, -1 : -(length(q_cat_OrdPred) + 1)]))==x)
                                })))

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

    # is.null(monoDir) | (sum(monoDir==0)==length(monoDir)) | sum(monoDir==0)==0
    if ( is.null(monoDir) | (sum(monoDir==0)==length(monoDir)) | sum(monoDir==0)==0 ) {
      orderAllModels <- 1:(dim(MDC)[1]+1)
    } else {
      orderAllModels <- c(match(names(OP.ID.aux),names(c(OP.ID.aux.monoDir,OP.ID.aux.monoDir0))),
                          dim(MDC)[1]+1)
    }

    allModels <- ModelsRes[,c(orderAllModels,reorderX+(dim(MDC)[1]+1))]
    colnames(allModels)[1:dim(MDC)[1]] <- names(OP.ID.aux)

    list(matY = matY,
         matX = matX[,match(names(estimates[-1:-(dim(matY)[2]-1)]),colnames(matX))],
         MDCproc = MDC,
         OPsNames = finalOPsNames,
         estimates = estimates,
         log.lik = ModelsRes[best_model, (length(q_cat_OrdPred) + 1)],
         allModels = allModels,
         constrOptimRes = constrOptimRes,
         UMLE = UMLE,
         UMLE_SE = UMLE_SE,
         UMLE_logLik = (summary(vglm.result))@criterion$loglikelihood,
         q_cat_OrdPred = q_cat_OrdPred
    )
  }
}
