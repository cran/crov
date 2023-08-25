
#' @title Monotonicity test using confidence regions
#'
#' @description Tests the null hypothesis of monotonicity over a set of parameters associated to an ordinal predictor. The log-likelihood ratio test is used after imposing ordinal constraints on the parameter estimates of a single ordinal predictor and comparing its results against the unconstrained MLEs.
#' @param formula A \code{formula} to be fitted with ordinal response, one or more ordinal predictors, and possibly one or more other predictors.
#' For ordinal response and ordinal predictors use ordered factors.
#' @param data A data.frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data.frame), containing the
#' variables in \code{formula}. Neither a matrix nor an array will be accepted.
#' @param monoDir Vector with monotonicity directions for the ordinal predictors to be used as constraints. Possible values for \code{monoDir} are
#' \code{1}, \code{0} and \code{-1}. Use \code{1} for "isotonic", \code{-1} for "antitonic", and \code{0} o test monotonicity of the prameters of an ordinal
#' predictor. The order of the elements in \code{monoDir} must be the same as the order of the ordinal
#' predictors in the object \code{formula}, i.e., the j-th element of \code{monoDir} must correspond to
#' the monotonicity direction of the j-th ordinal predictor in \code{formula}. If \code{monoDir} is not used (default option),
#' the monotonicity of all ordinal predictors' effects are tested.
#' @param SignifLevel Numerical value for the significance level.
#' @return \code{resConfRegTest}: Data frame with columns:
#' \code{OPName}=Name of the ordinal predictor (OP),
#' \code{Num_Cat}=Number of categories of the OP,
#' \code{UMLE_logLik}=log-likelihood of the unconstrained model,
#' \code{CMLE_logLik}=log-likelihood of the constrained model using \code{\link[crov:mdcp]{mdcp}} assuming monotonicity for each OP,
#' \code{degreesOfFreedom}=degrees of freedom used in the hypothesis test,
#' \code{Statistic}=value of the statistic,
#' \code{CritValue}=critical value resulting from the statistic,
#' \code{SignifLevel}=significance level used in the test,
#' \code{P.Value}=p-value,
#' \code{RejectMonotonicity}=TRUE if the null hypothesis of monotonicity is rejected, FALSE otherwise.
#' @examples # Ordinal predictors: EduLevel, IncQuint and Health
#'   monoTestConfRegExample <- monoTestConfReg(QoL ~ EduLevel + Age + IncQuint + Gender +
#'   Health, data = crovData, monoDir=c(0,-1,-1), SignifLevel = 0.05)
#'   monoTestConfRegExample$resConfRegTest
#' @seealso \code{\link[crov:mdcp]{mdcp}},
#' \code{\link[crov:monoTestBonf]{monoTestBonf}},
#' \code{\link[crov:confRegUCRandUCCR]{confRegUCRandUCCR}},
#' \code{\link[crov:confRegCCR]{confRegCCR}},
#' \code{\link[crov:plotCMLE]{plotCMLE}},
#' \code{\link[VGAM:vglm]{vlgm}}.
#' @importFrom stats qnorm qchisq pchisq
#' @export
#'

monoTestConfReg<-function(formula, data = NULL, monoDir = NULL, SignifLevel=0.05){

  matmf.aux <- model.frame(formula, drop.unused.levels = TRUE, data = data)
  OP.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "ordered")[-1]

  if ( !is.null(monoDir) ) {
    names(monoDir) <- names(OP.ID.aux)
    monoDirText <- factor(levels = c("Isotonic", "Antitonic", "Estimate"))
    monoDirText[monoDir==-1] <- "Antitonic"
    monoDirText[monoDir==0]  <- "Estimate"
    monoDirText[monoDir==1]  <- "Isotonic"
    names(monoDirText) <- names(OP.ID.aux)
  }

  if ( !is.null(monoDir) ) { OP.ID.auxTest <- OP.ID.aux[!monoDir%in%c(-1,1)]
  } else {
    OP.ID.auxTest <- OP.ID.aux
  }

  resConfRegTest <- data.frame("OPName"=names(OP.ID.aux),
                               "Num_Cat"=if(is.null(dim(matmf.aux[,OP.ID.aux]))){
                                 length(unique(matmf.aux[,OP.ID.aux]))
                               } else {apply(matmf.aux[,OP.ID.aux],2,
                                            function(x){length(unique(x))})},
                               "UMLE_logLik"=0, "CMLE_logLik"=0,
                               "degreesOfFreedom"=0,"Statistic"=0,"CritValue"=0,
                               "SignifLevel"=0,"P.Value"=0,
                               "RejectMonotonicity"=FALSE)

  for (j in 1:length(OP.ID.auxTest)) {
    matmf.aux1OP <- matmf.aux
    for (j2 in OP.ID.auxTest[-j]) {
      matmf.aux1OP[,j2] <- factor(matmf.aux[,j2], ordered=FALSE)
    }
    matmf.aux1OP <- model.frame(formula, drop.unused.levels = TRUE, data = matmf.aux1OP)
    attr(attr(matmf.aux1OP,"terms"),"dataClasses")
    OP.ID.iter <- which(attr(attr(matmf.aux1OP, "terms"), "dataClasses") == "ordered")[-1]
    #resMDCPaux <- mdcp(formula=formula, tryAllMonoDir = TRUE, data=matmf.aux1OP, method="MDCS3")
    if ( !is.null(monoDir) ) { monoDirIter <- monoDir[match(names(OP.ID.iter),names(monoDir))]
    } else {
      monoDirIter <- NULL
    }
    resMDCPaux <- mdcp(formula=formula, tryAllMonoDir = TRUE, data=matmf.aux1OP,
                       method = NULL, monoDir = monoDirIter)
    rowID <- match(names(OP.ID.auxTest)[j],rownames(resConfRegTest))
    resConfRegTest[rowID,c("UMLE_logLik","CMLE_logLik")] <- c(resMDCPaux$UMLE_logLik,
                                                              resMDCPaux$log.lik)
    resConfRegTest[rowID,"degreesOfFreedom"] <- resConfRegTest[rowID,"Num_Cat"]-1
    resConfRegTest[rowID,"Statistic"] <- 2*(resConfRegTest[rowID,"UMLE_logLik"]-
                                          resConfRegTest[rowID,"CMLE_logLik"])
    resConfRegTest[rowID,"CritValue"] <- qchisq(1-SignifLevel,
                                                df = resConfRegTest[rowID,"degreesOfFreedom"])
    resConfRegTest[rowID,"SignifLevel"] <- SignifLevel
    resConfRegTest[rowID,"P.Value"] <- 1- pchisq(resConfRegTest[rowID,"Statistic"],
                                             resConfRegTest[rowID,"degreesOfFreedom"],
                                             lower.tail = TRUE)
    resConfRegTest[rowID,"RejectMonotonicity"] <- resConfRegTest[rowID,"Statistic"]>
      resConfRegTest[rowID,"CritValue"]
  }

  if ( !is.null(monoDir) ) {
    resConfRegTest <- resConfRegTest[match(names(monoDir[monoDir==0]),
                                           rownames(resConfRegTest)),]
  }
  list(resConfRegTest = resConfRegTest)
}
