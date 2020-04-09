#' @title Monotonicity test using confidence regions
#'
#' @description Tests the null hypothesis of monotonicity over a set of parameters associated to an ordinal predictor.
#' @param formula A \code{formula} to be fitted with ordinal response, one or more ordinal predictors, and possibly one or more other predictors.
#' For ordinal response and ordinal predictors use ordered factors.
#' @param data A data.frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data.frame), containing the
#' variables in \code{formula}. Neither a matrix nor an array will be accepted.
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
#' @examples # Ordinal predictors: EduLevel, IncQuint, Health,
#'   # Overcrowd, and NumChildren
#'   # monoTestConfRegExample <- monoTestConfReg(QoL ~ EduLevel + Age + IncQuint + Gender +
#'   # Health + Overcrowd + Activity + NumChildren, data = crovData, SignifLevel = 0.05)
#'   # monoTestConfRegExample$resConfRegTest
#' @seealso \code{\link[crov:mdcp]{mdcp}}, \code{\link[crov:monoTestBonf]{monoTestBonf}}, \code{\link[VGAM:vglm]{vlgm}}.
#' @importFrom stats qnorm qchisq pchisq
#' @export
#'

monoTestConfReg<-function(formula, data = NULL, SignifLevel=0.05){

  matmf.aux <- model.frame(formula, drop.unused.levels = TRUE, data = data)
  OP.ID.aux <- which(attr(attr(matmf.aux, "terms"), "dataClasses") == "ordered")[-1]

  resConfRegTest <- data.frame("OPName"=names(OP.ID.aux),
                               "Num_Cat"=if(is.null(dim(matmf.aux[,OP.ID.aux]))){
                                 length(unique(matmf.aux[,OP.ID.aux]))
                               }else {apply(matmf.aux[,OP.ID.aux],2,
                                            function(x){length(unique(x))})},
                               "UMLE_logLik"=0, "CMLE_logLik"=0,
                               "degreesOfFreedom"=0,"Statistic"=0,"CritValue"=0,
                               "SignifLevel"=0,"P.Value"=0,
                               "RejectMonotonicity"=FALSE)
  for (j in 1:length(OP.ID.aux)) {
    matmf.aux1OP <- matmf.aux
    for (j2 in OP.ID.aux[-j]) {
      matmf.aux1OP[,j2] <- factor(matmf.aux[,j2], ordered=FALSE)
    }
    matmf.aux1OP <- model.frame(formula, drop.unused.levels = TRUE, data = matmf.aux1OP)
    attr(attr(matmf.aux1OP,"terms"),"dataClasses")
    resMDCPaux <- mdcp(formula=formula, tryAllMonoDir = TRUE, data=matmf.aux1OP,
                       method="MDCS3")
    resConfRegTest[j,c("UMLE_logLik","CMLE_logLik")] <- c(resMDCPaux$UMLE_logLik,resMDCPaux$log.lik)
    resConfRegTest[j,"degreesOfFreedom"] <- resConfRegTest[j,"Num_Cat"]-1
    resConfRegTest[j,"Statistic"] <- 2*(resConfRegTest[j,"UMLE_logLik"]-
                                          resConfRegTest[j,"CMLE_logLik"])
    resConfRegTest[j,"CritValue"] <- qchisq(1-SignifLevel,df = resConfRegTest[j,"degreesOfFreedom"])
    resConfRegTest[j,"SignifLevel"] <- SignifLevel
    resConfRegTest[j,"P.Value"] <- 1- pchisq(resConfRegTest[j,"Statistic"],
                                             resConfRegTest[j,"degreesOfFreedom"], lower.tail = TRUE)
    resConfRegTest[j,"RejectMonotonicity"] <- resConfRegTest[j,"Statistic"]>
      resConfRegTest[j,"CritValue"]
  }
  resConfRegTest
  list(resConfRegTest = resConfRegTest)
}
