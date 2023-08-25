#' @title Parameter Vector in Confidence Region CCR
#'
#' @description Determines whether a parameter vector is in the confidence region CCR, according to the definitions in Espinosa and Hennig (2023) <https://doi.org/10.48550/arXiv.2107.04946>.
#' @param CMLE A vector with the constrained maximum likelihood estimates.
#' @param paramVals A vector with the parameter values for which it is needed to
#' assess whether it is part of the confidence region or not.
#' The order of the parameters must be the same as the one of \code{CMLE}.
#' As in Espinosa and Hennig (2023), the parameter vector contains the parameters
#' of interest, beta_{0r}, and the remaining ones are the constrained MLEs given beta_{0r}.
#' @param paramIDs A vector indicating the positions of the parameter values of
#' beta_{0r} in \code{paramVals}, which are those of interest, usually the ones
#' corresponding to some ordinal predictor. For instance, \code{paramIDs=7:11} indicates
#' that the 7th to the 11th parameter values in \code{paramVals} are the ones of interest
#' and correspond to some ordinal predictor.
#' @param SignifLevel A decimal number indicating the significant level. Usually, 0.05.
#' @param df Degrees of freedom to be used.
#' @param matY matY resulting from mdcp().
#' @param matX matX resulting from mdcp().
#' @return \code{confRegions}: Data frame with columns:
#' \code{CMLE_logLik}=log-likelihood of the constrained model,
#' \code{param_logLik}=log-likelihood of the model using \code{paramVals},
#' \code{monotonicBeta0}=logical value, \code{TRUE} if the set of parameters
#' of \code{paramVals} indicated by \code{paramIDs} are monotonic,
#' \code{df}=degrees of freedom used to calculate the critical value,
#' \code{StatCCR}=value of the statistic used for \code{CCR},
#' \code{CritValue}=critical value, chi-squared with \code{df} and \code{1-SignifLevel},
#' \code{SignifLevel}=significance level used to calculate the critical value,
#' \code{inCCR}=logical value, \code{TRUE} if \code{paramVals} belongs to the confidence region \code{CCR},
#' @examples resAux <- mdcp(QoL ~ EduLevel + Age + IncQuint + Gender + Health, data = crovData)
#' plotCMLE(resAux)
#' myVector <- resAux$estimates
#' myVectorID <- 10:12
#' myVector[myVectorID]
#'
#' # non-monotonic beta_{0r}
#' myVector[myVectorID] <- seq(0.195,0.185,length.out=3)
#' confRegCCR(CMLE=resAux$estimates, paramVals=myVector, paramIDs=myVectorID,SignifLevel=0.05, df=3,
#' matY= resAux$matY, matX= resAux$matX)
#'
#' # monotonic beta_{0r} and paramVals in CCR
#' myVector[myVectorID] <- seq(0.048,0.049,length.out=3)
#' confRegCCR(CMLE=resAux$estimates, paramVals=myVector, paramIDs=myVectorID,SignifLevel=0.05, df=3,
#' matY= resAux$matY, matX= resAux$matX)
#'
#' # monotonic beta_{0r} and paramVals out of CCR
#' myVector[myVectorID] <- seq(0.047,0.048,length.out=3)
#' confRegCCR(CMLE=resAux$estimates, paramVals=myVector, paramIDs=myVectorID,SignifLevel=0.05, df=3,
#' matY= resAux$matY, matX= resAux$matX)
#' @seealso \code{\link[crov:confRegUCRandUCCR]{confRegUCRandUCCR}},
#' \code{\link[crov:mdcp]{mdcp}},
#' \code{\link[crov:monoTestBonf]{monoTestBonf}},
#' \code{\link[crov:monoTestConfReg]{monoTestConfReg}},
#' \code{\link[crov:plotCMLE]{plotCMLE}},
#' \code{\link[VGAM:vglm]{vlgm}}.
#' @importFrom stats qnorm qchisq pchisq
#' @references Espinosa, J., and Hennig, C. "Inference for the proportional odds cumulative logit model with monotonicity constraints for ordinal predictors and ordinal response." Arxiv (2023).
#' <https://doi.org/10.48550/arXiv.2107.04946>.
#' @export
#'

confRegCCR <- function(CMLE=NULL, paramVals=NULL, paramIDs=NULL,
                       SignifLevel=0.05, df, matY, matX) {

  firstDifParams <- paramVals[paramIDs]-c(0,paramVals[paramIDs][-length(paramIDs)])
  monotonic <- (sum(firstDifParams>=0)==length(paramIDs) |
                  sum(firstDifParams<=0)==length(paramIDs)) &
    (sum(paramVals[paramIDs]>=0)==length(paramIDs) | sum(paramVals[paramIDs]<=0)==length(paramIDs))

  confRegions <- data.frame("CMLE_logLik"=0,"param_logLik"=0,"monotonicBeta0"=FALSE,
                            "df"=0, "StatCCR"=0,
                            "CritValue"=0, "SignifLevel"=0, "inCCR"=FALSE)

  confRegions$CMLE_logLik    <- loglikOrdPred(paramInit=CMLE,matY=matY,matX=matX)
  confRegions$param_logLik   <- loglikOrdPred(paramInit=paramVals,matY=matY,matX=matX)
  confRegions$monotonicBeta0 <- monotonic
  confRegions$df             <- df
  confRegions$StatCCR        <-ifelse(monotonic==FALSE,NA,
                                      2*(confRegions$CMLE_logLik-confRegions$param_logLik))
  confRegions$CritValue      <- qchisq(1-SignifLevel,df = df)
  confRegions$SignifLevel    <- SignifLevel
  confRegions$inCCR          <- ifelse(monotonic==FALSE,FALSE,
                                       confRegions$StatCCR<=confRegions$CritValue)

  list(confRegions)

}
