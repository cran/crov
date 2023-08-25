
#' @title Parameter Vector in Confidence Regions UCR and/or UCCR
#'
#' @description Determines whether a parameter vector is in the confidence region UCR and/or UCCR, according to the definitions in Espinosa and Hennig (2023) <https://doi.org/10.48550/arXiv.2107.04946>.
#' @param UMLE A vector with the unconstrained maximum likelihood estimates.
#' @param paramVals A vector with the parameter values for which it is needed to
#' assess whether it is part of one of the confidence regions or not.
#' The order of the parameters must be the same as the one of \code{UMLE}.
#' As in Espinosa and Hennig (2023), the parameter vector contains the parameters
#' of interest, beta_{0r}, and the remaining ones are the unconstrained MLEs given beta_{0r}.
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
#' \code{UMLE_logLik}=log-likelihood of the unconstrained model,
#' \code{param_logLik}=log-likelihood of the model using \code{paramVals},
#' \code{monotonicBeta0}=logical value, \code{TRUE} if the set of parameters
#' of \code{paramVals} indicated by \code{paramIDs} are monotonic,
#' \code{df}=degrees of freedom used to calculate the critical value,
#' \code{StatUCR}=value of the statistic used for \code{UCR},
#' \code{StatUCCR}=value of the statistic used for \code{UCCR},
#' \code{CritValue}=critical value, chi-squared with \code{df} and \code{1-SignifLevel},
#' \code{SignifLevel}=significance level used to calculate the critical value,
#' \code{inUCR}=logical value, \code{TRUE} if \code{paramVals} belongs to the confidence region \code{UCR},
#' \code{inUCCR}=logical value, \code{TRUE} if \code{paramVals} belongs to the confidence region \code{UCCR},
#' @examples resAux <- mdcp(QoL ~ EduLevel + Age + IncQuint + Gender + Health, data = crovData)
#' plotCMLE(resAux)
#' myVector <- resAux$estimates
#' myVectorID <- 10:12
#' myVector[myVectorID]
#'
#' # non-monotonic beta_{0r}, paramVals in UCR but not in UCCR
#' myVector[myVectorID] <- seq(0.195,0.185,length.out=3)
#' confRegUCRandUCCR(UMLE=resAux$UMLE, paramVals=myVector, paramIDs=myVectorID,SignifLevel=0.05, df=3,
#' matY= resAux$matY, matX= resAux$matX)
#'
#' # monotonic beta_{0r}, paramVals in UCR and UCCR
#' myVector[myVectorID] <- seq(0.073,0.074,length.out=3)
#' confRegUCRandUCCR(UMLE=resAux$UMLE, paramVals=myVector, paramIDs=myVectorID,SignifLevel=0.05, df=3,
#' matY= resAux$matY, matX= resAux$matX)
#'
#' # monotonic beta_{0r}, paramVals out of UCR and UCCR
#' myVector[myVectorID] <- seq(0.072,0.073,length.out=3)
#' confRegUCRandUCCR(UMLE=resAux$UMLE, paramVals=myVector, paramIDs=myVectorID,SignifLevel=0.05, df=3,
#' matY= resAux$matY, matX= resAux$matX)
#' @seealso \code{\link[crov:confRegCCR]{confRegCCR}},
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

confRegUCRandUCCR <- function(UMLE=NULL, paramVals=NULL, paramIDs=NULL,
                              SignifLevel=0.05, df, matY, matX) {

  firstDifParams <- paramVals[paramIDs]-c(0,paramVals[paramIDs][-length(paramIDs)])
  monotonic <- (sum(firstDifParams>=0)==length(paramIDs) |
                  sum(firstDifParams<=0)==length(paramIDs)) &
    (sum(paramVals[paramIDs]>=0)==length(paramIDs) | sum(paramVals[paramIDs]<=0)==length(paramIDs))

  confRegions <- data.frame("UMLE_logLik"=0,"param_logLik"=0,"monotonicBeta0"=FALSE,
                            "df"=0, "StatUCR"=0, "StatUCCR"=0,
                            "CritValue"=0, "SignifLevel"=0, "inUCR"=FALSE, "inUCCR"=FALSE)

  confRegions$UMLE_logLik    <- loglikOrdPred(paramInit=UMLE,matY=matY,matX=matX)
  confRegions$param_logLik   <- loglikOrdPred(paramInit=paramVals,matY=matY,matX=matX)
  confRegions$monotonicBeta0 <- monotonic
  confRegions$df             <- df
  confRegions$StatUCR       <- 2*(confRegions$UMLE_logLik-confRegions$param_logLik)
  confRegions$StatUCCR       <-ifelse(monotonic==FALSE,NA,
                                      2*(confRegions$UMLE_logLik-confRegions$param_logLik))
  confRegions$CritValue      <- qchisq(1-SignifLevel,df = df)
  confRegions$SignifLevel    <- SignifLevel
  confRegions$inUCR          <- confRegions$StatUCR  <= confRegions$CritValue
  confRegions$inUCCR         <- ifelse(monotonic==FALSE,FALSE,
                                       confRegions$StatUCCR<=confRegions$CritValue)

  list(confRegions)

}
