#' @title Monotonicity test
#'
#' @description Tests the null hypothesis of monotonicity over a set of parameters associated to an ordinal predictor, according to Espinosa and Hennig (2019) <DOI:10.1007/s11222-018-9842-2>.
#' @param simultAlpha Numerical value for the simultaneous significance level.
#' @param OP_UMLE Vector with the unconstrained parameter estimates of an ordinal predictor's categories represented by dummy variables
#' in an unconstrained model for ordinal response (see \code{\link[VGAM:vglm]{vlgm}}).
#' @param OP_SE Vector with the standard error of the parameters of an ordinal predictor's categories represented by dummy variables
#' in an unconstrained model for ordinal response (see \code{\link[VGAM:vglm]{vlgm}}).
#' @return \code{testRes}: String value with outcomes either "Reject H_0" or "Not Reject H_0".
#' @return \code{simultAlpha}: Numerical value with the simultaneous significance level.
#' @return \code{indivAlphaA}: Numerical value with the individual significance level for each confidence interval.
#' @return \code{simultPvalue}: Numerical value with the p-value associated to the simultaneous significance level.
#' @examples monoTestBonf(simultAlpha=0.05, OP_UMLE = c(-0.352177095,-0.403928770,
#' -0.290875028,-0.769834449), OP_SE = c(0.246638339,0.247723681,0.267577633,0.300951441))
#' @seealso \code{\link[crov:mdcp]{mdcp}},
#' \code{\link[crov:monoTestConfReg]{monoTestConfReg}},
#' \code{\link[crov:plotCMLE]{plotCMLE}},
#' \code{\link[VGAM:vglm]{vlgm}}.
#' @importFrom stats qnorm
#' @references Espinosa, J., and Hennig, C. "A constrained regression model for an ordinal response
#' with ordinal predictors." Statistics and Computing 29.5 (2019): 869-890.
#' https://doi.org/10.1007/s11222-018-9842-2.
#' @export
#'
monoTestBonf <- function(simultAlpha = 0.05, OP_UMLE, OP_SE){

  aux.res <- monoTestBonfSingle(simultAlpha = simultAlpha, OP_UMLE = OP_UMLE, OP_SE = OP_SE)

  alpLB <- 0
  alpUB <- 1
  while (alpUB - alpLB > 0.00000001) {
    auxTestRes <- monoTestBonfSingle(simultAlpha = alpLB + (alpUB - alpLB) / 2, OP_UMLE = OP_UMLE, OP_SE = OP_SE)
    if (auxTestRes$testRes == "Reject H_0") {alpUB <- alpLB + (alpUB - alpLB) / 2} else {alpLB <- alpLB + (alpUB - alpLB) / 2}
  }

  list(testRes = aux.res$testRes, simultAlpha = aux.res$simultAlpha, indivAlpha = aux.res$indivAlpha ,
       simultPvalue = round(ifelse(alpLB<0.5,alpLB,alpUB),7))
}
