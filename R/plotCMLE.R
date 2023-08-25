#' @title Plot unconstrained and constrained proportional odds logit model
#'
#' @description Uses the results of function \code{mdcp} to produce a plot for the Maximum Likelihood Estimators of the parameters of both
#' the unconstrained and constrained proportional odds logit models (UMLE and CMLE
#' correspondingly). The UMLE includes confidence intervals. Parameter estimates of ordinal predictors
#' are graphically linked with segments.
#' @param mdcpResult An object of class \code{list} storing the results of function \code{mdcp}, which fits both the unconstrained
#'  and constrained proportional odds logit models.
#' @param SignifLevel Significance level to be used when constructing the confidence intervals for each parameter
#' of the unconstrained proportional odds logit model. Default value 0.05.
#' @param xposLegend Position of legend on the x-axis. If \code{xposLegend} or \code{yposLegend} are
#' not used, then the legend is located using \code{topleft} option.
#' @param yposLegend Position of legend on the y-axis. If \code{xposLegend} or \code{yposLegend} are
#' not used, then the legend is located using \code{topleft} option.
#' @param xcex.axis Size of \code{cex.axis} for the x-axis. Default value is 0.8.
#' @param ycex.axis Size of \code{cex.axis} for the y-axis. Default value is 0.8.
#' @param cexLegend Size of legend text to be assigned to \code{cex}. Default value is 1.
#' @param methodName Method name to be used in the main title of the plot.
#' @return Plot.
#' @examples # Ordinal predictors: EduLevel, IncQuint, Health,
#' # Overcrowd, and NumChildren
#' mdcpExample <- mdcp(QoL ~ EduLevel + Age + IncQuint + Gender + Health, data = crovData,
#' CLS1 = 0.95, TLBS2 = 0.90, TLNS2 = 0.99, StepSizeCLS2 = 0.0002)
#' plotCMLE(mdcpResult=mdcpExample,SignifLevel=0.05,xposLegend=14, yposLegend=4.8,
#' cexLegend=0.8, method="MDCS3")
#' @seealso \code{\link[crov:monoTestBonf]{monoTestBonf}},
#' \code{\link[crov:monoTestConfReg]{monoTestConfReg}},
#' \code{\link[crov:monoTestBonf]{monoTestBonf}},
#' \code{\link[stats:constrOptim]{constrOptim}}.
#' @importFrom stats qnorm
#' @importFrom graphics par abline segments axis box legend points
#' @export

plotCMLE <- function(mdcpResult=NULL, SignifLevel=0.05,
                     xposLegend=NULL, yposLegend=NULL,
                     xcex.axis=0.8, ycex.axis=0.8,
                     cexLegend=1, methodName="Not indicated") {
  LI <- mdcpResult$UMLE + qnorm(SignifLevel/2,0,1)*mdcpResult$UMLE_SE
  LS <- mdcpResult$UMLE + qnorm(1-SignifLevel/2,0,1)*mdcpResult$UMLE_SE
  alphas.ID <- 1:(dim(mdcpResult$matY)[2]-1)
  q_cat_OrdPred <- mdcpResult$q_cat_OrdPred
  matrixOP.ID <- matrix(0,nrow=length(mdcpResult$OPsNames),ncol=dim(mdcpResult$matX)[2])
  rowMat <- 0
  for (x in mdcpResult$OPsNames) {
    rowMat <- rowMat + 1
    matrixOP.ID[rowMat,] <- grepl(x,colnames(mdcpResult$matX))
  }
  id_from_OP <- NULL
  for (x in 1:length(mdcpResult$OPsNames)) {
    id_from_OP <- c(id_from_OP,match(1,t(apply(matrixOP.ID,1,cumsum))[x,]))
  }
  id_from_OP <- id_from_OP + length(alphas.ID)
  id_to_OP <- id_from_OP + (q_cat_OrdPred - 2)

  par(mar=c(5.6,3.8,1.5,0.5),no.readonly=TRUE)
  plot(mdcpResult$UMLE, axes=FALSE, pch=19, xlab="",ylab="", main=paste("Method:",methodName),
       ylim=range(c(LI,LS)), col="darkorange3")
  abline(h=0, lty=2, col="grey80")
  segments(1:length(mdcpResult$UMLE),LI,1:length(mdcpResult$UMLE),LS,col="darkorange3")
  for (i in 1:length(id_from_OP)) {
    for (j in id_from_OP[i]:(id_to_OP[i]-1)) {
      segments(j,mdcpResult$UMLE[j], j+1, mdcpResult$UMLE[j+1],
               col="darkorange3",  lty=1, lwd=1)
    }
  }
  axis(1,las=2,at=1:length(mdcpResult$UMLE),labels=names(mdcpResult$UMLE),cex.axis=xcex.axis)
  axis(2,las=1,at=round(seq(min(LI),max(LS),l=5),4),cex.axis=ycex.axis)
  points(mdcpResult$estimates, pch=13, col="darkblue")
  box()
  if (is.null(xposLegend) | is.null(yposLegend) ) {
    legend(x="topleft", legend=c(paste("UMLE and ", 100-SignifLevel*100,"% CIs",sep=""),"CMLE"),
           col=c("darkorange3","darkblue"),
           lwd=1, lty=c(1,1), pch=c(19,13),cex=cexLegend,bty='n',y.intersp = 2)
  } else {
    legend(x=xposLegend,y=yposLegend, legend=c(paste("UMLE and ", 100-SignifLevel*100,"% CIs",sep=""),"CMLE"),
           col=c("darkorange3","darkblue"),
           lwd=1, lty=c(1,1), pch=c(19,13),cex=cexLegend,bty='n',y.intersp = 2)
  }
}
