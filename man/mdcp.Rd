% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mdcp.R
\name{mdcp}
\alias{mdcp}
\title{Monotonicity Direction Classification (MDC) procedure}
\usage{
mdcp(
  formula,
  data = NULL,
  tryAllMonoDir = FALSE,
  monoDir = NULL,
  CLS1 = 0.95,
  TLBS2 = 0.85,
  TLNS2 = 0.999,
  StepSizeCLS2 = 1e-04,
  method = NULL,
  monoTestSignLevel = 0.05
)
}
\arguments{
\item{formula}{A \code{formula} to be fitted with ordinal response, one or more ordinal predictors, and possibly one or more other predictors.
For ordinal response and ordinal predictors use ordered factors.}

\item{data}{A data.frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data.frame), containing the
variables in \code{formula}. Neither a matrix nor an array will be accepted.}

\item{tryAllMonoDir}{A logical value that indicates whether one model should be fitted for each one of the possible combinations of monotonicity
directions. Use \code{TRUE} if none monotonicity direction is pre-specified using \code{monoDir} and the MDC procedure is not used.}

\item{monoDir}{Vector with monotonicity directions for the ordinal predictors to be used as constraints. Possible values for \code{monoDir} are
\code{TRUE} and \code{FALSE}. Use \code{TRUE} for "isotonic" and \code{FALSE} for "antitonic". The order of the elements in \code{monoDir} must be
the same as the order of the ordinal predictors in the object \code{formula}, i.e., the j-th element of \code{monoDir} must correspond to
the monotonicity direction of the j-th ordinal predictor in \code{formula}. If \code{tryAllMonoDir} and \code{monoDir} are not used (default option),
the monotonicity direction classification prodecure is executed to find the monotonicity directions associated to the model with the maximum log-likelihood.}

\item{CLS1}{Numerical value for the confidence level to be used in the first step of the MDC procedure. This parameter is active if
\code{tryAllMonoDir} and \code{monoDir} are not used.}

\item{TLBS2}{Numerical value for the tolerance level to be used in the second step of the MDC procedure over those ordinal
predictors classified as "Both" in the first step. This parameter is active if \code{tryAllMonoDir} and \code{monoDir} are not used.}

\item{TLNS2}{Numerical value for the tolerance level to be used in the second step of the MDC procedure over those ordinal
predictors classified as "None" in the first step. This parameter is active if \code{tryAllMonoDir} and \code{monoDir} are not used.}

\item{StepSizeCLS2}{Numerical value for the magnitude in which the confidence levels will be increased or decreased during the second step of
the MDC procedure. This parameter is active if \code{tryAllMonoDir} and \code{monoDir} are not used.}

\item{method}{The type of constrained method to be used among "MDCS1", "MDCS2", "MDCS3", "CMLEbonferroni", "CMLEconfReg", and "CMLEfiltered". Default value corresponds to "MDCS3".}

\item{monoTestSignLevel}{Significance level used when method is "CMLEbonferroni" or "CMLEconfReg". Default value 0.05.}
}
\value{
\code{MDCproc}: Data frame with the monotonicity direction classification (Isotonic, Antitonic, Both, or None) used for each
ordinal predictor in each one of the steps of the MDC procedure (S1, S2 and S3), together with their individual confidence levels (CL). If
\code{monoDir} is used, \code{MDCproc} shows the monotonicity directions in \code{monoDir}.

\code{estimates}: Vector of parameter estimates of the model.

\code{estimates_se}: Vector of standard errors of the parameter estimates of the model.

\code{log.lik}: Value of the log-likelihood of the model.

\code{allModels}: Data frame with monotonicity directions, log-likelihood and parameter estimates of all models involved in the third step of the MDC
procedure. If parameter \code{monoDir} is used, \code{allModels} shows these results from the model with monotonicity directions
used in \code{monoDir} only. If parameter \code{tryAllMonoDir} is used, \code{allModels} shows these results from all the models according to all possible
combinations of monotonicity directions.

\code{constrOptimRes}: List with the outcomes provided by the function \code{\link[stats:constrOptim]{constrOptim}}.

\code{UMLE}: Vector with the parameter estimates of the unconstrained version of the model.

\code{UMLE_SE}: Vector with the standard errors of the unconstrained version of the model.
}
\description{
Fits a constrained regression model for an ordinal response with ordinal predictors
and possibly others, Espinosa and Hennig (2019) <https://doi.org/10.1007/s11222-018-9842-2>. The parameter estimates associated with an ordinal
predictor are constrained to be monotonic. If a monotonicity direction (isotonic or antitonic) is
not specified for an ordinal predictor by the user, then the monotonicity direction classification
procedure establishes it.
}
\examples{
# Ordinal predictors: EduLevel, IncQuint, Health,
# Overcrowd, and NumChildren
# mdcpExample <- mdcp(QoL ~ EduLevel + Age + IncQuint + Gender + Health +
# Overcrowd + Activity + NumChildren, data = crovData,
# CLS1 = 0.95, TLBS2 = 0.90, TLNS2 = 0.99, StepSizeCLS2 = 0.0002)
# mdcpExample$MDCproc
# cbind("CMLE"=mdcpExample$estimates,"UMLE"=mdcpExample$UMLE)
# mdcpExample$UMLE_SE
# mdcpExample$log.lik
# mdcpExample$allModels[1:6]
}
\references{
Espinosa, J., Hennig, C. A constrained regression model for an ordinal response
with ordinal predictors. Stat Comput 29, 869-890 (2019). https://doi.org/10.1007/s11222-018-9842-2.
}
\seealso{
\code{\link[crov:monoTestBonf]{monoTestBonf}}, \code{\link[crov:monoTestConfReg]{monoTestConfReg}}, \code{\link[stats:constrOptim]{constrOptim}}.
}
