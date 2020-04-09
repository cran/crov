cmleConfReg <- function(formula, data = NULL, signLevelConfReg = 0.95 ) {

  resMonoTest.aux <- monoTestConfReg(formula=formula, data = data,
                                     SignifLevel=signLevelConfReg)

  namesOPs <-
    row.names(resMonoTest.aux$resConfRegTest)[
      which(resMonoTest.aux$resConfRegTest[,"RejectMonotonicity"]==FALSE)]

  newData <- data

  if(length(namesOPs)==0) {toBeUnordered <- row.names(resMonoTest.aux$resConfRegTest)} else {
    toBeUnordered <- row.names(resMonoTest.aux$resConfRegTest)[-which(row.names(resMonoTest.aux$resConfRegTest)%in%namesOPs)]
  }

  for (i in toBeUnordered) {newData[,i] <- factor(newData[,i],ordered=FALSE)}

  return(list(cmleConfRegresults=resMonoTest.aux,
              namesOPs=namesOPs,
              newData=newData))

}
