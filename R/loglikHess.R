# @title Log-likelihood Hessian matrix
#
# @description Computes the Hessian matrix of the log-likelihood to calculate the SE once the parameter estimates are obtained.
# @param paramInit Vector of parameter values for intitialisation.
# @param matY Matrix with binary values of the ordinal response variable.
# @param matX Matrix with binary values of the ordinal predictors variables and non-ordinal predictors.
# @return  Hessian matrix of the Log-likelihood.
loglikHess<- function(paramInit,matY,matX) {
  N<-dim(matY)[1]
  c<-dim(matY)[2]
  qbTot<-dim(matX)[2]
  paramIter<-c(-999999999,paramInit[1:(c-1)],999999999,paramInit[c:(c-1+qbTot)])

  z<-as.matrix(rep(1,N))%*%paramIter[1:(c+1)]+matX%*%paramIter[(c+2):(length(paramIter))]%*%t(as.matrix(rep(1,c+1)))
  z_j<-z[,2:(c+1)]
  z_jm1<-z[,1:c]
  CF_j<-ifelse(is.finite(exp(z_j)/(1+exp(z_j))), exp(z_j)/(1+exp(z_j)), 1)
  CF_jm1<-exp(z_jm1)/(1+exp(z_jm1))
  DF_j<-ifelse(is.finite(exp(z_j)/((1+exp(z_j))^2)), exp(z_j)/((1+exp(z_j))^2), 0)
  DF_jm1<-exp(z_jm1)/((1+exp(z_jm1))^2)
  zOwn_j<-ifelse(is.finite((exp(z_j)-1)/((1+exp(z_j)))), (exp(z_j)-1)/((1+exp(z_j))), 0)
  zOwn_jm1<-(exp(z_jm1)-1)/((1+exp(z_jm1)))

  HessianMat<-matrix(0,ncol=c-1+qbTot,nrow=c-1+qbTot)

  #### Computation of second partial derivatives of alpha_k,alpha_l
  k<-1; l<-1
  for (k in 1:(c-1)){
    for (l in 1:(c-1)) {
      kron_jk<-matrix(0,ncol=c,nrow=N)
      kron_jm1k<-matrix(0,ncol=c,nrow=N)
      kron_jk[,k]<-1
      kron_jm1k[,k+1]<-1

      kron_jl<-matrix(0,ncol=c,nrow=N)
      kron_jm1l<-matrix(0,ncol=c,nrow=N)
      kron_jl[,l]<-1
      kron_jm1l[,l+1]<-1

      SecParDer_a_k_a_l<-(CF_j-CF_jm1)*
        ((CF_j-CF_jm1)*(DF_jm1*zOwn_jm1*kron_jm1k*kron_jm1l-DF_j*zOwn_j*kron_jk*kron_jl)/((CF_j-CF_jm1)^2)
         -(DF_j*kron_jk-DF_jm1*kron_jm1k)*(DF_j*kron_jl-DF_jm1*kron_jm1l)/((CF_j-CF_jm1)^2))
      HessianMat[k,l]<-sum(sum(SecParDer_a_k_a_l))
      HessianMat[l,k]<-sum(sum(SecParDer_a_k_a_l))
    }
  }
  #HessianMat

  #### Computation of second partial derivatives of beta_k,alpha_l
  k<-1; l<-1
  for (k in 1:qbTot){
    for (l in 1:(c-1)) {
      kron_jl<-matrix(0,ncol=c,nrow=N)
      kron_jm1l<-matrix(0,ncol=c,nrow=N)
      kron_jl[,l]<-1
      kron_jm1l[,l+1]<-1

      SecParDer_b_k_a_l<-(CF_j-CF_jm1)*(matX[,k]%*%t(as.matrix(rep(1,c))))*
        ((DF_jm1-DF_j)*(DF_j*kron_jl-DF_jm1*kron_jm1l)/((CF_j-CF_jm1)^2)
         -(CF_j-CF_jm1)*(DF_j*zOwn_j*kron_jl-DF_jm1*zOwn_jm1*kron_jm1l)/((CF_j-CF_jm1)^2)
        )
      HessianMat[l,k+c-1]<-sum(sum(SecParDer_b_k_a_l))
      HessianMat[k+c-1,l]<-sum(sum(SecParDer_b_k_a_l))
    }
  }
  #HessianMat

  #### Computation of second partial derivatives of beta_k,beta_l
  k<-1; l<-1
  for (k in 1:qbTot){
    for (l in 1:qbTot) {
      SecParDer_b_k_b_l<-(CF_j-CF_jm1)*(matX[,k]%*%t(as.matrix(rep(1,c))))*(matX[,l]%*%t(as.matrix(rep(1,c))))*
        ((CF_j-CF_jm1)*(DF_jm1*zOwn_jm1-DF_j*zOwn_j)/((CF_j-CF_jm1)^2)
         -((DF_jm1-DF_j)^2)/((CF_j-CF_jm1)^2)
        )

      HessianMat[k+c-1,l+c-1]<-(sum(sum(SecParDer_b_k_b_l)))
      HessianMat[l+c-1,k+c-1]<-(sum(sum(SecParDer_b_k_b_l)))
    }
  }
  #HessianMat

  return(HessianMat)
}

