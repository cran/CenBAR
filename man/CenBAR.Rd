\name{CenBAR}
\alias{CenBAR}
\title{Broken Adaptive Ridge Estimator for Censored Data in AFT Model}
\usage{
CenBAR(X,Y,delta,lambda.path=NULL, enableScreening=FALSE)
}
\arguments{
  \item{X}{input matrix, of dimension nobs x nvars; each row is an observation vector.}
  \item{Y}{response variable.}
  \item{delta}{The status indicator, normally 0=alive, 1=dead.}
  \item{lambda.path}{A user supplied lambda sequence. One usage is to have the program
         compute its own lambda sequence based on nlambda and lambdaMax.
         lamdMax = max((t(x)*Y)^2/(4*t(x)*x)). The other usage is use the sequence depend
         on user's data.}
  \item{enableScreening}{If nobs > nvars, there is no need to do screening;
        If nobs <= nvars, it will do variable screening and then variable selection
        and estimate (defalt is FALSE).}
}
\value{
  \item{beta}{the coefficients estimation of the variables.}
}
\description{
Prints 'Broken adaptive ridge (BAR) method to the semi-parametric
     accelerated failure time (AFT) model for right-censored survival data
     by applying the Leurgan's synthetic data.'.
}
\examples{
  X=matrix(rnorm(10*2),10,2)
  Y=abs(rnorm(10))
  delta=sample(0:1,10,replace=TRUE)
  lambda.path <- seq(0.1, 10, l=5)
  fit = CenBAR(X,Y,delta,lambda.path)
}
\author{
  Zhihua Sun, Chunyu Yu
}
