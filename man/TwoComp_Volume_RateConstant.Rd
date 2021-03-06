% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TwoComp_Volume_RateConstant.R
\name{TwoComp_Volume_RateConstant}
\alias{TwoComp_Volume_RateConstant}
\title{Convert pharmacokinetic parameters for two compartment model}
\usage{
TwoComp_Volume_RateConstant(V1,k10,k12,k21,
              V1.sd=NA,k10.sd=NA,k12.sd=NA,k21.sd=NA,
              covar=c(V1k10=NA,V1k12=NA,V1k21=NA,
              k10k12=NA,k10k21=NA,k12k21=NA),...)
}
\arguments{
\item{V1}{The volume of distribution of compartment 1}

\item{k10}{elimination rate constant}

\item{k12}{transfer rate constants from compartment 1 to compartment 2}

\item{k21}{transfer rate constants from compartment 2 to compartment 1}

\item{V1.sd}{standard error of V1}

\item{k10.sd}{standard error of k10}

\item{k12.sd}{standard error of k12}

\item{k21.sd}{standard error of k21}

\item{covar}{covariances among parameters}

\item{...}{arguments to be passed to methods}
}
\description{
Calculate pharmacokinetic parameters with volume of distribution (V1),
elimination rate constant (k10), and transter rate constants (k12, k21)
}
\examples{
TwoComp_Volume_RateConstant(V1=5,k10=0.7,k12=0.5,k21=0.05,
         V1.sd=0.01,k10.sd=0.002,k12.sd=0.001,k21.sd=0.0005)
}
\references{
\url{http://www.nonmemcourse.com/convert.xls}
}
