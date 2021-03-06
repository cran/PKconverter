% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TwoComp_Volume_Clearance_HalfLife.R
\name{TwoComp_Volume_Clearance_HalfLife}
\alias{TwoComp_Volume_Clearance_HalfLife}
\title{Convert pharmacokinetic parameters for two compartment model}
\usage{
TwoComp_Volume_Clearance_HalfLife(V1,Cl1,t_alpha,t_beta,
 V1.sd=NA,Cl1.sd=NA,t_alpha.sd=NA,
 t_beta.sd=NA,covar=c(V1Cl1=NA,V1talpha=NA,V1tbeta=NA,Cl1talpha=NA,
   Cl1tbeta=NA,talphatbeta=NA),...)
}
\arguments{
\item{V1}{The volume of distribution of compartment 1}

\item{Cl1}{Clearance from compartment 1}

\item{t_alpha}{half life of compartment 1}

\item{t_beta}{half life of compartment 2}

\item{V1.sd}{standard error of V1}

\item{Cl1.sd}{standard error of Cl1}

\item{t_alpha.sd}{standard error of t_alpha}

\item{t_beta.sd}{standard error of t_beta}

\item{covar}{covariances among parameters}

\item{...}{arguments to be passed to methods}
}
\description{
Calculate pharmacokinetic parameters with volume of distribution(V1),
clearance (Cl1) and half-lives (t_alpha and t_beta)
}
\examples{
TwoComp_Volume_Clearance_HalfLife(V1=5,Cl1=3.5,t_alpha=0.568,t_beta=24.2,
V1.sd=0.01,Cl1.sd=0.01,t_alpha.sd=0.002,t_beta.sd=0.5)
}
\references{
\url{http://www.nonmemcourse.com/convert.xls}
}
