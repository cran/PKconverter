#' Convert pharmacokinetic parameters for one compartment model
#'
#' Calculate pharmacokinetic parameters with parameters (A and alpha)
#' in one compartment model "Aexp(-alpha)"
#' @usage OneComp_Coefficient_Exponent(A,alpha,A.sd=NA,alpha.sd=NA,
#'            covar=c(Aalpha=NA),...)
#' @param A parameter in one compartment model "Aexp(-alpha)"
#' @param alpha parameter in one compartment model "Aexp(-alpha)"
#' @param A.sd standard error of A
#' @param alpha.sd standard error of alpha
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' OneComp_Coefficient_Exponent(A=0.125,alpha=0.5,A.sd=0.002,alpha.sd=0.009)

OneComp_Coefficient_Exponent<-function(A,alpha,A.sd=NA,
  alpha.sd=NA,covar=c(Aalpha=NA),...){

  if(is.na(covar[1])) covar<-0
  A.var = (A.sd)^2
  alpha.var = (alpha.sd)^2

  V1<-1/A
  Vdss<-V1

  V1_deriv<-as.matrix(attr(eval(stats::deriv(~1/A,"A")),"gradient"))
  V1.sd<-sqrt(V1_deriv * A.var * V1_deriv)
  Vdss.sd<-V1.sd

  Cl1<-alpha/A
  Cl1_deriv<-as.matrix(attr(eval(stats::deriv(~alpha/A,c("A","alpha"))),
                            "gradient"))
  sigma2<-matrix(as.numeric(c(A.var,covar[1],covar[1],alpha.var)),2,2,byrow=T)
  Cl1.sd<-sqrt(Cl1_deriv %*% sigma2 %*% t(Cl1_deriv))

  k10<-alpha
  k10.sd<-alpha.sd

  t_alpha<-log(2)/alpha
  t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/alpha,"alpha")),
                                "gradient"))
  t_alpha.sd<-sqrt(t_alpha_deriv * alpha.var *t_alpha_deriv)

  Frac_A<-1
  Frac_A.sd<-ifelse(is.na(A.sd),NA,0)
  if(is.na(A[1])){
    param = rep(NA,8)
    sd = rep(NA,8)
  } else{
    param = c(A,alpha,V1,Cl1,k10,t_alpha,Vdss,Frac_A)
    sd = c(A.sd,alpha.sd,V1.sd,Cl1.sd,k10.sd,t_alpha.sd,Vdss.sd,Frac_A.sd)
  }

  result = data.frame(Parameter=c("True_A","alpha","V1","Cl1","k10","t_alpha",
                       "Vdss","Frac_A"),
                       Estimate=param, Std.err=sd)
  row.names(result) <- c("True_A","alpha","V1","Cl1","k10","t_alpha","Vdss",
    "Frac_A")
  result<-result[c("Vdss","V1","Cl1","k10","alpha",
                   "t_alpha","True_A","Frac_A"),]
  return(result)
}
