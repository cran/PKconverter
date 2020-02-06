#' Convert pharmacokinetic parameters for one compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distribution(V1) and
#' parameter (alpha) in the model "Aexp(-alpha)"
#' @usage OneComp_Volume_Exponent(V1,alpha,V1.sd=NA,alpha.sd=NA,
#'          covar=c(V1alpha=NA),...)
#' @param V1 The volume of distribution of compartment 1
#' @param alpha parameter in one compartment model "Aexp(-alpha)"
#' @param V1.sd standard error of V1
#' @param alpha.sd standard error of A
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' OneComp_Volume_Exponent(V1=8,alpha=0.5,V1.sd=0.01,alpha.sd=0.001)
#'
OneComp_Volume_Exponent<-function(V1,alpha,V1.sd=NA,alpha.sd=NA,
  covar=c(V1alpha=NA),...){
  if(is.na(covar[1])) covar<-0
  V1.var = (V1.sd)^2
  alpha.var = (alpha.sd)^2

  Vdss<-V1
  Vdss.sd<-V1.sd

  k10<-alpha
  k10.sd<-alpha.sd

  Cl1<-V1*alpha
  sigma2<-matrix(as.numeric(c(V1.var,covar[1],covar[1],alpha.var)),2,2,byrow=T)
  Cl1_deriv<-as.matrix(attr(eval(stats::deriv(~V1*alpha,c("V1","alpha"))),
                            "gradient"))
  Cl1.sd<-sqrt(Cl1_deriv %*% sigma2 %*% t(Cl1_deriv))

  True_A<-1/V1
  True_A_deriv<-as.matrix(attr(eval(stats::deriv(~1/V1,"V1")),"gradient"))
  True_A.sd<-sqrt(True_A_deriv * V1.var *True_A_deriv)

  Frac_A<-1
  Frac_A.sd<-ifelse(is.na(V1.sd),NA,0)

  t_alpha<-log(2)/alpha
  t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/alpha,"alpha")),
                                "gradient"))
  t_alpha.sd<-sqrt(t_alpha_deriv *alpha.var * t_alpha_deriv)

  if(is.na(V1[1])){
    param = rep(NA,8)
    sd = rep(NA,8)
  }else{
    param =c(V1,alpha,Cl1,k10,t_alpha,Vdss,True_A,Frac_A)
    sd = c(V1.sd,alpha.sd,Cl1.sd,k10.sd,t_alpha.sd,Vdss.sd,True_A.sd,Frac_A.sd)
  }

  result = data.frame(Parameter=c("V1","alpha","Cl1","k10","t_alpha",
                      "Vdss","True_A","Frac_A"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("V1","alpha","Cl1","k10","t_alpha","Vdss",
                         "True_A","Frac_A")
  result<-result[c("Vdss","V1","Cl1","k10","alpha",
                   "t_alpha","True_A","Frac_A"),]

  return(result)
}
