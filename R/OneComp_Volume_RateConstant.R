#' Convert pharmacokinetic parameters for one compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distribution (V1) and
#' elimination rate constant (k10)
#' @usage OneComp_Volume_RateConstant(V1,k10,
#'              V1.sd=NA,k10.sd=NA,covar=c(V1k10=0),...)
#' @param V1 The volume of distribution of compartment 1
#' @param k10 elimination rate constant
#' @param V1.sd standard error of V1
#' @param k10.sd standard error of k10
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' OneComp_Volume_RateConstant(V1=8,k10=0.5,V1.sd=0.01,k10.sd=0.002)
#'
OneComp_Volume_RateConstant<-function(V1,k10,
    V1.sd=NA,k10.sd=NA,covar=c(V1k10=0),...){
  if(is.na(covar[1])) covar<-0
  V1.var = (V1.sd)^2
  k10.var = (k10.sd)^2

  Vdss<-V1
  Vdss.sd<-V1.sd

  Cl1<-V1*k10
  sigma<-matrix(as.numeric(c(V1.var,covar[1],covar[1],k10.var)),2,2,byrow=T)
  Cl1_deriv<-as.matrix(attr(eval(stats::deriv(~V1*k10,c("V1","k10"))),
                            "gradient"))
  Cl1.sd<-sqrt(Cl1_deriv %*% sigma %*% t(Cl1_deriv))

  t_alpha<-log(2)/k10
  t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/k10,"k10")),
                                "gradient"))
  t_alpha.sd<-sqrt(t_alpha_deriv * k10.var *t_alpha_deriv)

  true_A<-1/V1
  true_A_deriv<-attr(eval(stats::deriv(~1/V1,"V1")),"gradient")
  true_A.sd<-sqrt(true_A_deriv* V1.var * true_A_deriv)

  frac_A<-1
  frac_A.sd<-ifelse(is.na(V1.sd),NA,0)

  alpha<-k10
  alpha.sd<-k10.sd
  if(is.na(V1[1])){
    param = rep(NA,8)
    sd = rep(NA,8)
  } else{
    param =  c(V1,k10,Cl1,t_alpha,Vdss,true_A,frac_A,alpha)
    sd = c(V1.sd,k10.sd,Cl1.sd,t_alpha.sd,Vdss.sd,true_A.sd,frac_A.sd,alpha.sd)
  }

  result = data.frame(Parameter=c("V1","k10","Cl1","t_alpha","Vdss",
                      "True_A","Frac_A","alpha"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("V1","k10","Cl1","t_alpha","Vdss",
                         "True_A","Frac_A","alpha")
  result<-result[c("Vdss","V1","Cl1","k10","alpha","t_alpha",
                 "True_A","Frac_A"),]
  return(result)
}
