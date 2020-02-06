#' Convert pharmacokinetic parameters for one compartment model
#'
#' Calculate pharmacokinetic parameters with clearance (Cl1) and
#' half-life (t_alpha)
#'
#' @usage OneComp_Volume_Clearance_HalfLife(Cl1,t_alpha,
#'          Cl1.sd=NA,t_alpha.sd=NA,covar=c(Cl1talpha=NA),...)
#' @param Cl1 Clearance from compartment 1
#' @param t_alpha half life of compartment 1
#' @param Cl1.sd standard error of Cl1
#' @param t_alpha.sd standard error of t_alpha
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' OneComp_Volume_Clearance_HalfLife(Cl1=4,t_alpha=0.568,
#'       Cl1.sd=0.01,t_alpha.sd=0.0003)

OneComp_Volume_Clearance_HalfLife<-function(Cl1,t_alpha,
  Cl1.sd=NA,t_alpha.sd=NA,covar=c(Cl1talpha=NA),...){
  if(is.na(covar[1])) covar<-0
  Cl1.var = (Cl1.sd)^2
  t_alpha.var = (t_alpha.sd)^2

  V1<-(Cl1*t_alpha)/log(2)
  Vdss<-V1

  sigma2<-matrix(as.numeric(c(Cl1.var,covar[1],covar[1],t_alpha.var)),
                 2,2,byrow=T)
  V1_deriv<-as.matrix(attr(eval(stats::deriv(~(Cl1*t_alpha)/log(2),
                                          c("Cl1","t_alpha"))),"gradient"))
  V1.sd<-sqrt(V1_deriv %*% sigma2 %*% t(V1_deriv))
  Vdss.sd<-V1.sd

  k10<-log(2)/t_alpha
  k10_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/t_alpha,"t_alpha")),
                              "gradient"))
  k10.sd<-sqrt(k10_deriv * t_alpha.var *k10_deriv)

  true_A<-log(2)/(Cl1*t_alpha)
  true_A_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/(Cl1*t_alpha),
                                         c("Cl1","t_alpha"))),"gradient"))
  true_A.sd<-sqrt(true_A_deriv %*% sigma2 %*% t(true_A_deriv))

  frac_A<-1
  frac_A.sd<-ifelse(is.na(Cl1.sd),NA,0)

  alpha<-k10
  alpha.sd<-k10.sd
  if(is.na(t_alpha[1])){
    param=rep(NA,8)
    sd=rep(NA,8)
  } else{
    param = c(Cl1,t_alpha,V1,k10,Vdss,true_A,frac_A,alpha)
    sd = c(Cl1.sd,t_alpha.sd,V1.sd,k10.sd,Vdss.sd,true_A.sd,frac_A.sd,alpha.sd)
  }

  result = data.frame(Parameter=c("Cl1","t_alpha","V1","k10","Vdss",
                      "True_A","Frac_A","alpha"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("Cl1","t_alpha","V1","k10","Vdss",
                         "True_A","Frac_A","alpha")
  result<-result[c("Vdss","V1","Cl1","k10","alpha",
                   "t_alpha","True_A","Frac_A"),]

  return(result)
}
