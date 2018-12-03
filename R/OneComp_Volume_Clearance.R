#' Convert pharmacokinetic parameters for one compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distribution (V1) and
#' clearance (Cl1)
#' @usage OneComp_Volume_Clearance(V1,Cl1,V1.sd=NA,Cl1.sd=NA,covar=c(V1Cl1=NA))
#' @param V1 The volume of distribution of compartment 1
#' @param Cl1 Clearance from compartment 1
#' @param V1.sd standard error of V1
#' @param Cl1.sd standard error of Cl1
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{www.nonmemcourse.com/convert.xls}
#' @references Gabrielsson and Weiner(2006) Pharmacokinetic and
#'             Pharmacodynamic Data Analysis: Concepts and Applications,
#'             Swedish Academy of Pharmaceutical Sciences.
#' @export
#' @examples
#' OneComp_Volume_Clearance(V1=8,Cl1=4,V1.sd=0.01,Cl1.sd=0.01)

OneComp_Volume_Clearance<-function(V1,Cl1,V1.sd=NA,Cl1.sd=NA,covar=c(V1Cl1=NA)){
  if(is.na(covar[1])) covar<-0
  Vdss = V1
  k10 = Cl1/V1
  t_alpha = (log(2)*V1)/Cl1
  true_A = 1/V1
  frac_A = ifelse(is.na(V1[1]),NA,1)
  frac_A.sd = ifelse(is.na(V1.sd[1]),NA,0)

  alpha = k10
  V1.var = (V1.sd)^2
  Cl1.var = (Cl1.sd)^2

  k10_deriv<-as.matrix(attr(eval(stats::deriv(~Cl1/V1,c("V1","Cl1"))),
                            "gradient"))
  t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(~(log(2)*V1)/Cl1,
                                               c("V1","Cl1"))),"gradient"))
  true_A_deriv<-as.matrix(attr(eval(stats::deriv(~1/V1,"V1")),"gradient"))
  sigma<-matrix(as.numeric(c(V1.var,covar[1],covar[1],Cl1.var)),2,2,byrow=T)

  Vdss.sd = V1.sd
  k10.sd = sqrt(k10_deriv %*% sigma %*% t(k10_deriv))
  t_alpha.sd = sqrt(t_alpha_deriv %*% sigma %*% t(t_alpha_deriv))
  true_A.sd = sqrt(true_A_deriv * as.matrix(V1.sd) *true_A_deriv)

  alpha.sd = k10.sd

  param = c(V1,Cl1,k10,t_alpha,Vdss,true_A,frac_A,alpha)
  sd = c(V1.sd,Cl1.sd,k10.sd,t_alpha.sd,Vdss.sd,true_A.sd,frac_A.sd,alpha.sd)

  result = data.frame(Parameter=c("V1","Cl1","k10","t_alpha","Vdss",
    "True_A","Frac_A","alpha"),Estimate=param, Std.err=sd)
  row.names(result) <- c("V1","Cl1","k10","t_alpha","Vdss",
    "True_A","Frac_A","alpha")
  result<-result[c("Vdss","V1","Cl1","k10","alpha",
                    "t_alpha","True_A","Frac_A"),]
  return(result)
}
