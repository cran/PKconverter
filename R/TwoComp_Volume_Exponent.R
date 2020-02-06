#' Convert pharmacokinetic parameters for two compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distribution(V1),
#' transfer rate constant (k12), and parameters (alpha and beta)
#' in the model "Aexp(-alpha)+Bexp(-beta)"
#' @usage TwoComp_Volume_Exponent(V1,alpha,beta,k21,V1.sd=NA,
#'  alpha.sd=NA,beta.sd=NA,k21.sd=NA,
#'  covar=c(V1alpha=NA,V1beta=NA,V1k21=NA,alphabeta=NA,
#'  alphak21=NA,betak21=NA),...)
#' @param V1 The volume of distribution of compartment 1
#' @param alpha parameter in one compartment model "Aexp(-alpha)"
#' @param beta parameter in two compartment model "Aexp(-alpha)+Bexp(-beta)"
#' @param k21 transfer rate constants from compartment 2 to compartment 1
#' @param V1.sd standard error of V1
#' @param alpha.sd standard error of alpha
#' @param beta.sd standard error of beta
#' @param k21.sd standard error of k21
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' TwoComp_Volume_Exponent(V1=5,alpha=1.221, beta=0.029, k21=0.05,
#' V1.sd=0.01,alpha.sd=0.01,beta.sd=0.00005,k21.sd=0.0006)


TwoComp_Volume_Exponent<-function(V1,alpha,beta,k21,V1.sd=NA,alpha.sd=NA,
  beta.sd=NA,k21.sd=NA,
  covar=c(V1alpha=NA,V1beta=NA,V1k21=NA,
          alphabeta=NA,alphak21=NA,betak21=NA),...){
  if(is.na(covar[1])) covar<-rep(0,6)
  V1.var = (V1.sd)^2;            alpha.var = (alpha.sd)^2
  beta.var = (beta.sd)^2;        k21.var = (k21.sd)^2

  f.V2<-quote(quote(V1*(alpha+beta-k21-(alpha*beta/k21))/k21))
  V2<-eval(eval(f.V2))
  ff.V2<-stats::as.formula(paste("~",as.character(f.V2[2],"")))
  f.Vdss<-quote(quote(V1+(V1*(alpha+beta-k21-(alpha*beta/k21))/k21)))
  Vdss<-eval(eval(f.Vdss))
  ff.Vdss<-stats::as.formula(paste("~",as.character(f.Vdss[2],"")))

  V2_deriv<-as.matrix(attr(eval(stats::deriv(ff.V2,
                                    c("V1","alpha","beta","k21"))),"gradient"))
  Vdss_deriv<-as.matrix(attr(eval(stats::deriv(ff.Vdss,
                                    c("V1","alpha","beta","k21"))),"gradient"))
  sigma4<-matrix(as.numeric(c(V1.var,covar[1],covar[2],covar[3],
                           covar[1],alpha.var,covar[4],covar[5],
                           covar[2],covar[4],beta.var,covar[6],
                           covar[3],covar[5],covar[6],k21.var)),4,4,byrow=T)
  V2.sd<-sqrt(V2_deriv %*% sigma4 %*% t(V2_deriv))
  Vdss.sd<-sqrt(Vdss_deriv %*% sigma4 %*% t(Vdss_deriv))

  Cl1<-V1*(alpha*beta/k21)
  Cl1_deriv<-as.matrix(attr(eval(stats::deriv(~V1*(alpha*beta/k21),
                                    c("V1","alpha","beta","k21"))),"gradient"))
  Cl1.sd<-sqrt(Cl1_deriv %*% sigma4 %*% t(Cl1_deriv))


  f.Cl2<-quote(quote(V1*(alpha+beta-k21-(alpha*beta/k21))))
  Cl2<-eval(eval(f.Cl2))
  ff.Cl2<-stats::as.formula(paste("~",as.character(f.Cl2[2],"")))
  Cl2_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl2,
                                    c("V1","alpha","beta","k21"))),"gradient"))
  Cl2.sd<-sqrt(Cl2_deriv %*% sigma4 %*% t(Cl2_deriv))

  k10<-alpha*beta/k21
  k12<-alpha+beta-k21-(alpha*beta/k21)
  sigma3<-matrix(as.numeric(c(alpha.var,covar[4],covar[5],
                               covar[4],beta.var,covar[6],
                               covar[5],covar[6],k21.var)),3,3,byrow=T)

  k10_deriv<-as.matrix(attr(eval(stats::deriv(~alpha*beta/k21,
        c("alpha","beta","k21"))),"gradient"))
  k12_deriv<-as.matrix(attr(eval(stats::deriv(~alpha+beta-k21-
        (alpha*beta/k21),c("alpha","beta","k21"))),"gradient"))

  k10.sd<-sqrt(k10_deriv %*% sigma3 %*% t(k10_deriv))
  k12.sd<-sqrt(k12_deriv %*% sigma3 %*% t(k12_deriv))

  t_alpha<-log(2)/alpha
  t_beta<-log(2)/beta
  t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/alpha,"alpha")),
                                "gradient"))
  t_beta_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/beta,"beta")),
                               "gradient"))
  t_alpha.sd<-sqrt(t_alpha_deriv * alpha.var * t_alpha_deriv)
  t_beta.sd<-sqrt(t_beta_deriv * beta.var * t_beta_deriv)

  True_A<-(k21-alpha)/(beta-alpha)/V1;
  True_B<-(k21-beta)/(alpha-beta)/V1
  True_A_deriv<-as.matrix(attr(eval(stats::deriv(~(k21-alpha)/(beta-alpha)/V1,
        c("V1","alpha","beta","k21"))),"gradient"))
  True_B_deriv<-as.matrix(attr(eval(stats::deriv(~(k21-beta)/(alpha-beta)/V1,
        c("V1","alpha","beta","k21"))),"gradient"))
  True_A.sd<-sqrt(True_A_deriv %*% sigma4 %*% t(True_A_deriv))
  True_B.sd<-sqrt(True_B_deriv %*% sigma4 %*% t(True_B_deriv))

  Frac_A<-(k21-alpha)/(beta-alpha)
  Frac_B<-(k21-beta)/(alpha-beta)
  Frac_A_deriv<-as.matrix(attr(eval(stats::deriv(~(k21-alpha)/(beta-alpha),
        c("V1","alpha","beta","k21"))),"gradient"))
  Frac_B_deriv<-as.matrix(attr(eval(stats::deriv(~(k21-beta)/(alpha-beta),
        c("V1","alpha","beta","k21"))),"gradient"))
  Frac_A.sd<-sqrt(Frac_A_deriv %*% sigma4 %*% t(Frac_A_deriv))
  Frac_B.sd<-sqrt(Frac_B_deriv %*% sigma4 %*% t(Frac_B_deriv))
  if(is.na(V1[1])){
    param = rep(NA,16)
    sd = rep(NA,16)
  } else{
    param = c(V1,alpha,beta,k21,V2,Cl1,Cl2,k10,k12,t_alpha,t_beta,Vdss,
              True_A,True_B, Frac_A,Frac_B)
    sd = c(V1.sd,alpha.sd,beta.sd,k21.sd,V2.sd,Cl1.sd,Cl2.sd,k10.sd,k12.sd,
           t_alpha.sd,t_beta.sd,Vdss.sd,True_A.sd,True_B.sd,Frac_A.sd,Frac_B.sd)
  }
  result = data.frame(Parameter=c("V1","alpha","beta","k21","V2","Cl1","Cl2",
                      "k10","k12","t_alpha","t_beta","Vdss","True_A","True_B",
                      "Frac_A","Frac_B"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("V1","alpha","beta","k21","V2","Cl1","Cl2","k10",
                         "k12","t_alpha","t_beta","Vdss","True_A","True_B",
                         "Frac_A","Frac_B")
  result<-result[c("Vdss","V1","V2","Cl1","Cl2","k10","k12","k21",
                   "alpha","beta","t_alpha","t_beta",
                   "True_A","True_B","Frac_A","Frac_B"),]
  return(result)
}
