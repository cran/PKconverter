#' Convert pharmacokinetic parameters for two compartment model
#'
#' Calculate pharmacokinetic parameters with parameters (A, B, alpha and beta)
#' in two compartment model "Aexp(-alpha)+Bexp(-beta)"
#'
#' @usage TwoComp_Coefficient_Exponent(A,B,alpha,beta,
#'  A.sd=NA,B.sd=NA,alpha.sd=NA,beta.sd=NA,
#'  covar=c(AB=NA,Aalpha=NA,Abeta=NA,Balpha=NA,Bbeta=NA,alphabeta=NA),...)
#' @param A parameter in one compartment model "Aexp(-alpha)"
#' @param B parameter in two compartment model "Aexp(-alpha)+Bexp(-beta)"
#' @param alpha parameter in one compartment model "Aexp(-alpha)"
#' @param beta parameter in two compartment model "Aexp(-alpha)+Bexp(-beta)"
#' @param A.sd standard error of A
#' @param B.sd standard error of B
#' @param alpha.sd standard error of alpha
#' @param beta.sd standard error of beta
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' TwoComp_Coefficient_Exponent(A=0.196,B=0.0036,alpha=1.221,beta=0.0287,
#' A.sd=0.002,B.sd=0.00005,alpha.sd=0.09,beta.sd=0.0006)

TwoComp_Coefficient_Exponent<-function(A,B,alpha,beta,
  A.sd=NA,B.sd=NA,alpha.sd=NA,beta.sd=NA,
  covar=c(AB=NA,Aalpha=NA,Abeta=NA,Balpha=NA,Bbeta=NA,alphabeta=NA),...){
  if(is.na(covar[1])) covar<-rep(0,6)
  A.var = (A.sd)^2;             B.var = (B.sd)^2
  alpha.var = (alpha.sd)^2;     beta.var = (beta.sd)^2

  V1<-1/(A+B)
  f.V2<-quote(quote((1/(A+B))*(alpha+beta-((A*beta+B*alpha)/(A+B))-
             (alpha*beta/((A*beta+B*alpha)/(A+B))))/((A*beta+B*alpha)/(A+B))))

  V2<-eval(eval(f.V2))
  ff.V2<-stats::as.formula(paste("~",as.character(f.V2[2],"")))

  f.Vdss<-quote(quote((1/(A+B))+((1/(A+B))*(alpha+beta-((A*beta+B*alpha)/(A+B))-
             (alpha*beta/((A*beta+B*alpha)/(A+B))))/((A*beta+B*alpha)/(A+B)))))

  Vdss<-eval(eval(f.Vdss))
  ff.Vdss<-stats::as.formula(paste("~",as.character(f.Vdss[2],"")))

  V1_deriv<-as.matrix(attr(eval(stats::deriv(~1/(A+B),c("A","B"))),"gradient"))
  V2_deriv<-as.matrix(attr(eval(stats::deriv(ff.V2,
                                c("A","B","alpha","beta"))),"gradient"))
  Vdss_deriv<-as.matrix(attr(eval(stats::deriv(ff.Vdss,
                                  c("A","B","alpha","beta"))),"gradient"))


  sigma2<-matrix(as.numeric(c(A.var,covar[1],covar[1],B.var)),2,2,byrow=T)

  sigma4<-matrix(as.numeric(c(A.var,covar[1],covar[2],covar[3],
                              covar[1],B.var,covar[4],covar[5],
                              covar[2],covar[4],alpha.var,covar[6],
                              covar[3],covar[5],covar[6],beta.var)),4,4,byrow=T)

  V1.sd<-sqrt(V1_deriv %*% sigma2 %*% t(V1_deriv))
  V2.sd<-sqrt(V2_deriv %*% sigma4 %*% t(V2_deriv))
  Vdss.sd<-sqrt(Vdss_deriv %*% sigma4 %*% t(Vdss_deriv))

  f.Cl1<-quote(quote((1/(A+B))*(alpha*beta/((A*beta+B*alpha)/(A+B)))))

  Cl1<-eval(eval(f.Cl1))
  ff.Cl1<-stats::as.formula(paste("~",as.character(f.Cl1[2],"")))
  f.Cl2<-quote(quote((1/(A+B))*(alpha+beta-((A*beta+B*alpha)/(A+B))-
                                   (alpha*beta/((A*beta+B*alpha)/(A+B))))))
  Cl2<-eval(eval(f.Cl2))
  ff.Cl2<-stats::as.formula(paste("~",as.character(f.Cl2[2],"")))
  Cl1_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl1,
    c("A","B","alpha","beta"))),"gradient"))
  Cl2_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl2,
    c("A","B","alpha","beta"))),"gradient"))

  Cl1.sd<-sqrt(Cl1_deriv %*% sigma4 %*% t(Cl1_deriv))
  Cl2.sd<-sqrt(Cl2_deriv %*% sigma4 %*% t(Cl2_deriv))

  f.k10<-quote(quote(alpha*beta/((A*beta+B*alpha)/(A+B))))

  k10<-eval(eval(f.k10))
  ff.k10<-stats::as.formula(paste("~",as.character(f.k10[2],"")))
  f.k12<-quote(quote(alpha+beta-((A*beta+B*alpha)/(A+B))-
                        (alpha*beta/((A*beta+B*alpha)/(A+B)))))
  k12<-eval(eval(f.k12))
  ff.k12<-stats::as.formula(paste("~",as.character(f.k12[2],"")))
  f.k21<-quote(quote((A*beta+B*alpha)/(A+B)))
  k21<-eval(eval(f.k21))
  ff.k21<-stats::as.formula(paste("~",as.character(f.k21[2],"")))
  k10_deriv<-as.matrix(attr(eval(stats::deriv(ff.k10,
                                       c("A","B","alpha","beta"))),"gradient"))
  k12_deriv<-as.matrix(attr(eval(stats::deriv(ff.k12,
                                       c("A","B","alpha","beta"))),"gradient"))
  k21_deriv<-as.matrix(attr(eval(stats::deriv(ff.k21,
                                       c("A","B","alpha","beta"))),"gradient"))

  k10.sd<-sqrt(k10_deriv %*% sigma4 %*% t(k10_deriv))
  k12.sd<-sqrt(k12_deriv %*% sigma4 %*% t(k12_deriv))
  k21.sd<-sqrt(k21_deriv %*% sigma4 %*% t(k21_deriv))

  t_alpha<-log(2)/alpha
  t_beta<-log(2)/beta
  t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/alpha,"alpha")),
                                "gradient"))
  t_beta_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/beta,"beta")),
                               "gradient"))
  t_alpha.sd<-sqrt(t_alpha_deriv * alpha.var * t_alpha_deriv)
  t_beta.sd<-sqrt(t_beta_deriv * beta.var * t_beta_deriv)

  Frac_A<-A*V1
  Frac_B<-B*V1
  Frac_A_deriv<-as.matrix(attr(eval(stats::deriv(~A/(A+B),c("A","B"))),
                               "gradient"))
  Frac_B_deriv<-as.matrix(attr(eval(stats::deriv(~B/(A+B),c("A","B"))),
                               "gradient"))

  Frac_A.sd<-sqrt(Frac_A_deriv %*% sigma2 %*% t(Frac_A_deriv))
  Frac_B.sd<-sqrt(Frac_B_deriv%*% sigma2 %*% t(Frac_B_deriv))
  if(is.na(A[1])){
    param = rep(NA,16)
    sd = rep(NA,16)
  } else{
    param = c(A,B,alpha,beta,V1,V2,Cl1,Cl2,k10,k12,k21,t_alpha,t_beta,Vdss,
              Frac_A,Frac_B)
    sd = c(A.sd,B.sd,alpha.sd,beta.sd,V1.sd,V2.sd,Cl1.sd,Cl2.sd,k10.sd,k12.sd,
           k21.sd,t_alpha.sd,t_beta.sd,Vdss.sd,Frac_A.sd,Frac_B.sd)
  }
  result = data.frame(Parameter=c("True_A","True_B","alpha","beta","V1","V2",
                      "Cl1","Cl2","k10","k12","k21","t_alpha","t_beta","Vdss",
                      "Frac_A","Frac_B"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("True_A","True_B","alpha","beta","V1","V2","Cl1",
                         "Cl2","k10","k12","k21","t_alpha","t_beta","Vdss",
                         "Frac_A","Frac_B")
  result<-result[c("Vdss","V1","V2", "Cl1","Cl2","k10","k12","k21",
                   "alpha","beta","t_alpha","t_beta","True_A","True_B",
                   "Frac_A","Frac_B"),]
  return(result)
}
