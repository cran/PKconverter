#' Convert pharmacokinetic parameters for two compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distribution(V1),
#' clearance (Cl1) and half-lives (t_alpha and t_beta)
#'
#' @usage TwoComp_Volume_Clearance_HalfLife(V1,Cl1,t_alpha,t_beta,
#'  V1.sd=NA,Cl1.sd=NA,t_alpha.sd=NA,
#'  t_beta.sd=NA,covar=c(V1Cl1=NA,V1talpha=NA,V1tbeta=NA,Cl1talpha=NA,
#'    Cl1tbeta=NA,talphatbeta=NA),...)
#' @param V1 The volume of distribution of compartment 1
#' @param Cl1 Clearance from compartment 1
#' @param t_alpha half life of compartment 1
#' @param t_beta half life of compartment 2
#' @param V1.sd standard error of V1
#' @param Cl1.sd standard error of Cl1
#' @param t_alpha.sd standard error of t_alpha
#' @param t_beta.sd standard error of t_beta
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' TwoComp_Volume_Clearance_HalfLife(V1=5,Cl1=3.5,t_alpha=0.568,t_beta=24.2,
#' V1.sd=0.01,Cl1.sd=0.01,t_alpha.sd=0.002,t_beta.sd=0.5)
TwoComp_Volume_Clearance_HalfLife<-function(V1,Cl1,t_alpha,t_beta,
  V1.sd=NA,Cl1.sd=NA,t_alpha.sd=NA,
  t_beta.sd=NA,covar=c(V1Cl1=NA,V1talpha=NA,V1tbeta=NA,Cl1talpha=NA,
    Cl1tbeta=NA,talphatbeta=NA),...){
  if(is.na(covar[1])) covar<-rep(0,6)
  V1.var = (V1.sd)^2;                  Cl1.var = (Cl1.sd)^2
  t_alpha.var = (t_alpha.sd)^2;        t_beta.var = (t_beta.sd)^2

  f.V2<-quote(quote((V1)*((log(2)/t_alpha)+(log(2)/t_beta)-((log(2)/t_alpha)*
              (log(2)/t_beta)/(Cl1/V1))-(Cl1/V1))/((log(2)/t_alpha)*
              (log(2)/t_beta)/(Cl1/V1))))
  V2<-eval(eval(f.V2))
  ff.V2<-stats::as.formula(paste("~",as.character(f.V2[2],"")))

  f.Vd<-quote(quote((V1)+(V1*((log(2)/t_alpha)+
                 (log(2)/t_beta)-((log(2)/t_alpha)*(log(2)/t_beta)/(Cl1/V1))-
                 (Cl1/V1))/((log(2)/t_alpha)*(log(2)/t_beta)/(Cl1/V1)))))

  Vdss<-eval(eval(f.Vd))
  ff.Vd<-stats::as.formula(paste("~",as.character(f.Vd[2],"")))

  V2_deriv<-as.matrix(attr(eval(stats::deriv(ff.V2,
                            c("V1","Cl1","t_alpha","t_beta"))),"gradient"))
  Vdss_deriv<-as.matrix(attr(eval(stats::deriv(ff.Vd,
                            c("V1","Cl1","t_alpha","t_beta"))),"gradient"))

  sigma4<-matrix(as.numeric(c(V1.var,covar[1],covar[2],covar[3],covar[1],
                              Cl1.var,covar[4],covar[5], covar[2],covar[4],
                              t_alpha.var,covar[6],covar[3],covar[5],covar[6],
                              t_beta.var)),4,4,byrow=T)

  V2.sd<-sqrt(V2_deriv %*% sigma4 %*% t(V2_deriv))
  Vdss.sd<-sqrt(Vdss_deriv %*% sigma4 %*% t(Vdss_deriv))

  f.Cl2<-quote(quote(V1*((log(2)/t_alpha)+(log(2)/t_beta)-
               ((log(2)/t_alpha)*(log(2)/t_beta)/(Cl1/V1))-(Cl1/V1))))
  Cl2<-eval(eval(f.Cl2))
  ff.Cl2<-stats::as.formula(paste("~",as.character(f.Cl2[2],"")))
  Cl2_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl2,
                          c("V1","Cl1","t_alpha","t_beta"))),"gradient"))

  Cl2.sd<-sqrt(Cl2_deriv %*% sigma4 %*% t(Cl2_deriv))


  k10<-Cl1/V1
  f.k12<-quote(quote((log(2)/t_alpha)+(log(2)/t_beta)-
                        ((log(2)/t_alpha)*(log(2)/t_beta)/(Cl1/V1))-(Cl1/V1)))
  k12<-eval(eval(f.k12))
  ff.k12<-stats::as.formula(paste("~",as.character(f.k12[2],"")))
  f.k21<-quote(quote((log(2)/(t_alpha))*(log(2)/(t_beta))/(Cl1/V1)))
  k21<-eval(eval(f.k21))
  ff.k21<-stats::as.formula(paste("~",as.character(f.k21[2],"")))

  sigma2<-matrix(as.numeric(c(V1.var,covar[1],covar[1],Cl1.var)),2,2,byrow=T)

  k10_deriv<-as.matrix(attr(eval(stats::deriv(~Cl1/V1,c("V1","Cl1"))),
                            "gradient"))
  k12_deriv<-as.matrix(attr(eval(stats::deriv(ff.k12,
                          c("V1","Cl1","t_alpha","t_beta"))),"gradient"))
  k21_deriv<-as.matrix(attr(eval(stats::deriv(ff.k21,
                          c("V1","Cl1","t_alpha","t_beta"))),"gradient"))

  k10.sd<-sqrt(k10_deriv %*% sigma2 %*% t(k10_deriv))
  k12.sd<-sqrt(k12_deriv %*% sigma4 %*% t(k12_deriv))
  k21.sd<-sqrt(k21_deriv %*% sigma4 %*% t(k21_deriv))

  f.true_A<-quote(quote((((log(2)/t_alpha)*(log(2)/t_beta)/(Cl1/V1))-
                 (log(2)/t_alpha))/((log(2)/t_beta)-(log(2)/t_alpha))/V1))
  true_A<-eval(eval(f.true_A))
  ff.true_A<-stats::as.formula(paste("~",as.character(f.true_A[2],"")))
  f.true_B<-quote(quote((((log(2)/t_alpha)*(log(2)/t_beta)/(Cl1/V1))-
                 (log(2)/t_beta))/(-((log(2)/t_beta)-(log(2)/t_alpha)))/V1))
  true_B<-eval(eval(f.true_B))
  ff.true_B<-stats::as.formula(paste("~",as.character(f.true_B[2],"")))

  true_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_A,
                               c("V1","Cl1","t_alpha","t_beta"))),"gradient"))
  true_A.sd<-sqrt(true_A_deriv %*% sigma4 %*% t(true_A_deriv))
  true_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_B,
                               c("V1","Cl1","t_alpha","t_beta"))),"gradient"))
  true_B.sd<-sqrt(true_B_deriv %*% sigma4 %*% t(true_B_deriv))

  f.frac_A<-quote(quote((((log(2)/t_alpha)*(log(2)/t_beta)/(Cl1/V1))-
             (log(2)/t_alpha))/((log(2)/t_beta)-(log(2)/t_alpha))))
  frac_A<-eval(eval(f.frac_A))
  ff.frac_A<-stats::as.formula(paste("~",as.character(f.frac_A[2],"")))
  f.frac_B<-quote(quote((((log(2)/t_alpha)*(log(2)/t_beta)/(Cl1/V1))-
             (log(2)/t_beta))/(-((log(2)/t_beta)-(log(2)/t_alpha)))))
  frac_B<-eval(eval(f.frac_B))
  ff.frac_B<-stats::as.formula(paste("~",as.character(f.frac_B[2],"")))

  frac_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_A,
                                c("V1","Cl1","t_alpha","t_beta"))),"gradient"))
  frac_A.sd<-sqrt(frac_A_deriv %*% sigma4 %*% t(frac_A_deriv))
  frac_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_B,
                                c("V1","Cl1","t_alpha","t_beta"))),"gradient"))
  frac_B.sd<-sqrt(frac_B_deriv %*% sigma4 %*% t(frac_B_deriv))

  alpha<-log(2)/(t_alpha)
  beta<-log(2)/(t_beta)

  alpha_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/t_alpha,"t_alpha")),
                              "gradient"))
  beta_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/t_beta,"t_beta")),
                             "gradient"))

  alpha.sd<-sqrt(alpha_deriv * t_alpha.var * alpha_deriv)
  beta.sd<-sqrt(beta_deriv * t_beta.var *beta_deriv)
  if(is.na(Cl1[1])){
    param = rep(NA,16)
    sd = rep(NA,16)
  } else{
    param = c(V1,Cl1,t_alpha,t_beta,V2,Cl2,k10,k12,k21,Vdss,true_A,true_B,
              frac_A,frac_B,alpha,beta)
    sd = c(V1.sd,Cl1.sd,t_alpha.sd,t_beta.sd,V2.sd,Cl2.sd,k10.sd,k12.sd,
           k21.sd,Vdss.sd,true_A.sd,true_B.sd,frac_A.sd,frac_B.sd,alpha.sd,
           beta.sd)
  }
  result = data.frame(Parameter=c("V1","Cl1","t_alpha","t_beta","V2","Cl2",
                      "k10","k12","k21","Vdss","True_A","True_B","Frac_A",
                      "Frac_B","alpha","beta"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("V1","Cl1","t_alpha","t_beta","V2","Cl2","k10","k12",
     "k21","Vdss","True_A","True_B","Frac_A","Frac_B","alpha","beta")
  result<-result[c("Vdss","V1","V2","Cl1","Cl2","k10","k12","k21",
                   "alpha","beta","t_alpha","t_beta",
                   "True_A","True_B","Frac_A","Frac_B"),]
  return(result)
}
