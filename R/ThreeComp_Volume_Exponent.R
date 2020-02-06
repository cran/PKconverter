
#' Convert pharmacokinetic parameters for three compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distribution(V1),
#' transfer rate constant (k12 and k31), and parameters
#' (alpha, beta and gamma) in the model "Aexp(-alpha)+Bexp(-beta)+Cexp(-gamma)"
#' @usage ThreeComp_Volume_Exponent(V1,alpha,beta,gamma,k21,k31,
#'  V1.sd=NA,alpha.sd=NA,beta.sd=NA,gamma.sd=NA,k21.sd=NA,k31.sd=NA,
#'  covar=c(V1alpha=NA,V1beta=NA,V1gamma=NA,V1k21=NA,V1k31=NA,
#'    alphabeta=NA,alphagamma=NA,alphak21=NA,alphak31=NA,
#'    betagamma=NA,betak21=NA,betak31=NA,gammak21=NA,gammak31=NA,
#'    k21k31=NA),...)
#' @param V1 The volume of distribution of compartment 1
#' @param k21 transfer rate constants from compartment 2 to compartment 1
#' @param k31 transfer rate constants from compartment 3 to compartment 1
#' @param alpha parameter in one compartment model "Aexp(-alpha)"
#' @param beta parameter in two compartment model "Aexp(-alpha)+Bexp(-beta)"
#' @param gamma parameter in three compartment model
#'              "Aexp(-alpha)+Bexp(-beta)+Cexp(-gamma)"
#' @param V1.sd standard error of V1
#' @param k21.sd standard error of k21
#' @param k31.sd standard error of k31
#' @param alpha.sd standard error of alpha
#' @param beta.sd standard error of beta
#' @param gamma.sd standard error of gamma
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' ThreeComp_Volume_Exponent(V1=10,alpha=0.6, beta=0.013, gamma=0.00074,
#'    k21=0.02, k31=0.001, V1.sd=0.01,alpha.sd=0.01,beta.sd=0.00005,
#'    gamma.sd=0.000002, k21.sd=0.0006,k31.sd=0.0000005)

ThreeComp_Volume_Exponent<-function(V1,alpha,beta,gamma,k21,k31,
  V1.sd=NA,alpha.sd=NA,beta.sd=NA,gamma.sd=NA,k21.sd=NA,k31.sd=NA,
  covar=c(V1alpha=NA,V1beta=NA,V1gamma=NA,V1k21=NA,V1k31=NA,
    alphabeta=NA,alphagamma=NA,alphak21=NA,alphak31=NA,
    betagamma=NA,betak21=NA,betak31=NA,gammak21=NA,gammak31=NA,
    k21k31=NA),...){
  if(is.na(covar[1])) covar<-rep(0,15)
  V1.var = (V1.sd)^2;       alpha.var = (alpha.sd)^2;  beta.var = (beta.sd)^2;
  gamma.var = (gamma.sd)^2; k21.var = (k21.sd)^2;      k31.var = (k31.sd)^2;

  f.V2<-quote(quote(V1*(((beta*gamma+alpha*beta+alpha*gamma)-k21*
            (alpha+beta+gamma)-(alpha*beta*gamma/k21/k31)*k31+k21*k21)/
            (k31-k21))/k21))

  V2<-eval(eval(f.V2))
  ff.V2<-stats::as.formula(paste("~",as.character(f.V2[2],"")))
  f.V3<-quote(quote(V1*(alpha+beta+gamma-((alpha*beta*gamma/k21/k31)+
        (((beta*gamma+alpha*beta+alpha*gamma)-k21*(alpha+beta+gamma)-
        (alpha*beta*gamma/k21/k31)*k31+k21*k21)/(k31-k21))+k21+k31))/k31))

  V3<-eval(eval(f.V3))
  ff.V3<-stats::as.formula(paste("~",as.character(f.V3[2],"")))

  V2_deriv<-as.matrix(attr(eval(stats::deriv(ff.V2,
                      c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))
  V3_deriv<-as.matrix(attr(eval(stats::deriv(ff.V3,
                      c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))

  f.Vdss<-quote(quote((V1)+(V1*(((beta*gamma+alpha*beta+alpha*gamma)-k21*
            (alpha+beta+gamma)-(alpha*beta*gamma/k21/k31)*k31+k21*k21)/
            (k31-k21))/k21)+((V1*(alpha+beta+gamma-((alpha*beta*gamma/k21/k31)+
        (((beta*gamma+alpha*beta+alpha*gamma)-k21*(alpha+beta+gamma)-
        (alpha*beta*gamma/k21/k31)*k31+k21*k21)/(k31-k21))+k21+k31))/k31))))

   Vdss<-V1+V2+V3
   ff.Vdss<-stats::as.formula(paste("~",as.character(f.Vdss[2],"")))

   Vdss_deriv<-as.matrix(attr(eval(stats::deriv(ff.Vdss,
                     c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))

   sigma6<-matrix(as.numeric(c(V1.var,covar[1],covar[2],covar[3],covar[4],
                  covar[5],covar[1],alpha.var,covar[6],covar[7],covar[8],
                  covar[9],covar[2],covar[6],beta.var,covar[9],covar[11],
                  covar[12], covar[3],covar[7],covar[10],gamma.var,covar[13],
                  covar[14],covar[4],covar[8],covar[11],covar[13],k21.var,
                  covar[15],covar[5],covar[9],covar[12],covar[14],covar[15],
                  k31.var)),6,6,byrow=T)

   V2.sd<-sqrt(V2_deriv %*% sigma6 %*% t(V2_deriv))
   V3.sd<-sqrt(V3_deriv %*% sigma6 %*% t(V3_deriv))
   Vdss.sd<-sqrt(Vdss_deriv %*% sigma6 %*% t(Vdss_deriv))


   f.Cl1<-quote(quote(V1*(alpha*beta*gamma/k21/k31)))
   Cl1<-eval(eval(f.Cl1))
   ff.Cl1<-stats::as.formula(paste("~",as.character(f.Cl1[2],"")))
   f.Cl2<-quote(quote(V1*(((beta*gamma+alpha*beta+alpha*gamma)-
        k21*(alpha+beta+gamma)-
        (alpha*beta*gamma/k21/k31)*k31+k21*k21)/(k31-k21))))
   Cl2<-eval(eval(f.Cl2))
   ff.Cl2<-stats::as.formula(paste("~",as.character(f.Cl2[2],"")))
   f.Cl3<-quote(quote(V1*(alpha+beta+gamma-((alpha*beta*gamma/k21/k31)+
        (((beta*gamma+alpha*beta+alpha*gamma)-k21*(alpha+beta+gamma)-
        (alpha*beta*gamma/k21/k31)*k31+k21*k21)/(k31-k21))+k21+k31))))
   Cl3<-eval(eval(f.Cl3))
   ff.Cl3<-stats::as.formula(paste("~",as.character(f.Cl3[2],"")))

   Cl1_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl1,
        c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))
   Cl2_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl2,
        c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))
   Cl3_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl3,
        c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))

   Cl1.sd<-sqrt(Cl1_deriv %*% sigma6 %*% t(Cl1_deriv))
   Cl2.sd<-sqrt(Cl2_deriv %*% sigma6 %*% t(Cl2_deriv))
   Cl3.sd<-sqrt(Cl3_deriv %*% sigma6 %*% t(Cl3_deriv))

   t_alpha<-log(2)/alpha
   t_beta<-log(2)/beta
   t_gamma<-log(2)/gamma

   t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/alpha,"alpha")),
                                 "gradient"))
   t_beta_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/beta,"beta")),
                                "gradient"))
   t_gamma_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/gamma,"gamma")),
                                 "gradient"))

   t_alpha.sd<-sqrt(t_alpha_deriv * alpha.var * t_alpha_deriv)
   t_beta.sd<-sqrt(t_beta_deriv * beta.var * t_beta_deriv)
   t_gamma.sd<-sqrt(t_gamma_deriv * gamma.var * t_gamma_deriv)


   f.k10<-quote(quote(alpha*beta*gamma/k21/k31))
   k10<-eval(eval(f.k10))
   ff.k10<-stats::as.formula(paste("~",as.character(f.k10[2],"")))
   f.k12<-quote(quote(((beta*gamma+alpha*beta+alpha*gamma)-
                        k21*(alpha+beta+gamma)-
                        (alpha*beta*gamma/k21/k31)*k31+k21*k21)/(k31-k21)))
   k12<-eval(eval(f.k12))
   ff.k12<-stats::as.formula(paste("~",as.character(f.k12[2],"")))
   f.k13<-quote(quote(alpha+beta+gamma-((alpha*beta*gamma/k21/k31)+
        (((beta*gamma+alpha*beta+alpha*gamma)-k21*(alpha+beta+gamma)-
        (alpha*beta*gamma/k21/k31)*k31+k21*k21)/(k31-k21))+k21+k31)))
   k13<-eval(eval(f.k13))
   ff.k13<-stats::as.formula(paste("~",as.character(f.k13[2],"")))

   sigma5<-matrix(as.numeric(c(alpha.var,covar[6],covar[7],covar[8],covar[9],
                              covar[6],beta.var,covar[10],covar[11],covar[12],
                              covar[7],covar[10],gamma.var,covar[13],covar[14],
                              covar[8],covar[11],covar[13],k21.var,covar[14],
                              covar[9],covar[12],covar[14],covar[15],k31.var)),
                              5,5,byrow=T)

   k10_deriv<-as.matrix(attr(eval(stats::deriv(ff.k10,
                          c("alpha","beta","gamma","k21","k31"))),"gradient"))
   k12_deriv<-as.matrix(attr(eval(stats::deriv(ff.k12,
                          c("alpha","beta","gamma","k21","k31"))),"gradient"))
   k13_deriv<-as.matrix(attr(eval(stats::deriv(ff.k13,
                          c("alpha","beta","gamma","k21","k31"))),"gradient"))

   k10.sd<-sqrt(k10_deriv %*% sigma5 %*% t(k10_deriv))
   k12.sd<-sqrt(k12_deriv %*% sigma5 %*% t(k12_deriv))
   k13.sd<-sqrt(k13_deriv %*% sigma5 %*% t(k13_deriv))


   True_A<-(k21-alpha)*(k31-alpha)/(alpha-beta)/(alpha-gamma)/V1
   True_B<-(k21-beta)*(k31-beta)/(beta-alpha)/(beta-gamma)/V1
   True_C<-(k21-gamma)*(k31-gamma)/(gamma-beta)/(gamma-alpha)/V1

   True_A_deriv<-as.matrix(attr(eval(stats::deriv(~(k21-alpha)*(k31-alpha)/
        (alpha-beta)/(alpha-gamma)/V1,
        c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))
   True_B_deriv<-as.matrix(attr(eval(stats::deriv(~(k21-beta)*(k31-beta)/
        (beta-alpha)/(beta-gamma)/V1,
        c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))
   True_C_deriv<-as.matrix(attr(eval(stats::deriv(~(k21-gamma)*(k31-gamma)/
        (gamma-beta)/(gamma-alpha)/V1,
        c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))

   True_A.sd<-sqrt(True_A_deriv %*% sigma6 %*% t(True_A_deriv))
   True_B.sd<-sqrt(True_B_deriv %*% sigma6 %*% t(True_B_deriv))
   True_C.sd<-sqrt(True_C_deriv %*% sigma6 %*% t(True_C_deriv))

   f.Frac_A<-quote(quote((k21-alpha)*(k31-alpha)/(alpha-beta)/(alpha-gamma)))

   Frac_A<-eval(eval(f.Frac_A))
   ff.Frac_A<-stats::as.formula(paste("~",as.character(f.Frac_A[2],"")))
   f.Frac_B<-quote(quote((k21-beta)*(k31-beta)/(beta-alpha)/(beta-gamma)))
   Frac_B<-eval(eval(f.Frac_B))
   ff.Frac_B<-stats::as.formula(paste("~",as.character(f.Frac_B[2],"")))
   f.Frac_C<-quote(quote((k21-gamma)*(k31-gamma)/(gamma-beta)/(gamma-alpha)))
   Frac_C<-eval(eval(f.Frac_C))
   ff.Frac_C<-stats::as.formula(paste("~",as.character(f.Frac_C[2],"")))

   Frac_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.Frac_A,
        c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))
   Frac_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.Frac_B,
        c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))
   Frac_C_deriv<-as.matrix(attr(eval(stats::deriv(ff.Frac_C,
        c("V1","alpha","beta","gamma","k21","k31"))),"gradient"))

   Frac_A.sd<-sqrt(Frac_A_deriv %*% sigma6 %*% t(Frac_A_deriv))
   Frac_B.sd<-sqrt(Frac_B_deriv %*% sigma6 %*% t(Frac_B_deriv))
   Frac_C.sd<-sqrt(Frac_C_deriv %*% sigma6 %*% t(Frac_C_deriv))
   if(is.na(V1[1])){
     param = rep(NA,24)
     sd = rep(NA,24)
   } else{
     param = c(V1,alpha,beta,gamma,k21,k31,V2,V3,Cl1,Cl2,Cl3,k10,k12,k13,
               t_alpha,t_beta,t_gamma,Vdss,True_A,True_B,True_C,
               Frac_A,Frac_B,Frac_C)
     sd = c(V1.sd,alpha.sd,beta.sd,gamma.sd,k21.sd,k31.sd,V2.sd,V3.sd,Cl1.sd,
            Cl2.sd,Cl3.sd,k10.sd,k12.sd,k13.sd,t_alpha.sd,t_beta.sd,t_gamma.sd,
            Vdss.sd,True_A.sd,True_B.sd,True_C.sd,Frac_A.sd,Frac_B.sd,Frac_C.sd)
   }
   result = data.frame(Parameter=c("V1","alpha","beta","gamma","k21","k31",
                       "V2","V3","Cl1","Cl2","Cl3","k10","k12","k13",
                       "t_alpha","t_beta","t_gamma","Vdss","True_A","True_B",
                       "True_C","Frac_A","Frac_B","Frac_C"),
                       Estimate=param, Std.err=sd)
   row.names(result) <- c("V1","alpha","beta","gamma","k21","k31","V2","V3",
                          "Cl1","Cl2","Cl3","k10","k12","k13","t_alpha",
                          "t_beta","t_gamma","Vdss","True_A","True_B",
                          "True_C","Frac_A","Frac_B","Frac_C")
  result<-result[c("Vdss","V1","V2","V3","Cl1","Cl2","Cl3",
                   "k10","k12","k21","k13","k31","alpha","beta","gamma",
                   "t_alpha","t_beta","t_gamma","True_A","True_B","True_C",
                   "Frac_A","Frac_B","Frac_C"),]
  return(result)
}
