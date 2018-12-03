#' Convert pharmacokinetic parameters for two compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distribution (V1),
#' elimination rate constant (k10), and transter rate constants (k12, k21)
#' @usage TwoComp_Volume_RateConstant(V1,k10,k12,k21,
#'               V1.sd=NA,k10.sd=NA,k12.sd=NA,k21.sd=NA,
#'               covar=c(V1k10=NA,V1k12=NA,V1k21=NA,
#'               k10k12=NA,k10k21=NA,k12k21=NA))
#' @param V1 The volume of distribution of compartment 1
#' @param k10 elimination rate constant
#' @param k12 transfer rate constants from compartment 1 to compartment 2
#' @param k21 transfer rate constants from compartment 2 to compartment 1
#' @param V1.sd standard error of V1
#' @param k10.sd standard error of k10
#' @param k12.sd standard error of k12
#' @param k21.sd standard error of k21
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' TwoComp_Volume_RateConstant(V1=5,k10=0.7,k12=0.5,k21=0.05,
#'          V1.sd=0.01,k10.sd=0.002,k12.sd=0.001,k21.sd=0.0005)
TwoComp_Volume_RateConstant<-function(V1,k10,k12,k21,
                  V1.sd=NA,k10.sd=NA,k12.sd=NA,k21.sd=NA,
                  covar=c(V1k10=NA,V1k12=NA,V1k21=NA,
                  k10k12=NA,k10k21=NA,k12k21=NA)){
  if(is.na(covar[1])) covar<-rep(0,6)
  V1.var <- (V1.sd)^2;
  k10.var <- (k10.sd)^2;  k12.var <- (k12.sd)^2;  k21.var <- (k21.sd)^2

  V2<-V1*k12/k21
  Vdss<-V1+(V1*k12/k21)

  sigma3<-matrix(as.numeric(c(V1.var,covar[2],covar[3],
                              covar[2],k12.var,covar[6],
                              covar[3],covar[6],k21.var)),3,3,byrow=T)

  V2_deriv<-attr(eval(stats::deriv(~V1*k12/k21,c("V1","k12","k21"))),
                 "gradient")
  V2.sd<-sqrt(V2_deriv %*% sigma3 %*% t(V2_deriv))
  Vdss_deriv<-attr(eval(stats::deriv(~V1+(V1*k12/k21),c("V1","k12","k21"))),
                   "gradient")
  Vdss.sd<-sqrt(Vdss_deriv %*% sigma3 %*% t(Vdss_deriv))

  Cl1<-V1*k10
  Cl2<-V1*k12

  sigma2<-matrix(as.numeric(c(V1.var,covar[1],covar[1],k10.var)),2,2,byrow=T)
  Cl1_deriv<-as.matrix(attr(eval(stats::deriv(~V1*k10,c("V1","k10"))),
                            "gradient"))
  Cl1.sd<-sqrt(Cl1_deriv %*% sigma2 %*% t(Cl1_deriv))

  sigma22<-matrix(as.numeric(c(V1.var,covar[2],covar[2],k21.var)),2,2,byrow=T)
  Cl2_deriv<-as.matrix(attr(eval(stats::deriv(~V1*k12,c("V1","k12"))),
                            "gradient"))
  Cl2.sd<-sqrt(Cl2_deriv %*% sigma22 %*% t(Cl2_deriv))

  f.t_alpha<-quote(quote(log(2)/((-(-k10-k12-k21)+sqrt(((-k10-k12-k21))^2-
                                                      4*(k10*k21)))/2)))
  t_alpha<-eval(eval(f.t_alpha))
  ff.t_alpha<-stats::as.formula(paste("~",as.character(f.t_alpha[2],"")))

  f.t_beta<-quote(quote(log(2)/((-(-k10-k12-k21)-sqrt(((-k10-k12-k21))^2-
                                                      4*(k10*k21)))/2)))
  t_beta<-eval(eval(f.t_beta))
  ff.t_beta<-stats::as.formula(paste("~",as.character(f.t_beta[2],"")))

  sigma33<-matrix(as.numeric(c(k10.var,covar[4],covar[5],
                               covar[4],k12.var,covar[6],
                               covar[5],covar[6],k21.var)),3,3,byrow=T)

  t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(ff.t_alpha,
                                         c("k10","k12","k21"))),"gradient"))
  t_alpha.sd<-sqrt(t_alpha_deriv %*% sigma33 %*% t(t_alpha_deriv))

  t_beta_deriv<-as.matrix(attr(eval(stats::deriv(ff.t_beta,
                                          c("k10","k12","k21"))),"gradient"))
  t_beta.sd<-sqrt(t_beta_deriv %*% sigma33 %*% t(t_beta_deriv))

  f.true_A<-quote(quote((k21-((-(-k10-k12-k21)+sqrt(((-k10-k12-k21))^2-
              4*(k10*k21)))/2))/(((-(-k10-k12-k21)-sqrt(((-k10-k12-k21))^2-
              4*(k10*k21)))/2)-((-(-k10-k12-k21)+sqrt(((-k10-k12-k21))^2-
              4*(k10*k21)))/2))/V1))
  true_A<-eval(eval(f.true_A))
  ff.true_A<-stats::as.formula(paste("~",as.character(f.true_A[2],"")))

  f.true_B<-quote(quote((k21-((-(-k10-k12-k21)-sqrt(((-k10-k12-k21))^2-
             4*(k10*k21)))/2))/(((-(-k10-k12-k21)+sqrt(((-k10-k12-k21))^2-
             4*(k10*k21)))/2)-((-(-k10-k12-k21)-sqrt(((-k10-k12-k21))^2-
             4*(k10*k21)))/2))/V1))
  true_B<-eval(eval(f.true_B))
  ff.true_B<-stats::as.formula(paste("~",as.character(f.true_B[2],"")))

  sigma4<-matrix(as.numeric(c(V1.var,covar[1],covar[2],covar[3],
                             covar[4],k10.var,covar[4],covar[5],
                             covar[2],covar[4],k12.var,covar[6],
                             covar[3],covar[5],covar[6],k21.var)),4,4,byrow=T)

  true_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_A,
                   c("V1","k10","k12","k21"))),"gradient"))
  true_A.sd<-sqrt(true_A_deriv %*% sigma4 %*% t(true_A_deriv))
  true_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_B,
                   c("V1","k10","k12","k21"))),"gradient"))
  true_B.sd<-sqrt(true_B_deriv %*% sigma4 %*% t(true_B_deriv))

  f.frac_A<-quote(quote((k21-((-(-k10-k12-k21)+sqrt(((-k10-k12-k21))^2-
            4*(k10*k21)))/2))/(((-(-k10-k12-k21)-sqrt(((-k10-k12-k21))^2-
            4*(k10*k21)))/2)-((-(-k10-k12-k21)+sqrt(((-k10-k12-k21))^2-
            4*(k10*k21)))/2))))
  frac_A<-eval(eval(f.frac_A))
  ff.frac_A<-stats::as.formula(paste("~",as.character(f.frac_A[2],"")))
  f.frac_B<-quote(quote((k21-((-(-k10-k12-k21)-sqrt(((-k10-k12-k21))^2-
                  4*(k10*k21)))/2))/(((-(-k10-k12-k21)+sqrt(((-k10-k12-k21))^2-
                  4*(k10*k21)))/2)-((-(-k10-k12-k21)-sqrt(((-k10-k12-k21))^2-
                  4*(k10*k21)))/2))))
  frac_B<-eval(eval(f.frac_B))
  ff.frac_B<-stats::as.formula(paste("~",as.character(f.frac_B[2],"")))

  frac_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_A,
    c("k10","k12","k21"))),"gradient"))
  frac_A.sd<-sqrt(frac_A_deriv %*% sigma33 %*% t(frac_A_deriv))
  frac_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_B,
    c("k10","k12","k21"))),"gradient"))
  frac_B.sd<-sqrt(frac_B_deriv %*% sigma33 %*% t(frac_B_deriv))


  f.alpha<-quote(quote((k10+k12+k21+sqrt((k10+k12+k21)^2-4*k10*k21))/2))
  alpha<-eval(eval(f.alpha))
  ff.alpha<-stats::as.formula(paste("~",as.character(f.alpha[2],"")))
  f.beta<-quote(quote((k10+k12+k21-sqrt((k10+k12+k21)^2-4*k10*k21))/2))
  beta<-eval(eval(f.beta))
  ff.beta<-stats::as.formula(paste("~",as.character(f.beta[2],"")))

  alpha_deriv<-as.matrix(attr(eval(stats::deriv(ff.alpha,
    c("k10","k12","k21"))),"gradient"))
  alpha.sd<-sqrt(alpha_deriv %*% sigma33 %*% t(alpha_deriv))
  beta_deriv<-as.matrix(attr(eval(stats::deriv(ff.beta,
    c("k10","k12","k21"))),"gradient"))
  beta.sd<-sqrt(beta_deriv %*% sigma33 %*% t(beta_deriv))
  if(is.na(V1[1])){
    param = rep(NA,16)
    sd = rep(NA,16)
  } else{
  param = c(V1,k10,k12,k21,V2,Cl1,Cl2,t_alpha,t_beta,Vdss,true_A,true_B,
    frac_A,frac_B,alpha,beta)
  sd = c(V1.sd,k10.sd,k12.sd,k21.sd,V2.sd,Cl1.sd,Cl2.sd,t_alpha.sd,t_beta.sd,
    Vdss.sd,true_A.sd,true_B.sd,frac_A.sd,frac_B.sd,alpha.sd,beta.sd)
  }
  result = data.frame(Parameter=c("V1","k10","k12","k21","V2","Cl1","Cl2",
                      "t_alpha","t_beta","Vdss","True_A","True_B","Frac_A",
                      "Frac_B","alpha","beta"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("V1","k10","k12","k21","V2","Cl1","Cl2","t_alpha",
    "t_beta","Vdss","True_A","True_B","Frac_A","Frac_B","alpha","beta")
  result<-result[c("Vdss","V1","V2", "Cl1","Cl2","k10","k12","k21",
                   "alpha","beta","t_alpha","t_beta",
                   "True_A","True_B","Frac_A","Frac_B"),]
  return(result)

}
