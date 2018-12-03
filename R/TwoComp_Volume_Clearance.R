#' Convert pharmacokinetic parameters for two compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distributions (V1 and V2)
#'  and clearances (Cl1 and Cl2)
#' @usage TwoComp_Volume_Clearance(V1,V2,Cl1,Cl2,
#'                   V1.sd=NA,V2.sd=NA,Cl1.sd=NA,Cl2.sd=NA,
#'                   covar=c(V1V2=0,V1Cl1=0,V1Cl2=0,
#'                     V2Cl1=0,V2Cl2=0,Cl1Cl2=0))
#' @param V1 The volume of distribution of compartment 1
#' @param V2 The volume of distribution of compartment 2
#' @param Cl1 Clearance from compartment 1
#' @param Cl2 Clearance from compartment 2
#' @param V1.sd standard error of V1
#' @param V2.sd standard error of V2
#' @param Cl1.sd standard error of Cl1
#' @param Cl2.sd standard error of Cl2
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' TwoComp_Volume_Clearance(V1=5,V2=50,Cl1=3.5,Cl2=2.5,
#'          V1.sd=0.01,V2.sd=0.1,Cl1.sd=0.01,Cl2.sd=0.01)

TwoComp_Volume_Clearance<-function(V1,V2,Cl1,Cl2,
                   V1.sd=NA,V2.sd=NA,Cl1.sd=NA,Cl2.sd=NA,
                   covar=c(V1V2=0,V1Cl1=0,V1Cl2=0,
                     V2Cl1=0,V2Cl2=0,Cl1Cl2=0)){
  if(is.na(covar[1])) covar<-rep(0,6)
  Vdss <- V1+V2
  V1.var <- (V1.sd)^2;    V2.var <- (V2.sd)^2
  Cl1.var <- (Cl1.sd)^2;  Cl2.var <- (Cl2.sd)^2

  Vd_sig<-matrix(as.numeric(c(V1.var,covar[1],covar[1],V2.var)),2,2,byrow=T)
  Vd_deriv<-attr(eval(stats::deriv(~V1+V2,c("V1","V2"))),"gradient")
  Vdss.sd<-sqrt(Vd_deriv %*% Vd_sig %*% t(Vd_deriv))

  k10 <- Cl1/V1
  k12 <- Cl2/V1
  k21 <- Cl2/V2

  k10_sig<-matrix(as.numeric(c(V1.var,covar[2],covar[2],Cl1.var)),2,2,byrow=T)
  k10_deriv<-as.matrix(attr(eval(stats::deriv(~Cl1/V1,c("V1","Cl1"))),
                            "gradient"))
  k10.sd<-sqrt(k10_deriv %*% k10_sig %*% t(k10_deriv))

  k12_sig<-matrix(as.numeric(c(V1.var,covar[3],covar[3],Cl2.var)),2,2,byrow=T)
  k12_deriv<-as.matrix(attr(eval(stats::deriv(~Cl2/V1,c("V1","Cl2"))),
                            "gradient"))
  k12.sd<-sqrt(k12_deriv %*% k12_sig %*% t(k12_deriv))

  k21_sig<-matrix(as.numeric(c(Cl2.var,covar[5],covar[5],V2.var)),2,2,byrow=T)
  k21_deriv<-as.matrix(attr(eval(stats::deriv(~Cl2/V2,c("Cl2","V2"))),
                            "gradient"))
  k21.sd<-sqrt(k21_deriv %*% k21_sig %*% t(k21_deriv))

  sigma<-matrix(as.numeric(c(V1.var,covar[1],covar[2],covar[3],
                            covar[1],V2.var,covar[4],covar[5],
                            covar[2],covar[4],Cl1.var,covar[6],
                            covar[3],covar[5],covar[6],Cl2.var)),4,4,byrow=T)

  f.t_alpha<-quote(quote((log(2)*2)/((Cl1+Cl2)/V1+Cl2/V2+
             sqrt(((Cl1+Cl2)/V1+Cl2/V2)^2-4*(Cl1/V1)*(Cl2/V2)))))
  t_alpha<-eval(eval(f.t_alpha))
  ff.t_alpha<-stats::as.formula(paste("~",as.character(f.t_alpha[2],"")))
  t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(ff.t_alpha,
                                    c("V1","V2","Cl1","Cl2"))),"gradient"))
  t_alpha.sd<-sqrt(t_alpha_deriv %*% sigma %*% t(t_alpha_deriv))

  f.t_beta<-quote(quote(log(2)/((-(-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))-
            sqrt((-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))^2-4*((Cl1/V1)*(Cl2/V2))))/2)))
  t_beta<-eval(eval(f.t_beta))
  ff.t_beta<-stats::as.formula(paste("~",as.character(f.t_beta[2],"")))
  t_beta_deriv<-as.matrix(attr(eval(stats::deriv(ff.t_beta,
                                    c("V1","V2","Cl1","Cl2"))),"gradient"))
  t_beta.sd<-sqrt(t_beta_deriv %*% sigma %*% t(t_beta_deriv))

  f.true_A<-quote(quote(((Cl2/V2)-((Cl1+Cl2)/V1+Cl2/V2+
            sqrt(((Cl1+Cl2)/V1+Cl2/V2)^2-
            4*(Cl1/V1)*(Cl2/V2)))/2)/(-sqrt(((Cl1+Cl2)/V1+Cl2/V2)^2-
            4*(Cl1/V1)*(Cl2/V2)))/V1))
  true_A<-eval(eval(f.true_A))
  ff.true_A<-stats::as.formula(paste("~",as.character(f.true_A[2],"")))
  true_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_A,
                                    c("V1","V2","Cl1","Cl2"))),"gradient"))
  true_A.sd<-sqrt(true_A_deriv %*% sigma %*% t(true_A_deriv))

  f.true_B<-quote(quote(((Cl2/V2)-((-(-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))-
           sqrt((-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))^2-4*((Cl1/V1)*(Cl2/V2))))/2))/
           (((-(-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))+sqrt((-((Cl1/V1)+(Cl2/V1)+
           (Cl2/V2)))^2-4*((Cl1/V1)*(Cl2/V2))))/2)-((-(-((Cl1/V1)+
           (Cl2/V1)+(Cl2/V2)))-sqrt((-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))^2-
           4*((Cl1/V1)*(Cl2/V2))))/2))/V1))
  true_B<-eval(eval(f.true_B))
  ff.true_B<-stats::as.formula(paste("~",as.character(f.true_B[2],"")))
  true_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_B,
                                   c("V1","V2","Cl1","Cl2"))),"gradient"))
  true_B.sd<-sqrt(true_B_deriv %*% sigma %*% t(true_B_deriv))

  f.frac_A<-quote(quote(((Cl2/V2)-((Cl1+Cl2)/V1+Cl2/V2+sqrt(((Cl1+Cl2)/V1+
      Cl2/V2)^2-4*(Cl1/V1)*(Cl2/V2)))/2)/(-sqrt(((Cl1+Cl2)/V1+Cl2/V2)^2-
      4*(Cl1/V1)*(Cl2/V2)))))
  frac_A<-eval(eval(f.frac_A))
  ff.frac_A<-stats::as.formula(paste("~",as.character(f.frac_A[2],"")))
  frac_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_A,
                                   c("V1","V2","Cl1","Cl2"))),"gradient"))
  frac_A.sd<-sqrt(frac_A_deriv %*% sigma %*% t(frac_A_deriv))

  f.frac_B<-quote(quote(((Cl2/V2)-((-(-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))-
          sqrt((-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))^2-4*((Cl1/V1)*(Cl2/V2))))/2))/
          (((-(-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))+sqrt((-((Cl1/V1)+(Cl2/V1)+
          (Cl2/V2)))^2-4*((Cl1/V1)*(Cl2/V2))))/2)-
          ((-(-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))-
          sqrt((-((Cl1/V1)+(Cl2/V1)+(Cl2/V2)))^2-4*((Cl1/V1)*(Cl2/V2))))/2))))
  frac_B<-eval(eval(f.frac_B))
  ff.frac_B<-stats::as.formula(paste("~",as.character(f.frac_B[2],"")))
  frac_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_B,
                                  c("V1","V2","Cl1","Cl2"))),"gradient"))
  frac_B.sd<-sqrt(frac_B_deriv %*% sigma %*% t(frac_B_deriv))

  f.alpha<-quote(quote(((Cl1+Cl2)/V1+Cl2/V2+sqrt(((Cl1+Cl2)/V1+Cl2/V2)^2-
           4*(Cl1/V1)*(Cl2/V2)))/2))
  alpha<-eval(eval(f.alpha))
  ff.alpha<-stats::as.formula(paste("~",as.character(f.alpha[2],"")))
  alpha_deriv<-as.matrix(attr(eval(stats::deriv(ff.alpha,
                                 c("V1","V2","Cl1","Cl2"))),"gradient"))
  alpha.sd<-sqrt(alpha_deriv %*% sigma %*% t(alpha_deriv))

  f.beta<-quote(quote(((Cl1+Cl2)/V1+Cl2/V2-sqrt(((Cl1+Cl2)/V1+Cl2/V2)^2-
           4*(Cl1/V1)*(Cl2/V2)))/2))
  beta<-eval(eval(f.beta))
  ff.beta<-stats::as.formula(paste("~",as.character(f.beta[2],"")))
  beta_deriv<-as.matrix(attr(eval(stats::deriv(ff.beta,
                                 c("V1","V2","Cl1","Cl2"))),"gradient"))
  beta.sd<-sqrt(beta_deriv %*% sigma %*% t(beta_deriv))

  if(is.na(V1[1])){
    param = rep(NA,16)
    sd = rep(NA,16)
  } else{
    param = c(V1,V2,Cl1,Cl2,k10,k12,k21,t_alpha,t_beta,Vdss,true_A,true_B,
              frac_A,frac_B,alpha,beta)
    sd = c(V1.sd,V2.sd,Cl1.sd,Cl2.sd,k10.sd,k12.sd,k21.sd,t_alpha.sd,t_beta.sd,
           Vdss.sd,true_A.sd,true_B.sd,frac_A.sd,frac_B.sd,alpha.sd,beta.sd)
  }


  result = data.frame(Parameter=c("V1","V2","Cl1","Cl2","k10","k12","k21",
                      "t_alpha","t_beta","Vdss","True_A","True_B","Frac_A",
                      "Frac_B","alpha","beta"),
                      Estimate=param,Std.err=sd)
  row.names(result) <- c("V1","V2","Cl1","Cl2","k10","k12","k21",
    "t_alpha","t_beta","Vdss","True_A","True_B","Frac_A","Frac_B",
    "alpha","beta")
  result<-result[c("Vdss","V1","V2","Cl1","Cl2","k10","k12","k21",
    "alpha","beta","t_alpha","t_beta","True_A","True_B","Frac_A","Frac_B"),]

  return(result)
}
