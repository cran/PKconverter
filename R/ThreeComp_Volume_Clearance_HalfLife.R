#' Convert pharmacokinetic parameters for three compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distributions (Vd and V1),
#' clearance (Cl1) and half-lives (t_alpha, t_beta, and t_gamma)
#'
#' @usage ThreeComp_Volume_Clearance_HalfLife(V1,Vd,Cl1,t_alpha,t_beta,t_gamma,
#'  V1.sd=NA,Vd.sd=NA,Cl1.sd=NA,t_alpha.sd=NA,t_beta.sd=NA,t_gamma.sd=NA,
#'  covar=c(V1Vd=NA,V1Cl1=NA,V1talpha=NA,V1tbeta=NA,V1tgamma=NA,VdCl1=NA,
#'    Vdtalpha=NA,Vdtbeta=NA,Vdtgamma=NA,Cl1talpha=NA,Cl1tbeta=NA,
#'    Cl1tgamma=NA,talphatbeta=NA,talphatgamma=NA,tbetatgamma=NA))
#' @param Vd Total volume of distributions
#' @param V1 The volume of distribution of compartment 1
#' @param Cl1 Clearance from compartment 1
#' @param t_alpha half life of compartment 1
#' @param t_beta half life of compartment 2
#' @param t_gamma half life of compartment 3
#' @param Vd.sd standard error of Vd
#' @param V1.sd standard error of V1
#' @param Cl1.sd standard error of Cl1
#' @param t_alpha.sd standard error of t_alpha
#' @param t_beta.sd standard error of t_beta
#' @param t_gamma.sd standard error of t_gamma
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' ThreeComp_Volume_Clearance_HalfLife(V1=5,Vd=1110,Cl1=3,
#' t_alpha=1.142,t_beta=52.2,t_gamma=931, V1.sd=0.01,Vd.sd=20,Cl1.sd=0.01,
#' t_alpha.sd=0.002,t_beta.sd=0.5,t_gamma.sd=5.6)

ThreeComp_Volume_Clearance_HalfLife<-function(V1,Vd,Cl1,t_alpha,t_beta,t_gamma,
  V1.sd=NA,Vd.sd=NA,Cl1.sd=NA,t_alpha.sd=NA,t_beta.sd=NA,t_gamma.sd=NA,
  covar=c(V1Vd=NA,V1Cl1=NA,V1talpha=NA,V1tbeta=NA,V1tgamma=NA,VdCl1=NA,
    Vdtalpha=NA,Vdtbeta=NA,Vdtgamma=NA,Cl1talpha=NA,Cl1tbeta=NA,
    Cl1tgamma=NA,talphatbeta=NA,talphatgamma=NA,tbetatgamma=NA)){
  if(is.na(covar[1])) covar<-rep(0,15)
  V1.var = (V1.sd)^2;                  Vd.var = (Vd.sd)^2
  Cl1.var = (Cl1.sd)^2;                t_alpha.var = (t_alpha.sd)^2;
  t_beta.var = (t_beta.sd)^2;          t_gamma.var = (t_gamma.sd)^2

  f.V2<-quote(quote((V1)*((log(2)/t_alpha)+(log(2)/t_beta)+(log(2)/t_gamma)-
      (Cl1/V1)-((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-((((Vd/V1-1)*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-((log(2)/t_alpha)+
        (log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1)))*((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))/
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4))))/((((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))+
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)))

  V2<-eval(eval(f.V2))
  ff.V2<-stats::as.formula(paste("~",as.character(f.V2[2],"")))

  f.V3<-quote(quote((V1)*((((Vd/V1-1)*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-((log(2)/t_alpha)+(log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-
        ((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1)))*((((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))/(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/
        ((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)))

  V3<-eval(eval(f.V3))
  ff.V3<-stats::as.formula(paste("~",as.character(f.V3[2],"")))

  V2_deriv<-as.matrix(attr(eval(stats::deriv(ff.V2,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  V3_deriv<-as.matrix(attr(eval(stats::deriv(ff.V3,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  f.Vdss<-quote(quote(((V1)*((log(2)/t_alpha)+
        (log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))
        -((((Vd/V1-1)*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        ((log(2)/t_alpha)+(log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-
        ((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1)))*((((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))/(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4))))/
        ((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*
        t_gamma)/(Cl1/V1)))/(Cl1/V1))+(sqrt(((((log(2)^2)/(t_alpha*t_beta))+
        ((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))+
        (V1*((((Vd/V1-1)*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        ((log(2)/t_alpha)+(log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-
        ((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1)))*((((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*
        t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1))*4)))/2))/(sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/((((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))+(V1)))

  Vdss<-eval(eval(f.Vdss))
  ff.Vdss<-stats::as.formula(paste("~",as.character(f.Vdss[2],"")))
  Vdss_deriv<-as.matrix(attr(eval(stats::deriv(ff.Vdss,
              c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  sigma6<-matrix(as.numeric(c(V1.var,covar[1],covar[2],covar[3],covar[4],
                 covar[5],covar[1],Vd.var,covar[6],covar[7],covar[8],covar[9],
                 covar[2],covar[6],Cl1.var,covar[9],covar[11],covar[12],
                 covar[3],covar[7],covar[10],t_alpha.var,covar[13],covar[14],
                 covar[4],covar[8],covar[11],covar[13],t_beta.var,covar[15],
                 covar[5],covar[9],covar[12],covar[14],covar[15],t_gamma.var)),
                 6,6,byrow=T)

  V2.sd<-sqrt(V2_deriv %*% sigma6 %*% t(V2_deriv))
  V3.sd<-sqrt(V3_deriv %*% sigma6 %*% t(V3_deriv))
  Vdss.sd<-sqrt(Vdss_deriv %*% sigma6 %*% t(Vdss_deriv))


  f.Cl2<-quote(quote(V1*((log(2)/t_alpha)+(log(2)/t_beta)+(log(2)/t_gamma)-
        (Cl1/V1)-((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-((((Vd/V1-1)*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-((log(2)/t_alpha)+
        (log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1)))*
        ((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))/
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4))))))

  Cl2<-eval(eval(f.Cl2))
  ff.Cl2<-stats::as.formula(paste("~",as.character(f.Cl2[2],"")))

  f.Cl3<-quote(quote(V1*((((Vd/V1-1)*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
      (Cl1/V1)))-((log(2)/t_alpha)+(log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-
        ((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1)))*((((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))/(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))))

  Cl3<-eval(eval(f.Cl3))
  ff.Cl3<-stats::as.formula(paste("~",as.character(f.Cl3[2],"")))

  Cl2_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl2,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))
  Cl2.sd<-sqrt(Cl2_deriv %*% sigma6 %*% t(Cl2_deriv))

  Cl3_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl3,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))
  Cl3.sd<-sqrt(Cl3_deriv %*% sigma6 %*% t(Cl3_deriv))

  k10<-Cl1/V1

  f.k12<-quote(quote((log(2)/t_alpha)+(log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-
        ((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))+(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        ((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        ((((Vd/V1-1)*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        ((log(2)/t_alpha)+(log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-
        ((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1)))*((((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))-(sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))/(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))))

  k12<-eval(eval(f.k12))
  ff.k12<-stats::as.formula(paste("~",as.character(f.k12[2],"")))

  f.k13<-quote(quote((((Vd/V1-1)*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-((log(2)/t_alpha)+(log(2)/t_beta)+(log(2)/t_gamma)-(Cl1/V1)-
        ((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1)))*((((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))/(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4))))

  k13<-eval(eval(f.k13))
  ff.k13<-stats::as.formula(paste("~",as.character(f.k13[2],"")))

  f.k21<-quote(quote((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))+(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))

  k21<-eval(eval(f.k21))

  ff.k21<-stats::as.formula(paste("~",as.character(f.k21[2],"")))

  f.k31<-quote(quote((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2))

  k31<-eval(eval(f.k31))
  ff.k31<-stats::as.formula(paste("~",as.character(f.k31[2],"")))

  sigma2<-matrix(as.numeric(c(V1.var,covar[1],covar[1],Cl1.var)),2,2,byrow=T)

  k10_deriv<-as.matrix(attr(eval(stats::deriv(~Cl1/V1,c("V1","Cl1"))),
                            "gradient"))
  k12_deriv<-as.matrix(attr(eval(stats::deriv(ff.k12,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))
  k13_deriv<-as.matrix(attr(eval(stats::deriv(ff.k13,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))
  k21_deriv<-as.matrix(attr(eval(stats::deriv(ff.k21,
        c("t_alpha","t_beta","t_gamma","V1","Vd","Cl1"))),"gradient"))
  k31_deriv<-as.matrix(attr(eval(stats::deriv(ff.k31,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  k10.sd<-sqrt(k10_deriv %*% sigma2 %*% t(k10_deriv))
  k12.sd<-sqrt(k12_deriv %*% sigma6 %*% t(k12_deriv))
  k13.sd<-sqrt(k13_deriv %*% sigma6 %*% t(k13_deriv))
  k21.sd<-sqrt(k21_deriv %*% sigma6 %*% t(k21_deriv))
  k31.sd<-sqrt(k31_deriv %*% sigma6 %*% t(k31_deriv))


  f.true_A<-quote(quote((((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
      (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
      (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
      (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))+(sqrt(((((log(2)^2)/
      (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
      (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
      (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
      (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
      (log(2)/t_alpha))*(((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
      (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
      (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
      (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-(sqrt(((((log(2)^2)/
      (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
      (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
      (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-
      (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
      (log(2)/t_alpha))/((log(2)/t_alpha)-(log(2)/t_beta))/
      ((log(2)/t_alpha)-(log(2)/t_gamma))/V1))

  true_A<-eval(eval(f.true_A))

  ff.true_A<-stats::as.formula(paste("~",as.character(f.true_A[2],"")))

  f.true_B<-quote(quote((((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))+(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_beta))*(((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_beta))/((log(2)/t_beta)-(log(2)/t_alpha))/((log(2)/t_beta)-
        (log(2)/t_gamma))/V1))

  true_B<-eval(eval(f.true_B))

  ff.true_B<-stats::as.formula(paste("~",as.character(f.true_B[2],"")))

  f.true_C<-quote(quote((((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))+(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_gamma))*(((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-((Vd-V1)/V1*
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-(((log(2))^3)/
        (t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-(sqrt(((((log(2)^2)/
        (t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/
        (t_beta*t_gamma))-((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/
        (Cl1/V1)))-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/
        (Cl1/V1))^2-(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_gamma))/((log(2)/t_gamma)-(log(2)/t_beta))/((log(2)/t_gamma)-
        (log(2)/t_alpha))/V1))

  true_C<-eval(eval(f.true_C))
  ff.true_C<-stats::as.formula(paste("~",as.character(f.true_C[2],"")))

  true_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_A,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  true_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_B,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  true_C_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_C,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  true_A.sd<-sqrt(true_A_deriv %*% sigma6 %*% t(true_A_deriv))
  true_B.sd<-sqrt(true_B_deriv %*% sigma6 %*% t(true_B_deriv))
  true_C.sd<-sqrt(true_C_deriv %*% sigma6 %*% t(true_C_deriv))

  f.frac_A<-quote(quote((((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))+
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_alpha))*(((((((log(2)^2)/(t_alpha*t_beta))+
        ((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_alpha))/((log(2)/t_alpha)-(log(2)/t_beta))/((log(2)/t_alpha)-
        (log(2)/t_gamma))))

  frac_A<-eval(eval(f.frac_A))

  ff.frac_A<-stats::as.formula(paste("~",as.character(f.frac_A[2],"")))
  f.frac_B<-quote(quote((((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))+
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_beta))*(((((((log(2)^2)/(t_alpha*t_beta))+
        ((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_beta))/((log(2)/t_beta)-(log(2)/t_alpha))/((log(2)/t_beta)-
        (log(2)/t_gamma))))

  frac_B<-eval(eval(f.frac_B))
  ff.frac_B<-stats::as.formula(paste("~",as.character(f.frac_B[2],"")))

  f.frac_C<-quote(quote((((((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/
        (t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))+
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_gamma))*(((((((log(2)^2)/(t_alpha*t_beta))+
        ((log(2)^2)/(t_alpha*t_gamma))+((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))-
        (sqrt(((((log(2)^2)/(t_alpha*t_beta))+((log(2)^2)/(t_alpha*t_gamma))+
        ((log(2)^2)/(t_beta*t_gamma))-
        ((Vd-V1)/V1*(((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1)))/(Cl1/V1))^2-
        (((log(2))^3)/(t_alpha*t_beta*t_gamma)/(Cl1/V1))*4)))/2)-
        (log(2)/t_gamma))/((log(2)/t_gamma)-(log(2)/t_beta))/((log(2)/t_gamma)-
        (log(2)/t_alpha))))

  frac_C<-eval(eval(f.frac_C))
  ff.frac_C<-stats::as.formula(paste("~",as.character(f.frac_C[2],"")))

  frac_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_A,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  frac_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_B,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  frac_C_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_C,
        c("V1","Vd","Cl1","t_alpha","t_beta","t_gamma"))),"gradient"))

  frac_A.sd<-sqrt(frac_A_deriv %*% sigma6 %*% t(frac_A_deriv))
  frac_B.sd<-sqrt(frac_B_deriv %*% sigma6 %*% t(frac_B_deriv))
  frac_C.sd<-sqrt(frac_C_deriv %*% sigma6 %*% t(frac_C_deriv))

  alpha<-log(2)/t_alpha; beta<-log(2)/t_beta; gamma<-log(2)/t_gamma

  alpha_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/t_alpha,"t_alpha")),
                               "gradient"))
  beta_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/t_beta,"t_beta")),
                               "gradient"))
  gamma_deriv<-as.matrix(attr(eval(stats::deriv(~log(2)/t_gamma,"t_gamma")),
                               "gradient"))

  alpha.sd<-sqrt(alpha_deriv * t_alpha.var * alpha_deriv)
  beta.sd<-sqrt(beta_deriv * t_beta.var * beta_deriv)
  gamma.sd<-sqrt(gamma_deriv * t_gamma.var * gamma_deriv)
  if(is.na(V1[1])){
    param = rep(NA,24)
    sd = rep(NA,24)
  } else{
    param = c(V1,Vdss,Cl1,t_alpha,t_beta,t_gamma,V2,V3,Cl2,Cl3,k10,k12,k13,
      k21,k31,true_A,true_B,true_C,frac_A,frac_B,frac_C,alpha,beta,gamma)
    sd = c(V1.sd,Vdss.sd,Cl1.sd,t_alpha.sd,t_beta.sd,t_gamma.sd,V2.sd,V3.sd,
      Cl2.sd,Cl3.sd,k10.sd,k12.sd,k13.sd,k21.sd,k31.sd,true_A.sd,true_B.sd,
      true_C.sd,frac_A.sd,frac_B.sd,frac_C.sd,alpha.sd,beta.sd,gamma.sd)
  }
  result = data.frame(Parameter=c("V1","Vdss","Cl1","t_alpha","t_beta",
                      "t_gamma","V2","V3","Cl2","Cl3","k10","k12","k13",
                      "k21","k31","True_A","True_B","True_C","Frac_A",
                      "Frac_B","Frac_C","alpha","beta","gamma"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("V1","Vdss","Cl1","t_alpha","t_beta","t_gamma","V2",
    "V3","Cl2","Cl3","k10","k12","k13","k21","k31","True_A","True_B","True_C",
    "Frac_A","Frac_B","Frac_C","alpha","beta","gamma")
  result<-result[c("Vdss","V1","V2","V3","Cl1","Cl2","Cl3",
                   "k10","k12","k21","k13","k31","alpha","beta","gamma",
                   "t_alpha","t_beta","t_gamma","True_A","True_B","True_C",
                   "Frac_A","Frac_B","Frac_C"),]
  return(result)
}
