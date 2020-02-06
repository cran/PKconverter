#' Convert pharmacokinetic parameters for three compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distributions
#' (V1, V2 and V3) and  clearances (Cl1, Cl2, and Cl3)
#' @usage ThreeComp_Volume_Clearance(V1,V2,V3,Cl1,Cl2,Cl3,
#'  V1.sd=NA,V2.sd=NA,V3.sd=NA, Cl1.sd=NA,Cl2.sd=NA,Cl3.sd=NA,
#'  covar=c(V1V2=NA,V1V3=NA,V1Cl1=NA,
#'    V1Cl2=NA,V1Cl3=NA,V2V3=NA,V2Cl1=NA,V2Cl2=NA,V2Cl3=NA,
#'    V3Cl1=NA,V3Cl2=NA,V3Cl3=NA,Cl1Cl2=NA,Cl1Cl3=NA,Cl2Cl3=NA),...)
#' @param V1 The volume of distribution of compartment 1
#' @param V2 The volume of distribution of compartment 2
#' @param V3 The volume of distribution of compartment 3
#' @param Cl1 Clearance from compartment 1
#' @param Cl2 Clearance from compartment 2
#' @param Cl3 Clearance from compartment 3
#' @param V1.sd standard error of V1
#' @param V2.sd standard error of V2
#' @param V3.sd standard error of V3
#' @param Cl1.sd standard error of Cl1
#' @param Cl2.sd standard error of Cl2
#' @param Cl3.sd standard error of Cl3
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' ThreeComp_Volume_Clearance(V1=10,V2=100,V3=1000,Cl1=3,Cl2=2,Cl3=1,
#'   V1.sd=0.01,V2.sd=0.1,V3.sd=1,Cl1.sd=0.01,Cl2.sd=0.01,Cl3.sd=0.01)
#'
ThreeComp_Volume_Clearance<-function(V1,V2,V3,Cl1,Cl2,Cl3,
  V1.sd=NA,V2.sd=NA,V3.sd=NA,Cl1.sd=NA,Cl2.sd=NA,Cl3.sd=NA,
  covar=c(V1V2=NA,V1V3=NA,V1Cl1=NA,
    V1Cl2=NA,V1Cl3=NA,V2V3=NA,V2Cl1=NA,V2Cl2=NA,V2Cl3=NA,
    V3Cl1=NA,V3Cl2=NA,V3Cl3=NA,Cl1Cl2=NA,Cl1Cl3=NA,Cl2Cl3=NA),...){
  if(is.na(covar[1])) covar<-rep(0,15)
  V1.var <- (V1.sd)^2;    V2.var <- (V2.sd)^2;   V3.var <- (V3.sd)^2
  Cl1.var <- (Cl1.sd)^2;  Cl2.var <- (Cl2.sd)^2; Cl3.var <- (Cl3.sd)^2

  Vdss <- V1+V2+V3
  Vd_sig<-matrix(as.numeric(c(V1.var,covar[1],covar[2],
                              covar[1],V2.var,covar[6],
                              covar[2],covar[6],V3.var)),3,3,byrow=T)
  Vd_deriv<-attr(eval(stats::deriv(~V1+V2+V3,c("V1","V2","V3"))),"gradient")
  Vdss.sd<-sqrt(Vd_deriv %*% Vd_sig %*% t(Vd_deriv))

  k10 <- Cl1/V1
  k12 <- Cl2/V1
  k13 <- Cl3/V1
  k21 <- Cl2/V2
  k31 <- Cl3/V3

  k10_sig<-matrix(as.numeric(c(V1.var,covar[3],covar[3],Cl1.var)),2,2,byrow=T)
  k10_deriv<-as.matrix(attr(eval(stats::deriv(~Cl1/V1,c("V1","Cl1"))),
                            "gradient"))
  k10.sd<-sqrt(k10_deriv %*% k10_sig %*% t(k10_deriv))

  k12_sig<-matrix(as.numeric(c(V1.var,covar[4],covar[4],Cl2.var)),2,2,byrow=T)
  k12_deriv<-as.matrix(attr(eval(stats::deriv(~Cl2/V1,c("V1","Cl2"))),
                            "gradient"))
  k12.sd<-sqrt(k12_deriv %*% k12_sig %*% t(k12_deriv))

  k13_sig<-matrix(as.numeric(c(V1.var,covar[5],covar[5],Cl3.var)),2,2,byrow=T)
  k13_deriv<-as.matrix(attr(eval(stats::deriv(~Cl3/V1,c("V1","Cl3"))),
                            "gradient"))
  k13.sd<-sqrt(k13_deriv %*% k13_sig %*% t(k13_deriv))

  k21_sig<-matrix(as.numeric(c(Cl2.var,covar[5],covar[5],V2.var)),2,2,byrow=T)
  k21_deriv<-as.matrix(attr(eval(stats::deriv(~Cl2/V2,c("Cl2","V2"))),
                            "gradient"))
  k21.sd<-sqrt(k21_deriv %*% k21_sig %*% t(k21_deriv))

  k31_sig<-matrix(as.numeric(c(Cl3.var,covar[12],covar[12],V3.var)),2,2,byrow=T)
  k31_deriv<-as.matrix(attr(eval(stats::deriv(~Cl3/V3,c("Cl3","V3"))),
                            "gradient"))
  k31.sd<-sqrt(k31_deriv %*% k31_sig %*% t(k31_deriv))

  root1<-(-(cos(acos((-((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/
        (sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
        (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
        (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
        (Cl2/V2)+(Cl3/V3))/3)))
  root2<-(-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/
        (sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
        (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
        (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*
        (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))
  root3<-(-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/
        (sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
        (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
        (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*(2*
        (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))
  root<-c(root1,root2,root3)
  if(is.na(root[1])){
    l1<-l2<-l3<-NA
  }else{
  l1<-max(root);  l2<-stats::median(root);   l3<-min(root)
  }
  alpha<-l1;      beta<-l2;           gamma<-l3
  t_alpha<-log(2)/l1; t_beta<-log(2)/l2;  t_gamma<-log(2)/l3

  f.t_alpha<-quote(quote(log(2)/(-(cos(acos((-((2*
      (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
      (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
      (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
      ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
      (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
      ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*
      (2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
      (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
      (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
      (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))
  ff.t_alpha<-stats::as.formula(paste("~",as.character(f.t_alpha[2],"")))
  t_alpha_deriv<-as.matrix(attr(eval(stats::deriv(ff.t_alpha,
                      c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  f.t_beta<-quote(quote(log(2)/(-(cos((acos((-((2*
      (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
      (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
      (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
      ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
      (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
      ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+
      2*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
      (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
      (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+
      (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))
  ff.t_beta<-stats::as.formula(paste("~",as.character(f.t_beta[2],"")))
  t_beta_deriv<-as.matrix(attr(eval(stats::deriv(ff.t_beta,
                     c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  f.t_gamma<-quote(quote(log(2)/(-(cos((acos((-((2*
      (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
      (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
      (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
      ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
      (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
      ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+
      4*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)
      +(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
      ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-
      (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))
  ff.t_gamma<-stats::as.formula(paste("~",as.character(f.t_gamma[2],"")))
  t_gamma_deriv<-as.matrix(attr(eval(stats::deriv(ff.t_gamma,
                       c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  sigma6<-matrix(as.numeric(c(V1.var,covar[1],covar[2],covar[3],covar[4],
                 covar[5],covar[1],V2.var,covar[6],covar[7],covar[8],covar[9],
                 covar[2],covar[6],V3.var,covar[9],covar[11],covar[12],
                 covar[3],covar[7],covar[10],Cl1.var,covar[13],covar[14],
                 covar[4],covar[8],covar[11],covar[13],Cl2.var,covar[15],
                 covar[5],covar[9],covar[12],covar[14],covar[15],Cl3.var)),
                 6,6,byrow=T)

  t_alpha.s<-sqrt(t_alpha_deriv %*% sigma6 %*% t(t_alpha_deriv))
  t_beta.s<-sqrt(t_beta_deriv %*% sigma6 %*% t(t_beta_deriv))
  t_gamma.s<-sqrt(t_gamma_deriv %*% sigma6 %*% t(t_gamma_deriv))
  half_sig<-c(t_alpha.s,t_beta.s,t_gamma.s)
  rt_ord<-order(root)
  t_alpha.sd<-half_sig[which(rt_ord==3)];
  t_beta.sd<-half_sig[which(rt_ord==2)];
  t_gamma.sd<-half_sig[which(rt_ord==1)]

  true_A<-(((Cl2/V2)-(l1))*((Cl3/V3)-(l1))/((l1)-(l2))/((l1)-(l3))/(V1))
  true_B<-(((Cl2/V2)-(l2))*((Cl3/V3)-(l2))/((l2)-(l1))/((l2)-(l3))/(V1))
  true_C<-(((Cl2/V2)-(l3))*((Cl3/V3)-(l3))/((l3)-(l2))/((l3)-(l1))/(V1))


  f.true_A<-quote(quote((((Cl2/V2)-((-(cos(acos((-((2*
    (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
    (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*
    (Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3)))))*((Cl3/V3)-((-(cos(acos((-((2*(((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*
    ((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
    (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos(acos((-((2*(((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
    (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*
    (Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3))))-((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*
    ((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
    (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos(acos((-((2*(((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
    (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3))))-((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
    (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3)))))/(V1))))

  ff.true_A<-stats::as.formula(paste("~",as.character(f.true_A[2],"")))
  true_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_A,
                   c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))
  f.true_B<-quote(quote((((Cl2/V2)-((-(cos((acos((-((2*
    (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
    (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*
    (Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3)))))*((Cl3/V3)-((-(cos((acos((-((2*(((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*
    ((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
    (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos((acos((-((2*(((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*
    ((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
    (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3))))-((-(cos(acos((-((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/
    (sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
    (Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/
    27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))/
    (((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
    (Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/
    (sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3))))-((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
    (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3)))))/(V1))))
  ff.true_B<-stats::as.formula(paste("~",as.character(f.true_B[2],"")))
  true_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_B,
                  c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  f.true_C<-quote(quote((((Cl2/V2)-((-(cos((acos((-((2*
    (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
    (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*
    (Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*
    (2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))*((Cl3/V3)-((-(cos((acos((-((2*(((Cl1/
    V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*
    ((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
    (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*
    ((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
    (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*(2*
    (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3))))-((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+
    (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/
    (sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*
    (Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
    (Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
    (Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/
    27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))/(((
    -(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-
    (((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*
    (Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
    (Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
    (Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/
    27)))/3)+4*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
    ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-
    (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3))))-((-(cos(acos((-
    ((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
    (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
    (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*
    (Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
    (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
    (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*(2*(
    exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
    (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
    (Cl2/V2)+(Cl3/V3))/3)))))/(V1))))

  ff.true_C<-stats::as.formula(paste("~",as.character(f.true_C[2],"")))
  true_C_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_C,
                  c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  true_sig1<-sqrt(true_A_deriv %*% sigma6 %*% t(true_A_deriv))
  true_sig2<-sqrt(true_B_deriv %*% sigma6 %*% t(true_B_deriv))
  true_sig3<-sqrt(true_C_deriv %*% sigma6 %*% t(true_C_deriv))
  true_sig<-c(true_sig1,true_sig2,true_sig3)

  true_A.sd<-true_sig[which(rt_ord==3)];
  true_B.sd<-true_sig[which(rt_ord==2)];
  true_C.sd<-true_sig[which(rt_ord==1)]

  frac_A<-(((Cl2/V2)-(l1))*((Cl3/V3)-(l1))/((l1)-(l2))/((l1)-(l3)))
  frac_B<-(((Cl2/V2)-(l2))*((Cl3/V3)-(l2))/((l2)-(l1))/((l2)-(l3)))
  frac_C<-(((Cl2/V2)-(l3))*((Cl3/V3)-(l3))/((l3)-(l2))/((l3)-(l1)))

  f.frac_A<-quote(quote(((Cl2/V2)-((-(cos(acos((-((2*
        (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
        ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*
        (2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))*((Cl3/V3)-((-(cos(acos((-((2*
        (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
        ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*
        (2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos(acos((-((2*(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/
        27)))/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3))))-((-(cos((acos((-((2*
        (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
        ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+
        2*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos(acos((-((2*
        (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
        ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*
        (2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3))))-((-(cos((acos((-((2*(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*
        (Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*
        (2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))))

  ff.frac_A<-stats::as.formula(paste("~",as.character(f.frac_A[2],"")))
  frac_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_A,
                           c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  f.frac_B<-quote(quote(((Cl2/V2)-((-(cos((acos((-((2*
        (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
        ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+
        2*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-
        (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))*((Cl3/V3)-
        ((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/
        (sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
        (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
        (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*
        (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos((acos((-((2*(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/
        27)))/3)+2*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
        ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-
        (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3))))-
        ((-(cos(acos((-((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/
        27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*
        (Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^2)/3))^3)/27)))/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
        (Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*
        ((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
        (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*
        (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3))))-((-(cos((acos((-((2*(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
        (Cl2/V2)+(Cl3/V3))/3)))))))
   ff.frac_B<-stats::as.formula(paste("~",as.character(f.frac_B[2],"")))
   frac_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_B,
                         c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

   f.frac_C<-quote(quote(((Cl2/V2)-((-(cos((acos((-((2*
        (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
        ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+
        4*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))*((Cl3/V3)-
        ((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/
        (sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+
        (Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
        (Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*(2*
        (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos((acos((-((2*(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+
        (Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+
        (Cl2/V2)+(Cl3/V3))/3))))-((-(cos((acos((-((2*(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))*
        ((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*(Cl2/V2)*
        (Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*(2*
        (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))/(((-(cos((acos((-((2*(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
        ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+
        4*pi/3)*(2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3))))-((-(cos(acos((-((2*
        (((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+
        (Cl3/V3)*(Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+
        ((Cl1/V1)*(Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-
        ((((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*
        (2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)))))))
  ff.frac_C<-stats::as.formula(paste("~",as.character(f.frac_C[2],"")))
  frac_C_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_C,
                         c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  frac_sig1<-sqrt(frac_A_deriv %*% sigma6 %*% t(frac_A_deriv))
  frac_sig2<-sqrt(frac_B_deriv %*% sigma6 %*% t(frac_B_deriv))
  frac_sig3<-sqrt(frac_C_deriv %*% sigma6 %*% t(frac_C_deriv))
  frac_sig<-c(frac_sig1,frac_sig2,frac_sig3)

  frac_A.sd<-frac_sig[which(rt_ord==3)];
  frac_B.sd<-frac_sig[which(rt_ord==2)];
  frac_C.sd<-frac_sig[which(rt_ord==1)]

  f.alpha<-quote(quote(-(cos(acos((-((2*(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*
        (Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)*(2*
        (exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3))))
  ff.alpha<-stats::as.formula(paste("~",as.character(f.alpha[2],"")))
  alpha_deriv<-as.matrix(attr(eval(stats::deriv(ff.alpha,
                     c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  f.beta<-quote(quote(-(cos((acos((-((2*(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*
        (Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+2*pi/3)*
        (2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3))))
  ff.beta<-stats::as.formula(paste("~",as.character(f.beta[2],"")))
  beta_deriv<-as.matrix(attr(eval(stats::deriv(ff.beta,
                      c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  f.gamma<-quote(quote(-(cos((acos((-((2*(((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^3)/27)-(((Cl1/V1)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V3)+(Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*
        (Cl2/V1))*((Cl1/V1)+(Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3)+((Cl1/V1)*
        (Cl2/V2)*(Cl3/V3)))/2)/(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+
        (Cl2/V2)*(Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+
        (Cl2/V1)+(Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27)))/3)+4*pi/3)*
        (2*(exp(log(sqrt(-((((Cl1/V1)*(Cl3/V3)+(Cl2/V2)*(Cl3/V3)+(Cl2/V2)*
        (Cl3/V1)+(Cl1/V1)*(Cl2/V2)+(Cl3/V3)*(Cl2/V1))-((((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))^2)/3))^3)/27))/3)))-(((Cl1/V1)+(Cl2/V1)+
        (Cl3/V1)+(Cl2/V2)+(Cl3/V3))/3))))
  ff.gamma<-stats::as.formula(paste("~",as.character(f.gamma[2],"")))
  gamma_deriv<-as.matrix(attr(eval(stats::deriv(ff.gamma,
                      c("V1","V2","V3","Cl1","Cl2","Cl3"))),"gradient"))

  exp_sig1<-sqrt(alpha_deriv %*% sigma6 %*% t(alpha_deriv))
  exp_sig2<-sqrt(beta_deriv %*% sigma6 %*% t(beta_deriv))
  exp_sig3<-sqrt(gamma_deriv %*% sigma6 %*% t(gamma_deriv))
  exp_sig<-c(exp_sig1,exp_sig2,exp_sig3)

  alpha.sd<-exp_sig[which(rt_ord==3)];
  beta.sd<-exp_sig[which(rt_ord==2)];
  gamma.sd<-exp_sig[which(rt_ord==1)]
  if(is.na(V1[1])){
    param = rep(NA,24)
    sd = rep(NA,24)
  } else{
    param = c(V1,V2,V3,Cl1,Cl2,Cl3,k10,k12,k13,k21,k31,t_alpha,
             t_beta,t_gamma,Vdss,true_A,true_B,true_C,
             frac_A,frac_B,frac_C,alpha,beta,gamma)
    sd = c(V1.sd,V2.sd,V3.sd,Cl1.sd,Cl2.sd,Cl3.sd,k10.sd,k12.sd,k13.sd,k21.sd,
           k31.sd,t_alpha.sd,t_beta.sd,t_gamma.sd,Vdss.sd,true_A.sd,true_B.sd,
           true_C.sd,frac_A.sd,frac_B.sd,frac_C.sd,alpha.sd,beta.sd,gamma.sd)
  }
  result = data.frame(Parameter=c("V1","V2","V3","Cl1","Cl2","Cl3","k10",
                      "k12","k13","k21","k31","t_alpha","t_beta","t_gamma",
                      "Vdss","True_A","True_B","True_C","Frac_A","Frac_B",
                      "Frac_C","alpha","beta","gamma"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("V1","V2","V3","Cl1","Cl2","Cl3","k10","k12","k13",
                         "k21", "k31","t_alpha","t_beta","t_gamma","Vdss",
                         "True_A","True_B","True_C","Frac_A","Frac_B",
                         "Frac_C","alpha","beta","gamma")
  result<-result[c("Vdss","V1","V2","V3","Cl1","Cl2","Cl3",
                   "k10","k12","k21","k13","k31","alpha","beta","gamma",
                   "t_alpha","t_beta","t_gamma","True_A","True_B","True_C",
                   "Frac_A","Frac_B","Frac_C"),]
  return(result)
}
