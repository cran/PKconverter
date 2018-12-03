#' Convert pharmacokinetic parameters for three compartment model
#'
#' Calculate pharmacokinetic parameters with volume of distribution (V1),
#' elimination rate constant (k10), and transter rate constants
#' (k12, k13, k21, and k31)
#' @usage ThreeComp_Volume_RateConstant(V1,k10,k12,k13,k21,k31,
#'  V1.sd=NA,k10.sd=NA,k12.sd=NA,
#'  k13.sd=NA,k21.sd=NA,k31.sd=NA,covar=c(V1k10=NA,V1k12=NA,V1k13=NA,
#'    V1k21=NA,V1k31=NA,k10k12=NA,k10k13=NA,k10k21=NA,k10k31=NA,
#'    k12k13=NA,k12k21=NA,k12k31=NA,k13k21=NA,k13k31=NA,k21k31=NA))
#' @param V1 The volume of distribution of compartment 1
#' @param k10 elimination rate constant
#' @param k12 transfer rate constants from compartment 1 to compartment 2
#' @param k13 transfer rate constants from compartment 1 to compartment 3
#' @param k21 transfer rate constants from compartment 2 to compartment 1
#' @param k31 transfer rate constants from compartment 3 to compartment 1
#' @param V1.sd standard error of V1
#' @param k10.sd standard error of k10
#' @param k12.sd standard error of k12
#' @param k13.sd standard error of k13
#' @param k21.sd standard error of k21
#' @param k31.sd standard error of k31
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' ThreeComp_Volume_RateConstant(V1=10,k10=0.3,k12=0.2,k13=0.1,k21=0.02,
#'          k31=0.001,V1.sd=0.1,k10.sd=0.002,k12.sd=0.001,
#'          k13.sd=0.0005,k21.sd=0.0005,k31.sd=0.000005)

ThreeComp_Volume_RateConstant<-function(V1,k10,k12,k13,k21,k31,
    V1.sd=NA,k10.sd=NA,k12.sd=NA,
    k13.sd=NA,k21.sd=NA,k31.sd=NA,covar=c(V1k10=NA,V1k12=NA,V1k13=NA,
    V1k21=NA,V1k31=NA,k10k12=NA,k10k13=NA,k10k21=NA,k10k31=NA,
    k12k13=NA,k12k21=NA,k12k31=NA,k13k21=NA,k13k31=NA,k21k31=NA)){
  if(is.na(covar[1])) covar<-rep(0,15)
  V1.var <- (V1.sd)^2;
  k10.var <- (k10.sd)^2;  k12.var <- (k12.sd)^2;  k13.var <- (k13.sd)^2
  k21.var <- (k21.sd)^2;  k31.var <- (k31.sd)^2


  V2<-V1*k12/k21
  V3<-V1*k13/k31
  Vdss<-V1+(V1*k12/k21)+(V1*k13/k31)

  sigma31<-matrix(as.numeric(c(V1.var,covar[2],covar[4],
                             covar[2],k12.var,covar[11],
                             covar[4],covar[11],k21.var)),3,3,byrow=T)
  V2_deriv<-as.matrix(attr(eval(stats::deriv(~V1*k12/k21,c("V1","k12","k21"))),
                           "gradient"))
  V2.sd<-sqrt(V2_deriv %*% sigma31 %*% t(V2_deriv))

  sigma32<-matrix(as.numeric(c(V1.var,covar[3],covar[5],
                               covar[3],k13.var,covar[14],
                               covar[5],covar[14],k31.var)),3,3,byrow=T)
  V3_deriv<-as.matrix(attr(eval(stats::deriv(~V1*k13/k31,c("V1","k13","k31"))),
                           "gradient"))
  V3.sd<-sqrt(V3_deriv %*% sigma32 %*% t(V3_deriv))

  sigma5<-matrix(as.numeric(c(V1.var,covar[2],covar[3],covar[4],covar[5],
                              covar[2],k12.var,covar[10],covar[11],covar[12],
                              covar[3],covar[10],k13.var,covar[13],covar[14],
                              covar[4],covar[8],covar[13],k21.var,covar[15],
                              covar[5],covar[12],covar[13],covar[15],k31.var)),
    5,5,byrow=T)

  Vdss_deriv<-as.matrix(attr(eval(stats::deriv(~V1+(V1*k12/k21)+(V1*k13/k31),
        c("V1","k12","k13","k21","k31"))),"gradient"))
  Vdss.sd<-sqrt(Vdss_deriv %*% sigma5 %*% t(Vdss_deriv))

  Cl1<-V1*k10;   Cl2<-V1*k12;   Cl3<-V1*k13

  sigma21<-matrix(as.numeric(c(V1.var,covar[1],covar[1],k10.var)),2,2,byrow=T)
  Cl1_deriv<-as.matrix(attr(eval(stats::deriv(~V1*k10,c("V1","k10"))),
                            "gradient"))
  Cl1.sd<-sqrt(Cl1_deriv %*% sigma21 %*% t(Cl1_deriv))

  sigma22<-matrix(as.numeric(c(V1.var,covar[2],covar[2],k12.var)),2,2,byrow=T)
  Cl2_deriv<-as.matrix(attr(eval(stats::deriv(~V1*k12,c("V1","k12"))),
                            "gradient"))
  Cl2.sd<-sqrt(Cl2_deriv %*% sigma22 %*% t(Cl2_deriv))

  sigma23<-matrix(as.numeric(c(V1.var,covar[3],covar[3],k13.var)),2,2,byrow=T)
  Cl3_deriv<-as.matrix(attr(eval(stats::deriv(~V1*k13,c("V1","k13"))),
                            "gradient"))
  Cl3.sd<-sqrt(Cl3_deriv %*% sigma23 %*% t(Cl3_deriv))

  root1<-(-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
        ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3))*(2*exp(log(sqrt((-((
        (k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/
        3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))

  root2<-(-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
        ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*exp(log(sqrt((-
        (((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+
        k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))

  root3<-(-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
        ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+4*pi/3)*(2*exp(log(sqrt((-
        (((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+
        k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))

  root<-c(root1,root2,root3)
  if(is.na(root[1])){
    l1<-l2<-l3<-NA
  }else{
    l1<-max(root);  l2<-stats::median(root);   l3<-min(root)
  }
  t_alpha<-log(2)/l1;  t_beta<-log(2)/l2;    t_gamma<-log(2)/l3

  t_alpha_deriv<-as.matrix(attr(eval(
        stats::deriv(~log(2)/(-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
          ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
          (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3))*(2*exp(log(sqrt((-
          (((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)),
          c("V1","k10","k12","k21","k13","k31"))),"gradient"))

  t_beta_deriv<-as.matrix(attr(eval(
        stats::deriv(~log(2)/(-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
          ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
          (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*
          (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)),
          c("V1","k10","k12","k21","k13","k31"))),"gradient"))

  t_gamma_deriv<-as.matrix(attr(eval(
        stats::deriv(~log(2)/(-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
          ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
          (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+4*pi/3)*
          (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)),
          c("V1","k10","k12","k21","k13","k31"))),"gradient"))

  sigma6<-matrix(as.numeric(c(V1.var,covar[1],covar[2],covar[3],covar[4],
                 covar[5], covar[1],k10.var,covar[6],covar[7],covar[8],covar[9],
                 covar[2],covar[6],k12.var,covar[9],covar[11],covar[12],
                 covar[3],covar[7],covar[10],k13.var,covar[13],covar[14],
                 covar[4],covar[8],covar[11],covar[13],k21.var,covar[15],
                 covar[5],covar[9],covar[12],covar[14],covar[15],k31.var)),
                 6,6,byrow=T)

  half_sig1<-sqrt(t_alpha_deriv %*% sigma6 %*% t(t_alpha_deriv))
  half_sig2<-sqrt(t_beta_deriv %*% sigma6 %*% t(t_beta_deriv))
  half_sig3<-sqrt(t_gamma_deriv %*% sigma6 %*% t(t_gamma_deriv))

  half_sig<-c(half_sig1,half_sig2,half_sig3)
  rt_ord<-order(root)
  t_alpha.sd<-half_sig[which(rt_ord==3)]
  t_beta.sd<-half_sig[which(rt_ord==2)]
  t_gamma.sd<-half_sig[which(rt_ord==1)]

  true_A<-(k21-l1)*(k31-l1)/(l1-l2)/(l1-l3)/V1
  true_B<-(k21-l2)*(k31-l2)/(l2-l1)/(l2-l3)/V1
  true_C<-(k21-l3)*(k31-l3)/(l3-l2)/(l3-l1)/V1

  f.true_A<-quote(quote((k21-((-(cos((acos((-(((2*(
        (k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+
        k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-((
        (k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/
        3))^3))/27)))/3))*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+
        k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+
        k31)/3))))*(k31-((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
        ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3))*(2*exp(log(sqrt((-((
        (k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/
        3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))/(((-(cos((acos((-(((2*(
        (k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*
        k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+
        k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/
        27)))/3))*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*
        k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/
        3)))-((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+
        k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*
        (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+
        k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))/
        (((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/
        2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+
        k21+k31)^2)/3))^3))/27)))/3))*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*
        k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-
        (k10+k12+k13+k21+k31)/3)))-((-(cos((acos((-(((2*((k10+k12+k13+k21+
        k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+
        k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*
        k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+4*pi/3)*
        (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/
        3))))/V1))
   ff.true_A<-stats::as.formula(paste("~",as.character(f.true_A[2],"")))
   true_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_A,
                        c("V1","k10","k12","k21","k13","k31"))),"gradient"))

   f.true_B<-quote(quote((k21-((-(cos((acos((-(((2*((k10+k12+
        k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*
        (k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+
        2*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))*
        (k31-((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+
        k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*
        exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/
        3))))/(((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
        ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*exp(log(sqrt((-(
        ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/
        3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)))-((-(cos((acos((-(((2*(
        (k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*
        k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/
        27)))/3))*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))/
        (((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/
        (sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+
        k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-
        (k10+k12+k13+k21+k31)/3)))-((-(cos((acos((-(((2*((k10+k12+k13+k21+
        k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*
        (k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+
        4*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)))
        )/V1))
   ff.true_B<-stats::as.formula(paste("~",as.character(f.true_B[2],"")))
   true_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_B,
                          c("V1","k10","k12","k21","k13","k31"))),"gradient"))

   f.true_C<-quote(quote((k21-((-(cos((acos((-(((2*((k10+k12+
        k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+
        k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*
        k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+
        4*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))*
        (k31-((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/
        2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+
        k21+k31)^2)/3))^3))/27)))/3)+4*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/
        3))-(k10+k12+k13+k21+k31)/3))))/(((-(cos((acos((-(((2*((k10+k12+k13+
        k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*
        (k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/
        3)+4*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*
        k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+
        k31)/3)))-((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
        ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*exp(log(sqrt((-(
        ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/
        3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))/(((-(cos((acos((-(((2*
        ((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*
        k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/
        3)+4*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*
        k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/
        3)))-((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/
        2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+
        k21+k31)^2)/3))^3))/27)))/3))*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*
        k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+
        k12+k13+k21+k31)/3))))/V1))
   ff.true_C<-stats::as.formula(paste("~",as.character(f.true_C[2],"")))
   true_C_deriv<-as.matrix(attr(eval(stats::deriv(ff.true_C,
                      c("V1","k10","k12","k21","k13","k31"))),"gradient"))

   true_sig1<-sqrt(true_A_deriv %*% sigma6 %*% t(true_A_deriv))
   true_sig2<-sqrt(true_B_deriv %*% sigma6 %*% t(true_B_deriv))
   true_sig3<-sqrt(true_C_deriv %*% sigma6 %*% t(true_C_deriv))
   true_sig<-c(true_sig1,true_sig2,true_sig3)

   true_A.sd<-true_sig[which(rt_ord==3)];
   true_B.sd<-true_sig[which(rt_ord==2)];
   true_C.sd<-true_sig[which(rt_ord==1)]

   frac_A<-(k21-l1)*(k31-l1)/(l1-l2)/(l1-l3)
   frac_B<-(k21-l2)*(k31-l2)/(l2-l1)/(l2-l3)
   frac_C<-(k21-l3)*(k31-l3)/(l3-l2)/(l3-l1)

   f.frac_A<-quote(quote((k21-((-(cos((acos((-(((2*((k10+k12+
        k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*
        (k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/
        3))*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))*
        (k31-((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/
        2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+
        k21+k31)^2)/3))^3))/27)))/3))*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*
        k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-
        (k10+k12+k13+k21+k31)/3))))/(((-(cos((acos((-(((2*((k10+k12+k13+k21+
        k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+
        k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*
        k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3))*(2*
        exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+
        k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)))-
        ((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/
        2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+
        k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/
        3))-(k10+k12+k13+k21+k31)/3))))/(((-(cos((acos((-(((2*((k10+k12+k13+
        k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+
        k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+
        k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3))*
        (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)))-
        ((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/
        (sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+
        k21+k31)^2)/3))^3))/27)))/3)+4*pi/3)*(2*exp(log(sqrt((-(((k10*k31+
        k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/
        27))/3))-(k10+k12+k13+k21+k31)/3))))))

   ff.frac_A<-stats::as.formula(paste("~",as.character(f.frac_A[2],"")))
   frac_A_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_A,
                          c("V1","k10","k12","k21","k13","k31"))),"gradient"))
   f.frac_B<-quote(quote((k21-((-(cos((acos((-(((2*((k10+
        k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*
        (k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+
        2*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))*
        (k31-((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+
        k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*
        exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))/
        (((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/
        (sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+
        k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*exp(log(sqrt((-(((k10*k31+
        k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/
        27))/3))-(k10+k12+k13+k21+k31)/3)))-((-(cos((acos((-(((2*((k10+k12+
        k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*
        (k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3))*
        (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3))))/
        (((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/
        (sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+
        k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-
        (k10+k12+k13+k21+k31)/3)))-((-(cos((acos((-(((2*((k10+k12+k13+k21+
        k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+
        k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+
        k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+4*pi/3)*
        (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-
        (k10+k12+k13+k21+k31)/3))))))
    ff.frac_B<-stats::as.formula(paste("~",as.character(f.frac_B[2],"")))
    frac_B_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_B,
                          c("V1","k10","k12","k21","k13","k31"))),"gradient"))

    f.frac_C<-quote(quote((k21-((-(cos((acos((-(((2*((k10+k12+
        k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*
        (k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*
        k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/
        3)+4*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+
        k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+
        k21+k31)/3))))*(k31-((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/
        27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/
        3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+
        k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+4*pi/3)*
        (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/
        3))))/(((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
        ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
        (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+4*pi/3)*(2*exp(log(sqrt(
        (-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+
        k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)))-
        ((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/
        (sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+
        k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*(2*exp(log(sqrt((-(((k10*k31+
        k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/
        27))/3))-(k10+k12+k13+k21+k31)/3))))/(((-(cos((acos((-(((2*((k10+k12+
        k13+k21+k31)^3))/27)-((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*
        (k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+
        4*pi/3)*(2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
        (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)))-
        ((-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+(k10*k21*k31))/2)/
        (sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+
        k21+k31)^2)/3))^3))/27)))/3))*(2*exp(log(sqrt((-(((k10*k31+k21*k31+
        k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-
        (k10+k12+k13+k21+k31)/3))))))
   ff.frac_C<-stats::as.formula(paste("~",as.character(f.frac_C[2],"")))
   frac_C_deriv<-as.matrix(attr(eval(stats::deriv(ff.frac_C,
                       c("V1","k10","k12","k21","k13","k31"))),"gradient"))

   frac_sig1<-sqrt(frac_A_deriv %*% sigma6 %*% t(frac_A_deriv))
   frac_sig2<-sqrt(frac_B_deriv %*% sigma6 %*% t(frac_B_deriv))
   frac_sig3<-sqrt(frac_C_deriv %*% sigma6 %*% t(frac_C_deriv))

   frac_sig<-c(frac_sig1,frac_sig2,frac_sig3)

   frac_A.sd<-frac_sig[which(rt_ord==3)]
   frac_B.sd<-frac_sig[which(rt_ord==2)]
   frac_C.sd<-frac_sig[which(rt_ord==1)]

   alpha<-l1;   beta<-l2;  gamma<-l3

   alpha_deriv<-as.matrix(attr(eval(
        stats::deriv(~(-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
          ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
          (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3))*(2*exp(log(sqrt((-((
          (k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-(((k10+k12+k13+k21+
          k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)),
          c("V1","k10","k12","k21","k13","k31"))),"gradient"))

   beta_deriv<-as.matrix(attr(eval(
        stats::deriv(~(-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
          ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
          (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+2*pi/3)*
          (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)),
          c("V1","k10","k12","k21","k13","k31"))),"gradient"))

   gamma_deriv<-as.matrix(attr(eval(
        stats::deriv(~(-(cos((acos((-(((2*((k10+k12+k13+k21+k31)^3))/27)-
          ((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)*(k10+k12+k13+k21+k31)/3)+
          (k10*k21*k31))/2)/(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27)))/3)+4*pi/3)*
          (2*exp(log(sqrt((-(((k10*k31+k21*k31+k21*k13+k10*k21+k31*k12)-
          (((k10+k12+k13+k21+k31)^2)/3))^3))/27))/3))-(k10+k12+k13+k21+k31)/3)),
          c("V1","k10","k12","k21","k13","k31"))),"gradient"))

   exp_sig1<-sqrt(alpha_deriv %*% sigma6 %*% t(alpha_deriv))
   exp_sig2<-sqrt(beta_deriv %*% sigma6 %*% t(beta_deriv))
   exp_sig3<-sqrt(gamma_deriv %*% sigma6 %*% t(gamma_deriv))

   exp_sig<-c(exp_sig1,exp_sig2,exp_sig3)
   alpha.sd<-exp_sig[which(rt_ord==3)]
   beta.sd<-exp_sig[which(rt_ord==2)]
   gamma.sd<-exp_sig[which(rt_ord==1)]
   if(is.na(V1[1])){
     param = rep(NA,24)
     sd = rep(NA,24)
   } else{
     param = c(V1,k10,k12,k13,k21,k31,V2,V3,Cl1,Cl2,Cl3,t_alpha,t_beta,
               t_gamma,Vdss,true_A,true_B,true_C,frac_A,frac_B,frac_C,
               alpha,beta,gamma)
     sd = c(V1.sd,k10.sd,k12.sd,k13.sd,k21.sd,k31.sd,V2.sd,V3.sd,Cl1.sd,
            Cl2.sd,Cl3.sd,t_alpha.sd,t_beta.sd,t_gamma.sd,Vdss.sd,
            true_A.sd,true_B.sd,true_C.sd,frac_A.sd,frac_B.sd,frac_C.sd,
            alpha.sd,beta.sd,gamma.sd)
   }
   result = data.frame(Parameter=c("V1","k10","k12","k13","k21","k31",
                       "V2","V3","Cl1","Cl2","Cl3","t_alpha","t_beta",
                       "t_gamma","Vdss","True_A","True_B","True_C","Frac_A",
                       "Frac_B","Frac_C","alpha","beta","gamma"),
                       Estimate=param, Std.err=sd)
   row.names(result) <- c("V1","k10","k12","k13","k21","k31","V2","V3",
                          "Cl1","Cl2","Cl3","t_alpha","t_beta","t_gamma",
                          "Vdss","True_A","True_B","True_C","Frac_A","Frac_B",
                          "Frac_C","alpha","beta","gamma")
  result<-result[c("Vdss","V1","V2","V3", "Cl1","Cl2","Cl3",
                   "k10","k12","k21","k13","k31","alpha","beta","gamma",
                   "t_alpha","t_beta","t_gamma","True_A","True_B","True_C",
                   "Frac_A","Frac_B","Frac_C"),]
  return(result)

}
