#' Convert pharmacokinetic parameters for three compartment model
#'
#' Calculate pharmacokinetic parameters with parameters (A, B, C,
#' alpha, beta, and gamma) in two compartment model
#' "Aexp(-alpha)+Bexp(-beta)+Cexp(-gamma)"
#'
#' @usage ThreeComp_Coefficient_Exponent(A,B,C,alpha,beta,gamma,A.sd=NA,
#'    B.sd=NA,C.sd=NA,alpha.sd=NA,beta.sd=NA,gamma.sd=NA,
#'    covar=c(AB=NA,AC=NA,Aalpha=NA,Abeta=NA,Agamma=NA,BC=NA,Balpha=NA,
#'    Bbeta=NA,Bgamma=NA,Calpha=NA,Cbeta=NA,Cgamma=NA,alphabeta=NA,
#'    alphagamma=NA,betagamma=NA),...)
#' @param A parameter in one compartment model "Aexp(-alpha)"
#' @param B parameter in two compartment model "Aexp(-alpha)+Bexp(-beta)"
#' @param C parameter in three compartment model
#'          "Aexp(-alpha)+Bexp(-beta)+Cexp(-gamma)"
#' @param alpha parameter in one compartment model "Aexp(-alpha)"
#' @param beta parameter in two compartment model "Aexp(-alpha)+Bexp(-beta)"
#' @param gamma parameter in three compartment model
#'              "Aexp(-alpha)+Bexp(-beta)+Cexp(-gamma)"
#' @param A.sd standard error of A
#' @param B.sd standard error of B
#' @param C.sd standard error of C
#' @param alpha.sd standard error of alpha
#' @param beta.sd standard error of beta
#' @param gamma.sd standard error of gamma
#' @param covar covariances among parameters
#' @param ... arguments to be passed to methods
#' @references \url{http://www.nonmemcourse.com/convert.xls}
#' @export
#' @examples
#' ThreeComp_Coefficient_Exponent(A=12.2,B=3.76,C=1.44,
#' alpha=0.870,beta=0.12,gamma=0.013, A.sd=0.2,B.sd=0.005,C.sd=0.0005,
#' alpha.sd=0.009,beta.sd=0.006,gamma.sd=0.00005)


ThreeComp_Coefficient_Exponent<-function(A,B,C,alpha,beta,gamma,A.sd=NA,
    B.sd=NA,C.sd=NA,alpha.sd=NA,beta.sd=NA,gamma.sd=NA,
    covar=c(AB=NA,AC=NA,Aalpha=NA,Abeta=NA,Agamma=NA,BC=NA,Balpha=NA,
    Bbeta=NA,Bgamma=NA,Calpha=NA,Cbeta=NA,Cgamma=NA,alphabeta=NA,
    alphagamma=NA,betagamma=NA),...){
  if(is.na(covar[1])) covar<-rep(0,15)
  A.var = (A.sd)^2;           B.var = (B.sd)^2;        C.var = (C.sd)^2
  alpha.var = (alpha.sd)^2;   beta.var = (beta.sd)^2;  gamma.var = (gamma.sd)^2

  V1<-1/(A+B+C)

  f.V2<-quote(quote((1/(A+B+C))*(((beta*gamma+alpha*beta+alpha*gamma)-
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+sqrt((-(alpha*(C+B)+
        gamma*(A+B)+beta*(A+C))/(A+B+C))^2-4*((alpha*beta*C+alpha*gamma*B+
        beta*gamma*A)/(A+B+C))))/2)*(alpha+beta+gamma)-(alpha*beta*gamma/
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+sqrt((-(alpha*(C+B)+
        gamma*(A+B)+beta*(A+C))/(A+B+C))^2-4*((alpha*beta*C+alpha*gamma*B+
        beta*gamma*A)/(A+B+C))))/2))+((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/
        (A+B+C))+sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)^2)/
        (-sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C)))))/
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+sqrt((-
        (alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-4*
        ((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))

  V2<-eval(eval(f.V2))
  ff.V2<-stats::as.formula(paste("~",as.character(f.V2[2],"")))
  f.V3<-quote(quote((V1)*(alpha+beta+gamma-((alpha*beta*gamma/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))+
        (((beta*gamma+alpha*beta+alpha*gamma)-
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        (alpha+beta+gamma)-(alpha*beta*gamma/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))*
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))/
        (((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)-
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))

  V3<-eval(eval(f.V3))
  ff.V3<-stats::as.formula(paste("~",as.character(f.V3[2],"")))

  V1_deriv<-as.matrix(attr(eval(stats::deriv(~1/(A+B+C),c("A","B","C"))),
                           "gradient"))
  V2_deriv<-as.matrix(attr(eval(stats::deriv(ff.V2,
    c("A","B","C","alpha","beta","gamma"))),"gradient"))
  V3_deriv<-as.matrix(attr(eval(stats::deriv(ff.V3,
    c("A","B","C","alpha","beta","gamma"))),"gradient"))

  f.Vdss<-quote(quote((1/(A+B+C))+((1/(A+B+C))*(((beta*gamma+alpha*beta+
        alpha*gamma)-((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+
        sqrt((-(alpha*(C+B)+
        gamma*(A+B)+beta*(A+C))/(A+B+C))^2-4*((alpha*beta*C+alpha*gamma*B+
        beta*gamma*A)/(A+B+C))))/2)*(alpha+beta+gamma)-(alpha*beta*gamma/
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+sqrt((-(alpha*(C+B)+
        gamma*(A+B)+beta*(A+C))/(A+B+C))^2-4*((alpha*beta*C+alpha*gamma*B+
        beta*gamma*A)/(A+B+C))))/2))+((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/
        (A+B+C))+sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)^2)/
        (-sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C)))))/
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+sqrt((-
        (alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-4*
        ((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))+
        ((V1)*(alpha+beta+gamma-((alpha*beta*gamma/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))+
        (((beta*gamma+alpha*beta+alpha*gamma)-
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        (alpha+beta+gamma)-(alpha*beta*gamma/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))*
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))/
        (((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)-
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))))

  Vdss<-eval(eval(f.Vdss))
  ff.Vdss<-stats::as.formula(paste("~",as.character(f.Vdss[2],"")))
  Vdss_deriv<-as.matrix(attr(eval(stats::deriv(ff.Vdss,
                         c("A","B","C","alpha","beta","gamma"))),"gradient"))

  sigma3<-matrix(as.numeric(c(A.var,covar[1],covar[2],
                              covar[1],B.var,covar[6],
                              covar[2],covar[6],C.var)),3,3,byrow = T)

  sigma6<-matrix(as.numeric(c(A.var,covar[1],covar[2],covar[3],covar[4],
                 covar[5],covar[1],B.var,covar[6],covar[7],covar[8],covar[9],
                 covar[2],covar[6],C.var,covar[9],covar[11],covar[12],
                 covar[3],covar[7],covar[10],alpha.var,covar[13],covar[14],
                 covar[4],covar[8],covar[11],covar[13],beta.var,covar[15],
                 covar[5],covar[9],covar[12],covar[14],covar[15],gamma.var)),
                 6,6,byrow=T)

  V1.sd<-sqrt(V1_deriv %*% sigma3 %*% t(V1_deriv))
  V2.sd<-sqrt(V2_deriv %*% sigma6 %*% t(V2_deriv))
  V3.sd<-sqrt(V3_deriv %*% sigma6 %*% t(V3_deriv))
  Vdss.sd<-sqrt(Vdss_deriv %*% sigma6 %*% t(Vdss_deriv))

  f.Cl1<-quote(quote((1/(A+B+C))*((alpha*beta*gamma)/
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))-
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))))

  Cl1<-eval(eval(f.Cl1))
  ff.Cl1<-stats::as.formula(paste("~",as.character(f.Cl1[2],"")))
  f.Cl2<-quote(quote((1/(A+B+C))*(((beta*gamma+alpha*beta+alpha*gamma)-
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        (alpha+beta+gamma)-
        (alpha*beta*gamma/((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))+
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)^2)/
        (-sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C)))))))

  Cl2<-eval(eval(f.Cl2))
  ff.Cl2<-stats::as.formula(paste("~",as.character(f.Cl2[2],"")))
  f.Cl3<-quote(quote((1/(A+B+C))*(alpha+beta+gamma-((alpha*beta*gamma/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))+
        (((beta*gamma+alpha*beta+alpha*gamma)-
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        (alpha+beta+gamma)-(alpha*beta*gamma/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))*
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))/
        (((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)-
          ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))))

  Cl3<-eval(eval(f.Cl3))
  ff.Cl3<-stats::as.formula(paste("~",as.character(f.Cl3[2],"")))

  Cl1_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl1,
                 c("A","B","C","alpha","beta","gamma"))),"gradient"))
  Cl2_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl2,
                 c("A","B","C","alpha","beta","gamma"))),"gradient"))
  Cl3_deriv<-as.matrix(attr(eval(stats::deriv(ff.Cl3,
                 c("A","B","C","alpha","beta","gamma"))),"gradient"))

  Cl1.sd<-sqrt(Cl1_deriv %*% sigma6 %*% t(Cl1_deriv))
  Cl2.sd<-sqrt(Cl2_deriv %*% sigma6 %*% t(Cl2_deriv))
  Cl3.sd<-sqrt(Cl3_deriv %*% sigma6 %*% t(Cl3_deriv))

  f.k10<-quote(quote((alpha*beta*gamma)/((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/
      (A+B+C))+sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))-
         sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-4*
                ((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))

  k10<-eval(eval(f.k10))
  ff.k10<-stats::as.formula(paste("~",as.character(f.k10[2],"")))

  f.k12<-quote(quote(((beta*gamma+alpha*beta+alpha*gamma)-
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        (alpha+beta+gamma)-(alpha*beta*gamma/
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))+
        ((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)^2)/
        (-sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))))

  k12<-eval(eval(f.k12))
  ff.k12<-stats::as.formula(paste("~",as.character(f.k12[2],"")))

  f.k13<-quote(quote(alpha+beta+gamma-((alpha*beta*gamma/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))+
        (((beta*gamma+alpha*beta+alpha*gamma)-
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        (alpha+beta+gamma)-(alpha*beta*gamma/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)/
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))*
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)*
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))/
        (((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)-
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)))+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))+
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2)+
        ((-((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))-
        sqrt(((-(alpha*C+alpha*B+gamma*A+gamma*B+beta*A+beta*C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))))

  k13<-eval(eval(f.k13))
  ff.k13<-stats::as.formula(paste("~",as.character(f.k13[2],"")))

  f.k21<-quote(quote((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))+
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))

  k21<-eval(eval(f.k21))
  ff.k21<-stats::as.formula(paste("~",as.character(f.k21[2],"")))

  f.k31<-quote(quote((((alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))-
        sqrt((-(alpha*(C+B)+gamma*(A+B)+beta*(A+C))/(A+B+C))^2-
        4*((alpha*beta*C+alpha*gamma*B+beta*gamma*A)/(A+B+C))))/2))

  k31<-eval(eval(f.k31))
  ff.k31<-stats::as.formula(paste("~",as.character(f.k31[2],"")))

  k10_deriv<-as.matrix(attr(eval(stats::deriv(ff.k10,
        c("A","B","C","alpha","beta","gamma"))),"gradient"))
  k12_deriv<-as.matrix(attr(eval(stats::deriv(ff.k12,
        c("A","B","C","alpha","beta","gamma"))),"gradient"))
  k13_deriv<-as.matrix(attr(eval(stats::deriv(ff.k13,
        c("A","B","C","alpha","beta","gamma"))),"gradient"))
  k21_deriv<-as.matrix(attr(eval(stats::deriv(ff.k21,
        c("A","B","C","alpha","beta","gamma"))),"gradient"))
  k31_deriv<-as.matrix(attr(eval(stats::deriv(ff.k31,
        c("A","B","C","alpha","beta","gamma"))),"gradient"))

  k10.sd<-sqrt(k10_deriv %*% sigma6 %*% t(k10_deriv))
  k12.sd<-sqrt(k12_deriv %*% sigma6 %*% t(k12_deriv))
  k13.sd<-sqrt(k13_deriv %*% sigma6 %*% t(k13_deriv))
  k21.sd<-sqrt(k21_deriv %*% sigma6 %*% t(k21_deriv))
  k31.sd<-sqrt(k31_deriv %*% sigma6 %*% t(k31_deriv))

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

  Frac_A<-V1*A
  Frac_B<-V1*B
  Frac_C<-V1*C

  Frac_A_deriv<-as.matrix(attr(eval(stats::deriv(~A/(A+B+C),
                                          c("A","B","C"))),"gradient"))
  Frac_B_deriv<-as.matrix(attr(eval(stats::deriv(~B/(A+B+C),
                                          c("A","B","C"))),"gradient"))
  Frac_C_deriv<-as.matrix(attr(eval(stats::deriv(~C/(A+B+C),
                                          c("A","B","C"))),"gradient"))

  Frac_A.sd<-sqrt(Frac_A_deriv %*% sigma3 %*% t(Frac_A_deriv))
  Frac_B.sd<-sqrt(Frac_B_deriv %*% sigma3 %*% t(Frac_B_deriv))
  Frac_C.sd<-sqrt(Frac_C_deriv %*% sigma3 %*% t(Frac_C_deriv))
  if(is.na(A[1])){
    param = rep(NA,24)
    sd = rep(NA,24)
  } else{
    param = c(A,B,C,alpha,beta,gamma,V1,V2,V3,Cl1,Cl2,Cl3,k10,k12,k13,k21,k31,
              t_alpha,t_beta,t_gamma,Vdss,Frac_A,Frac_B,Frac_C)
    sd = c(A.sd,B.sd,C.sd,alpha.sd,beta.sd,gamma.sd,V1.sd,V2.sd,V3.sd,
           Cl1.sd,Cl2.sd,Cl3.sd,k10.sd,k12.sd,k13.sd,k21.sd,k31.sd,t_alpha.sd,
           t_beta.sd,t_gamma.sd,Vdss.sd,Frac_A.sd,Frac_B.sd,Frac_C.sd)
  }
  result = data.frame(Parameter=c("True_A","True_B","True_C","alpha","beta",
                      "gamma","V1","V2","V3","Cl1","Cl2","Cl3","k10","k12",
                      "k13","k21","k31","t_alpha","t_beta","t_gamma","Vdss",
                      "Frac_A","Frac_B","Frac_C"),
                      Estimate=param, Std.err=sd)
  row.names(result) <- c("True_A","True_B","True_C","alpha","beta","gamma",
    "V1","V2","V3","Cl1","Cl2","Cl3","k10","k12","k13","k21","k31",
    "t_alpha","t_beta","t_gamma","Vdss","Frac_A","Frac_B","Frac_C")
  result<-result[c("Vdss","V1","V2","V3","Cl1","Cl2","Cl3",
                   "k10","k12","k21","k13","k31","alpha","beta","gamma",
                   "t_alpha","t_beta","t_gamma","True_A","True_B","True_C",
                   "Frac_A","Frac_B","Frac_C"),]
  return(result)
}
