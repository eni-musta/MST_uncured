library(smcure) # for the data
library(bpcp)   # for the test for difference of survival probabilities 

#####################################################################
# comparison of mean survival times for the uncured (nonparametric) #
#####################################################################


# Input: Data1, Data2 - two datasets (corresponding to the two groups to be compared), 
#                       first column is survival time, second column is censoring indicator
#                       cure_comp - logical: if TRUE computes also parameters related to the cure component
# Output: m - difference in mean survival time for the uncured
#       var - asymptotic variance of a_n(hat_m-m)
#       If cure_comp==TRUE returns also:
#       p1  - cure fraction in sample 1
#       p2  - cure fraction in sample 2
#       CI1 - 95% CI for p1
#       CI2 - 95% CI for p2
#      pval - p-val for testing H0: p1=p2 against H1: p1!=p2

estimation=function(Data1,Data2,cure_comp){
  
  colnames(Data1)=c("Y","status")
  colnames(Data2)=c("Y","status")
  
  n1=dim(Data1)[1]
  n2=dim(Data2)[1]
  
  km_fit1 <- survfit(Surv(Y, status) ~ 1, data=Data1)
  km_fit2 <- survfit(Surv(Y, status) ~ 1, data=Data2)
  
  #cure fractions
  p1=min(km_fit1$surv)
  p2=min(km_fit2$surv)
  
  
  #observed times (excluding ties)
  times1=km_fit1$time
  times2=km_fit2$time
  
  #survival functions
  s1=km_fit1$surv
  surv1=stepfun(c(0,times1),c(0,1,s1),f=0,right=FALSE)
  s2=km_fit2$surv
  surv2=stepfun(c(0,times2),c(0,1,s2),f=0,right=FALSE)
   
  # integrands for computing the mean survival time of the uncured (piecewise constant)
  f1=(s1-p1)/(1-p1)
  f2=(s2-p2)/(1-p2)
  
  # mean survival time of the uncured 
  E1=sum(diff(c(0,times1))*(c(1,f1[-length(f1)])))
  E2=sum(diff(c(0,times2))*(c(1,f2[-length(f2)])))
  
  # difference of mean survival times of the uncured 
  m=E1-E2
  
  # estimate of tau_{0,i}
  tau_1=max(Data1$Y[which(Data1$status==1)])
  tau_2=max(Data2$Y[which(Data2$status==1)])
  
  # v function that appears in the asymptotic variance expression
  v1_val=n1*km_fit1$std.err^2
  v1=stepfun(c(0,times1),c(0,0,n1*km_fit1$std.err^2),f=0,right=FALSE)
  v2=stepfun(c(0,times2),c(0,0,n2*km_fit2$std.err^2),f=0,right=FALSE)
  v2_val=n2*km_fit2$std.err^2
  
  # computation of assymptotic variance
  
  # for the first sample
  integrand1=function(u){
    val=rep(NA,length(u))
    for(i in 1:length(u)){
      t_u=c(0,times1)[which(c(0,times1)<=u[i])]
      t_u=max(t_u)
      ind=which(c(0,times1)==t_u)
      int11=0
      if(t_u>0){
        int11=sum(diff(c(0,times1)[1:ind])*(c(0,s1[1:(ind-2)]*v1_val[1:(ind-2)])))
      }
      if(u[i]>t_u){
        int11=int11+(u[i]-t_u)*surv1(u[i])*v1(u[i])
      }
      ind_tau=which(c(0,times1)==tau_1)
      int12=0
      if(u[i]<tau_1){
        if(c(0,times1)[ind+1]<tau_1){
          int12=sum(diff(c(0,times1)[(ind+1):ind_tau])*(c(0,s1)[(ind+1):(ind_tau-1)]))
        }
        int12=int12+(c(0,times1)[ind+1]-u[i])*surv1(u[i])
      }
      val[i]=surv1(u[i])*(int11+v1(u[i])*int12)
    }
    return(val)
  }
  
  integral1=integrate(integrand1,0,tau_1,subdivisions=2000)$value
  ind_tau=which(c(0,times1)==tau_1)
  integral2=sum(diff(c(0,times1)[1:ind_tau])*(c(0,s1[1:(ind_tau-2)]*v1_val[1:(ind_tau-2)])))
  
  sigma_squared1=integral1/((1-p1)^2)+p1^2*((1-p1)^{-2})*((E1-tau_1)^2)*v1(tau_1)+2*p1*((1-p1)^{-2})*(E1-tau_1)*integral2
   
  # for the second sample
  integrand2=function(u){
    val=rep(NA,length(u))
    for(i in 1:length(u)){
      t_u=c(0,times2)[which(c(0,times2)<=u[i])]
      t_u=max(t_u)
      ind=which(c(0,times2)==t_u)
      int21=0
      if(t_u>0){
        int21=sum(diff(c(0,times2)[1:ind])*(c(0,s2[1:(ind-2)]*v2_val[1:(ind-2)])))
      }
      if(u[i]>t_u){
        int21=int21+(u[i]-t_u)*surv2(u[i])*v2(u[i])
      }
      ind_tau=which(c(0,times2)==tau_2)
      int22=0
      if(u[i]<tau_2){
        if(c(0,times2)[ind+1]<tau_2){
          int22=sum(diff(c(0,times2)[(ind+1):ind_tau])*(c(0,s2)[(ind+1):(ind_tau-1)]))
        }
        int22=int22+(c(0,times2)[ind+1]-u[i])*surv2(u[i])
      }
      val[i]=surv2(u[i])*(int21+v2(u[i])*int22)
    }
    return(val)
  }
  
  integral1=integrate(integrand2,0,tau_2,subdivisions=2000)$value
  ind_tau=which(c(0,times2)==tau_2)
  integral2=sum(diff(c(0,times2)[1:ind_tau])*(c(0,s2[1:(ind_tau-2)]*v2_val[1:(ind_tau-2)])))
  
  sigma_squared2=integral1/((1-p2)^2)+p2^2*((1-p2)^{-2})*((E2-tau_2)^2)*v2(tau_2)+2*p2*((1-p2)^{-2})*(E2-tau_2)*integral2
  
  # final variance
  sigma_squared=(n2/(n1+n2))*sigma_squared1+(n1/(n1+n2))*sigma_squared2
  
  if(cure_comp==FALSE){
    return(list(m=m,var=sigma_squared))
  }
  if(cure_comp==TRUE){
    # analysis of the cure fraction
    
    CI1=c(min(km_fit1$lower),min(km_fit1$upper))
    CI2=c(min(km_fit2$lower),min(km_fit2$upper))
    
    d=rbind(Data1,Data2)
    d=cbind(d,c(rep(0,n1),rep(1,n2)))
    colnames(d)=c("Y","status","TRT")
    pval=fixtdiff(d$Y,d$status,d$TRT,testtime=max(d$Y[which(d$status==1)]),trans='cloglog', varpooled=FALSE)$p2
    
    return(list(m=m,var=sigma_squared,p1=p1,p2=p2,CI1=CI1,CI2=CI2,pval=pval))
  }
}

#######################
# Example smcure data #
#######################


data(e1684)  
d=na.omit(e1684) 

Data1=d[which(d$TRT==0),c(2,3)]
Data2=d[which(d$TRT==1),c(2,3)]

# code below is standard, not depending on dataset
# preliminary analysis of the data

colnames(Data1)=c("Y","status")
colnames(Data2)=c("Y","status")

n1=dim(Data1)[1]
n2=dim(Data2)[1]

cens_rate_1=length(which(Data1$status==0))/n1 
cens_rate_2=length(which(Data2$status==0))/n2

Y_r1=max(Data1$Y[which(Data1$status==1)]) # maximum observed event time
Y_r2=max(Data2$Y[which(Data2$status==1)])

plateau_1=length(which(Data1$Y>Y_r1))/n1  # percentage of observations in the plateau
plateau_2=length(which(Data2$Y>Y_r2))/n2

##############
# Estimation #
##############

est_res=estimation(Data1,Data2,cure_comp = TRUE)

#cure fraction and 95% CI for sample 1
p1=est_res$p1
l1=est_res$CI1[1]
u1=est_res$CI1[2]

#cure fraction and 95% CI for sample 2
p2=est_res$p2
l2=est_res$CI2[1]
u2=est_res$CI2[2]

#p-val for testing H0: p1=p2 against H1: p1!=p2 (comparison of survival probabilities at the final point using test of Klein (2007)
p_val_cure=est_res$pval

# comparison of mean survival times for the uncured 
m=est_res$m
sds=sqrt(est_res$var)

k=sqrt(n1*n2/(n1+n2))

# asymptotic pvalue for testing H0: m=0 aginst H1: m!=0
as_pval=2*(1-pnorm(abs(m*k/sds)))

#asymptotic 95% CI for m
m+qnorm(0.975)*sds/k
m-qnorm(0.975)*sds/k

########################
# Permutation approach #
########################

K=500
n=n1+n2
pooled_data=rbind(Data1,Data2)

m_p=rep(0,K)
sd_p=rep(0,K)


for(j in 1:K){
  
  id<-sample(1:n,n,replace=FALSE)
  new_Data1=pooled_data[id[1:n1],]
  new_Data2=pooled_data[id[(n1+1):n],]
  
  perm_est=estimation(new_Data1,new_Data2,cure_comp = FALSE)
  m_p[j]=perm_est$m
  sd_p[j]=sqrt(perm_est$var)
  
}

q=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
quantiles=c(quantile(m_p/sd_p,q))
quantiles=c(quantile(na.omit(m_p/sd_p),q))


# pvalue with permutation approach for testing H0: m=0 aginst H1: m!=0
perm_pval=(length(which(m_p/sd_p<(-abs(m/sds))))+length(which(m_p/sd_p>abs(m/sds))))/K

# 95% CI for m with permutation approach
m-quantiles[1]*sds
m-quantiles[7]*sds
