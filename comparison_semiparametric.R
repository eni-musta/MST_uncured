library(smcure) 

###################################################################
# conditional comparison of mean survival times for the uncured   #
###################################################################

# load required functions 

# the four following functions are from the smcure package for the logistic-Cox cure model (with small modifications)

# survival function estimation for the Cox part
surv_cox = function(Time,Status,Z,beta,w){    
  death_point = sort(unique(subset(Time, Status==1)))
  coxexp = exp((beta)%*%t(Z))  
  lambda = numeric()
  event = numeric()
  for(i in 1: length(death_point)){
    event[i] = sum(Status*as.numeric(Time==death_point[i]))
    temp = sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
    lambda[i] = event[i]/temp
  }
  CumHazard = numeric()
  for(i in 1:length(Time)){
    CumHazard[i] = sum(as.numeric(Time[i]>=death_point)*lambda)
    if(Time[i]>max(death_point))CumHazard[i] = Inf
    if(Time[i]<min(death_point))CumHazard[i] = 0
  }
  survival = exp(-CumHazard)
  list(survival=survival,lambda=lambda) 
}

# EM algorithm as in the smcure package
em1 =  function(Time,Status,X,Z,gamma,beta,emmax,eps)
{     
  w = Status	
  n = length(Status)
  s = surv_cox(Time,Status,Z,beta,w)$survival
  convergence= 1000
  i = 1
  while (convergence > eps & i < emmax){  
    uncureprob = matrix(1/(1+exp(-(gamma)%*%t(X))),ncol=1)
    survival = drop(s^(exp((beta)%*%t(Z))))
    ## E step 
    w <- Status+(1-Status)*(uncureprob*survival)/((1-uncureprob)+uncureprob*survival)
    ## M step
    logistfit<- eval(parse(text = paste("glm", "(", "w~X[,-1]",",family = quasibinomial(link='", "logit", "'",")",")",sep = "")))
    update_gamma = logistfit$coef
    update_beta = coxph(Surv(Time, Status)~Z+offset(log(w)), subset=w!=0, method="breslow")$coef
    update_s = surv_cox(Time,Status,Z,beta,w)$survival
    lambda=surv_cox(Time,Status,Z,beta,w)$lambda
    convergence = sum(c(update_gamma-gamma,update_beta-beta)^2)+sum((s-update_s)^2)
    gamma = update_gamma
    beta = update_beta 
    s = update_s
    i = i+1
  }
  em1 = list(b=gamma, beta= beta,survival=s,tau=convergence,lambda=lambda)
}

# estimates a cure logistic model via the EM algorithm
smcure1=function (lvars, ivars, data, 
                  model = "ph", link = "logit",x=NULL, Var = TRUE, emmax = 100, 
                  eps = 1e-07, nboot = 500) # x is the vector of covariates on which you condition to compare cure probabilities, needed if Var=TRUE
{
  n <- dim(data)[1]
  X <- as.matrix(cbind(rep(1, n), data[, ivars]))
  colnames(X) <- c("(Intercept)", ivars)
  
  Z <- as.matrix(data[,lvars])
  colnames(Z) <- lvars
  
  Time <- data[, 1]
  Status <- data[, 2]
  bnm <- colnames(X)
  nb <- ncol(X)
 
  betanm <- colnames(Z)
  nbeta <- ncol(Z) 
  
  nx=dim(x)[1]
  
  w <- Status
  b <- eval(parse(text = paste("glm", "(", "w~X[,-1]", ",family = quasibinomial(link='", 
                               link, "'", ")", ")", sep = "")))$coef

  beta <- coxph(Surv(Time, Status) ~ Z + offset(log(w)), 
                  subset = w != 0, method = "breslow")$coef
  
  emfit <- em1(Time, Status, X, Z, b, beta, emmax, eps)
  b <- emfit$b
  beta <- emfit$beta
  s <- emfit$Survival
  lambda=emfit$lambda
  
  p=function(gamma,x){
    r=gamma%*%t(x)
    r=1/(1+exp(-r))
    return(r)}
  
  cprob=1-as.numeric(p(b,cbind(rep(1,nx),x)))
 
  if (Var) {
    
      b_boot <- matrix(rep(0, nboot * nb), nrow = nboot)
      beta_boot <- matrix(rep(0, nboot * (nbeta)), nrow = nboot)
      cure_prob=matrix(rep(0, nboot * nx), nrow = nboot)
      iter <- matrix(rep(0, nboot), ncol = 1)
    
    
    tempdata <- cbind(Time, Status, X, Z)
    data1 <- subset(tempdata, Status == 1)
    data0 <- subset(tempdata, Status == 0)
    n1 <- nrow(data1)
    n0 <- nrow(data0)
    i <- 1
    while (i <= nboot) {
      id1 <- sample(1:n1, n1, replace = TRUE)
      id0 <- sample(1:n0, n0, replace = TRUE)
      bootdata <- rbind(data1[id1, ], data0[id0, ])
      bootX <- bootdata[, bnm]
      bootZ<- as.matrix(bootdata[,betanm])
  
      bootfit <- em1(bootdata[, 1], bootdata[, 2], bootX, 
                    bootZ,  b, beta, emmax,eps)
      b_boot[i, ] <- bootfit$b
      beta_boot[i, ] <- bootfit$beta
      
      cure_prob[i,]=1-as.numeric(p(bootfit$b,cbind(rep(1,nx),x)))

      if (bootfit$tau < eps) 
        i <- i + 1
    }
    b_var <- apply(b_boot, 2, var)
    beta_var <- apply(beta_boot, 2, var)
    b_sd <- sqrt(b_var)
    beta_sd <- sqrt(beta_var)
    p_CI=matrix(0,nx,2)
    for(j in 1:nx){
      p_CI[j,1]=quantile(cure_prob[,j],0.025)
      p_CI[j,2]=quantile(cure_prob[,j],0.975)
    }
  }
  fit <- list()
  class(fit) <- c("smcure")
  fit$b <- b
  fit$beta <- beta
  fit$cure_prob=cprob
  if (Var) {
    fit$b_var <- b_var
    fit$b_sd <- b_sd
    fit$b_zvalue <- fit$b/b_sd
    fit$b_pvalue <- (1 - pnorm(abs(fit$b_zvalue))) * 2
    fit$beta_var <- beta_var
    fit$beta_sd <- beta_sd
    fit$beta_zvalue <- fit$beta/beta_sd
    fit$beta_pvalue <- (1 - pnorm(abs(fit$beta_zvalue))) *2
    fit$p_CI=p_CI
  }
  cat(" done.\n")
  fit$call <- call
  fit$bnm <- bnm
  fit$betanm <- betanm
  fit$s <- s
  fit$lambda=lambda
  fit$Time <- Time
  fit
  printsmcure1(fit, Var)
}

# prints results of smcure1
printsmcure1=function (x, Var = TRUE, ...) 
{
  if (is.null(Var)) 
    Var = TRUE
  if (!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nCure probability model:\n")
  if (Var) {
    b <- array(x$b, c(length(x$b), 4))
    rownames(b) <- x$bnm
    colnames(b) <- c("Estimate", "Std.Error", 
                     "Z value", "Pr(>|Z|)")
    b[, 2] <- x$b_sd
    b[, 3] <- x$b_zvalue
    b[, 4] <- x$b_pvalue
  }
  if (!Var) {
    b <- array(x$b, c(length(x$b), 1))
    rownames(b) <- x$bnm
    colnames(b) <- "Estimate"
  }
  print(b)
  cat("\n")
  cat("\nFailure time distribution model:\n")
  if (Var) {
    beta <- array(x$beta, c(length(x$beta), 4))
    rownames(beta) <- x$betanm
    colnames(beta) <- c("Estimate", "Std.Error", 
                        "Z value", "Pr(>|Z|)")
    beta[, 2] <- x$beta_sd
    beta[, 3] <- x$beta_zvalue
    beta[, 4] <- x$beta_pvalue
  }
  if (!Var) {
    beta <- array(x$beta, c(length(x$beta), 1))
    rownames(beta) <- x$betanm
    colnames(beta) <- "Estimate"
  }
  print(beta)
  cat("\n")
  cat("\nConditional cure probabilites :\n")
  if (Var) {
    b <- array(x$cure_prob, c(length(x$cure_prob), 3))
    rownames(b) <- c()
    colnames(b) <- c("Estimate", "lower bound CI", 
                     "upper bound CI")
    b[, 2] <- x$p_CI[,1]
    b[, 3] <- x$p_CI[,2]
  }
  if (!Var) {
    b <- array(x$cure_prob, c(length(x$cure_prob), 1))
    rownames(b) <- c()
    colnames(b) <- "cond cure prob"
  }
  print(b)
  invisible(x)
}

# function that computes CI for difference in cure rates
cure_CI=function (lvars, ivars, data1,data2, 
                  model = "ph", link = "logit",x, emmax = 100, 
                  eps = 1e-07, nboot = 500,g1,b1,g2,b2) # x is the vector of covariates on which you condition to compare cure probabilities, needed if Var=TRUE
{
  nx=dim(x)[1]
  cure_prob1=matrix(rep(0, nboot * nx), nrow = nboot)
  cure_prob2=matrix(rep(0, nboot * nx), nrow = nboot)
  boot_cure_diff=matrix(rep(0, nboot * nx), nrow = nboot)
  
  i=1
  while(i<=nboot){
    
    data11<-subset(data1,data1[,2]==1)
    data10<-subset(data1,data1[,2]==0)
    n11<-nrow(data11)
    n10<-nrow(data10) 
    id11<-sample(1:n11,n11,replace=TRUE)
    id10<-sample(1:n10,n10,replace=TRUE)
    bootdata1<-rbind(data11[id11,],data10[id10,])
    
    data21<-subset(data2,data2[,2]==1)
    data20<-subset(data2,data2[,2]==0)
    n21<-nrow(data21)
    n20<-nrow(data20) 
    id21<-sample(1:n21,n21,replace=TRUE)
    id20<-sample(1:n20,n20,replace=TRUE)
    bootdata2<-rbind(data21[id21,],data20[id20,])
    
    p=function(gamma,x){
      r=gamma%*%t(x)
      r=1/(1+exp(-r))
      return(r)}
    
    bootX1 <- as.matrix(cbind(rep(1,n11+n10),bootdata1[,ivars]))
    bootZ1<- as.matrix(bootdata1[,lvars])
    bootfit1 <- em1(bootdata1[, 1], bootdata1[, 2], bootX1, 
                   bootZ1,  g1, b1, emmax,eps)
    cure_prob1[i,]=1-as.numeric(p(bootfit1$b,cbind(rep(1,nx),x)))
  
    bootX2 <-as.matrix(cbind(rep(1,n21+n20),bootdata2[,ivars]))
    bootZ2<- as.matrix(bootdata2[,lvars])
    bootfit2 <- em1(bootdata2[, 1], bootdata2[, 2], bootX2, 
                    bootZ2,  g2, b2, emmax,eps)
    cure_prob2[i,]=1-as.numeric(p(bootfit2$b,cbind(rep(1,nx),x)))
    
    boot_cure_diff[i,]=cure_prob1[i,]-cure_prob2[i,]
    
    if(bootfit1$tau<1e-7 & bootfit2$tau<1e-7){
      i=i+1
    }
  }
  p_CI=matrix(NA,nx,2)
  for(j in 1:nx){
    p_CI[j,1]=quantile(boot_cure_diff[,j],0.025)
    p_CI[j,2]=quantile(boot_cure_diff[,j],0.975)
  }
  return(p_CI)

}

# main function

# Input: Data1, Data2 - two datasets (corresponding to the two groups to be compared), 
#                       first column is survival time ("Y"), second column is censoring indicator ("status"), then other covariates
#       ivars, lvars  - names of covariates that are used in the incidence and latency components
#                x,z  - matrix of covariate values (for incidence and latency) on which we condition, each row is a different choice of covariates, ncol is equal to the number of cov for the latency model
#           cure_comp - logical: if TRUE computes also parameters related to the cure component
#      beta.hat1,gamma.hat1,beta.hat2,gamma.hat2 - initial estimates for the parameters of the two models, can also be NULL

# Output:  m - a vector containing the difference in mean survival time for the uncured conditionally on each covariate
#         conv1, conv2 - logical showing whether EM algorithm converged
#         g1,b1,g2,b2 - estimates of the parameters
#         If cure_comp==TRUE returns also:
#            p1  - cure probability in sample 1 conditional on x
#            p2  - cure probability in sample 2 conditional on x
#            p_CI_1 - 95% CI for p1 using bootstrap
#            p_CI_2 - 95% CI for p2 using bootstrap

estimation=function(Data1,Data2,ivars,lvars,z,x,cure_comp,beta.hat1,gamma.hat1,beta.hat2,gamma.hat2){
  
  n1=dim(Data1)[1]
  n2=dim(Data2)[1]
  
  nz=dim(z)[1]  # number of covariates on which we condition
  
  # extract from the datasets the matrices of covariates for incidence and latency
  Z1= as.matrix(Data1[,lvars])
  X1= as.matrix(Data1[,ivars])
  
  Z2= as.matrix(Data2[,lvars])
  X2= as.matrix(Data2[,ivars])
  
  # computes initial estimates (needed for the EM algorithm) if they are not provided
  if(is.null(beta.hat1)){
    beta.hat1 <- coxph(Surv(Y,status)~Z1,subset = status!=0,data = Data1,method = "breslow")$coef 
    gamma.hat1 <- glm(status~X1,family=binomial(link=logit),data=Data1)$coefficients
    beta.hat2 <- coxph(Surv(Y,status)~Z2,subset = status!=0,data = Data2,method = "breslow")$coef 
    gamma.hat2 <- glm(status~X2,family=binomial(link=logit),data=Data2)$coefficients
  }
  
  # estimate parameters and mean survival times for each dataset
  
  MCM1=em1(Data1$Y,Data1$status,Z=Z1,X=cbind(rep(1,n1),X1),gamma.hat1,beta.hat1,emmax=100,eps=1e-7) 
  
  est_gamma1=MCM1$b
  est_beta1=MCM1$beta
  event_times1=sort(Data1$Y)
  est_s1=sort(MCM1$survival,decreasing = TRUE)
  
  E1=rep(NA,nz)
  for(j in 1:nz){
     E1[j]=sum(diff(c(0,event_times1))*(c(1,est_s1[-n1])^(exp(c(est_beta1%*%z[j,])))))
  }
  
  MCM2=em1(Data2$Y,Data2$status,Z=Z2,X=cbind(rep(1,n2),X2),gamma.hat2,beta.hat2,emmax=100,eps=1e-7) 
  
  est_gamma2=MCM2$b
  est_beta2=MCM2$beta
  event_times2=sort(Data2$Y)
  est_s2=sort(MCM2$survival,decreasing = TRUE)
  
  E2=rep(NA,nz)
  for(j in 1:nz){
    E2[j]=sum(diff(c(0,event_times2))*(c(1,est_s2[-n2])^(exp(c(est_beta2%*%z[j,])))))
  }
  
  #difference of conditional mean survival times
  m=E1-E2
  
  if(cure_comp==FALSE){
    return(list(m=m,conv1=MCM1$tau,conv2=MCM2$tau,g1=est_gamma1,b1=est_beta1,g2=est_gamma2,b2=est_beta2))
  }
  
  if(cure_comp==TRUE){
    
    p=function(gamma,x){
      r=gamma%*%t(x)
      r=1/(1+exp(-r))
      return(r)}
    
    nx=dim(x)[1]
    cure_prob1=1-as.numeric(p(est_gamma1,cbind(rep(1,nx),x)))
    cure_prob2=1-as.numeric(p(est_gamma2,cbind(rep(1,nx),x)))
    
    p_CI=cure_CI(lvars,ivars,data1=Data1,data2=Data2,model="ph",link = "logit",x=x,emmax=100,nboot=500,g1=est_gamma1,b1=est_beta1,g2=est_gamma2,b2=est_beta2)
    
    return(list(m=m,conv1=MCM1$tau,conv2=MCM2$tau,g1=est_gamma1,b1=est_beta1,g2=est_gamma2,b2=est_beta2,p1=cure_prob1,p2=cure_prob2,p_CI=p_CI))
  }
}

# performs estimation of m for B bootstrap samples separately on two datasets conditionally on covariates z (given initial parameters that are used for the estimation in the previous function)
# z is a matrix containing in each row the covariate values 

bootstrap=function(Data1,Data2,ivars,lvars,b1,g1,b2,g2,z,B){
  
  ncov=dim(z)[1]
  
  boot_m=matrix(data=NA,nrow=ncov,ncol=B)
   
  i=1
  while(i<=B){
    
    data11<-subset(Data1,Data1[,2]==1)
    data10<-subset(Data1,Data1[,2]==0)
    n11<-nrow(data11)
    n10<-nrow(data10) 
    id11<-sample(1:n11,n11,replace=TRUE)
    id10<-sample(1:n10,n10,replace=TRUE)
    bootdata1<-rbind(data11[id11,],data10[id10,])
    
    data21<-subset(Data2,Data2[,2]==1)
    data20<-subset(Data2,Data2[,2]==0)
    n21<-nrow(data21)
    n20<-nrow(data20) 
    id21<-sample(1:n21,n21,replace=TRUE)
    id20<-sample(1:n20,n20,replace=TRUE)
    bootdata2<-rbind(data21[id21,],data20[id20,])
    
    boot_results=estimation(bootdata1,bootdata2,ivars,lvars,x=NULL,z,cure_comp=FALSE,b1,g1,b2,g2)
   
    boot_m[,i]=boot_results$m
    
    if(boot_results$conv1<1e-7 & boot_results$conv2<1e-7){
      i=i+1
    }
  }
  return(boot_m)
}

# end function


#######################
# Example smcure data #
#######################

data(e1684)  
d=na.omit(e1684) 

Data1=d[which(d$TRT==0),c(2,3,4,5)]  # the two datasets with survival time, censoring indicator in the first two columns and other covariates afterwards
Data2=d[which(d$TRT==1),c(2,3,4,5)]

colnames(Data1)=c("Y","status","Age","Sex")  #keep the same names for the first two columns
colnames(Data2)=c("Y","status","Age","Sex")

ivars=c("Age","Sex") # covariates used in the incidence model
lvars=c("Age","Sex") # covariates used in the latency model

# matrix of covariate values on which we condition for the latency, first column is age (centered to the mean), second is sex
z=c(0,1)
z=rbind(z,c(0,0))
z=rbind(z,c(10,0))
z=rbind(z,c(10,1))
z=rbind(z,c(-10,0))
z=rbind(z,c(-10,1))
z=rbind(z,c(20,0))
z=rbind(z,c(20,1))
z=rbind(z,c(-20,0))
z=rbind(z,c(-20,1))

# matrix of covariate values on which we condition for the incidence, 
x=z

#############################################################################
# the code below is standard (not depending on the dataset)
# preliminary analysis of the data

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

est_res=estimation(Data1,Data2,ivars,lvars,z=z,x=x,cure_comp = TRUE,beta.hat1=NULL,gamma.hat1=NULL,beta.hat2=NULL,gamma.hat2=NULL)

est_gamma1=est_res$g1
est_beta1=est_res$b1

est_gamma2=est_res$g2
est_beta2=est_res$b2

cure_prob1=est_res$p1
cure_prob_CI1=est_res$p_CI_1
cure_prob2=est_res$p2
cure_prob_CI2=est_res$p_CI_2

cure_prob1-cure_prob2  # difference in conditional cure probabilities
cure_prob_CI1-cure_prob_CI2[,c(2,1)]  # confidence interval for difference

m=est_res$m

B=100  # number of bootstrap samples for estimating variance of m
boot_results=bootstrap(Data1,Data2,ivars,lvars,est_beta1,est_gamma1,est_beta2,est_gamma2,z,B)

# vector of standard errors for m conditional on each covariate value using bootstrap
sds=apply(boot_results,1,sd)

#######################
# Asymptotic approach #
#######################

# asymptotic pvalue for testing H0: m_z=0 aginst H1: m_z!=0

as_pval=rep(NA,(dim(z)[1]))
for(j in 1:(dim(z)[1])){
  as_pval[j]=2*(1-pnorm(abs(m[j]/sds[j])))
}

# asymptotic 95% CI for m
as_CI=matrix(0,dim(z)[1],2)
for(j in 1:(dim(z)[1])){
  as_CI[j,1]=m[j]-qnorm(0.975)*sds[j]
  as_CI[j,2]=m[j]+qnorm(0.975)*sds[j]
}

########################
# Permutation approach #
########################

K=500  #nr of permutations
n=n1+n2
pooled_data=rbind(Data1,Data2)

m_p=matrix(NA,nrow=dim(z)[1],ncol=K)
sd_p=matrix(NA,nrow=dim(z)[1],ncol=K)

for(j in 1:K){
  id<-sample(1:n,n,replace=FALSE)
  new_Data1=pooled_data[id[1:n1],]
  new_Data2=pooled_data[id[(n1+1):n],]
  
  perm_est=estimation(new_Data1,new_Data2,ivars,lvars,x=NULL,z,cure_comp=FALSE,NULL,NULL,NULL,NULL)
  
  m_p[,j]=perm_est$m
  
  boot_perm=bootstrap(new_Data1,new_Data2,ivars,lvars,perm_est$b1,perm_est$g1,perm_est$b2,perm_est$g2,z,B)
  sd_p[,j]=apply(boot_perm,1,sd)
}

# quantiles for the permutation distribution 
q=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975)
quantiles=matrix(0,dim(z)[1],7)
for(j in 1:(dim(z)[1])){
   quantiles[j,]=c(quantile(m_p[j,]/sd_p[j,],q))
}

# pvalue with permutation approach for testing H0: m_z=0 aginst H1: m_z!=0
perm_pval=rep(NA,dim(z)[1])
for(j in 1:dim(z)[1]){
  perm_pval[j]=(length(which(m_p[j,]/sd_p[j,]<(-abs(m[j]/sds[j]))))+length(which(m_p[j,]/sd_p[j,]>abs(m[j]/sds[j]))))/K
}

# 95% CI for m with permutation approach
perm_CI=matrix(0,dim(z)[1],2)
for(j in 1:(dim(z)[1])){
 perm_CI[j,1]=m[j]-quantiles[j,7]*sds[j]
 perm_CI[j,2]=m[j]-quantiles[j,1]*sds[j]
}

