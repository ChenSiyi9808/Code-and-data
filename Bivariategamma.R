#load(file='F:\\degradation analysis\\Bivariate exponential model\\Code\\SimulationG.Rdata')
#Data<-read.csv('F:\\degradation analysis\\Bivariate exponential model\\Code\\Crack.meeker.csv')
load(file='C:\\Users\\DELL\\Desktop\\Bivariate exponential model\\Bivariate exponential model\\Code\\SimulationG.Rdata')
Data<-read.csv('C:\\Users\\DELL\\Desktop\\Bivariate exponential model\\Bivariate exponential model\\Code\\Crack.meeker.csv')
library(tweedie)
lambda.t<-function(t,Gam){
  tau<-matrix(NA,ncol=length(Gam),nrow=length(t))
  for(j in 1:length(Gam)){
    tau[,j]<-t^Gam[j] 
  }
  
  return(tau)
}
Simu.deg<-(Data)
n<-6;m<-5
t<-matrix(Simu.deg$time,ncol=n)
t<-c(t[,1])
##data processing
for(i in 1:n){
  x<-Simu.deg[which(Simu.deg$unit==i),]
  x[2:nrow(x),3:(3+1)]<-apply(x[1:nrow(x),3:(3+1)],2,diff)
  x[1,3:(3+1)]<-NA
  if(i==1){
    Y<-x
  }else{
    y<-x
    Y<-rbind(Y,x)
  }
}
Y<-na.omit(Y)
eta.s<-c(1)
alpha.s<-c(0.03,0.02)
beta.s<-5
lambda.s<-c(0.01,0.05)
Gam.s<-c(1.3,1.5)
Delta<-0.1;L<-1
Para<-list()
Emu<-matrix(0,nrow=n,ncol=1)
Elnmu<-matrix(0,nrow=n,ncol=1)
library(parallel)
library(doParallel)
time_start<-Sys.time()
AIC<-0
L<-1
shape.s<-0
scale.s<-0
while(Delta>1e-5){
  Para[[L]]<-list(beta.s=beta.s,alpha.s=alpha.s,lambda.s=lambda.s,Gam.s=Gam.s,AIC=AIC)
  ##############E-step:
  delta.t.s<-apply(lambda.t(t,Gam.s),2,diff)
  t.s<-lambda.t(t,Gam.s)
  for(i in 1:n){
    x<-cumsum(Y[which(Y$unit==i),])#ith individual degradation
    y<-x[nrow(x),3:4]
    shape.s<-sum(lambda.s*t.s[nrow(t.s),])+beta.s
    scale.s<-sum((lambda.s*y)/alpha.s)+beta.s
    Emu[i]<-shape.s/scale.s
    Elnmu[i]<-digamma(shape.s)-log(scale.s)
  }
  
  # plot(par[,3])
  #sigma.s<-sqrt(mean(apply((par[501:1000,]-1)^2,2,mean)))
  
  #############M-step
  beta.f<-function(beta){
    z<-beta*log(beta)-log(gamma(beta))+(beta-1)*mean(Elnmu)-beta*mean(Emu)
    return(-z)
  }
  Fit<-optim(par=beta.s,beta.f,method='Brent',lower=0,upper=10)
  beta.s<-Fit$par
  ####################################alpha
  s1<-matrix(0,nrow=n,ncol=2)
  s2<-matrix(0,nrow=n,ncol=2)
  for(i in 1:n){
    x<-as.matrix(Simu.deg[which(Simu.deg$unit==i),3:4])
    s1[i,]<-(x[nrow(x),])*Emu[i]
    s2[i,]<-t.s[nrow(t.s),]
  }
  alpha.s<-c(sum(s1[,1])/sum(s2[,1]),sum(s1[,2])/sum(s2[,2]))
  ##########################
  lambda.f<-function(j){
    fit<-function(lambda){
      s<-matrix(0,nrow=n,ncol=1)
      for(i in 1:n){
        x<-Y[which(Y$unit==i),]
        tau<-apply(lambda.t(t,Gam.s[j]),2,diff)
        s[i]<-sum(lambda*tau*(Elnmu[i]+log(lambda)-log(alpha.s[j]))-log(gamma(lambda*tau))+(lambda*tau-1)*log(x[,j+2])-((Emu[i]*lambda)/(alpha.s[j]))*x[,j+2])
      }
      s<--sum(s)
      return(s)
    }
    Fit<-optim(par=lambda.s[j],fit,method='Brent',lower=0,upper=2)
    return(Fit$par)
  }
  library(doParallel)
  lambda.s<-unlist(foreach(i=1:length(alpha.s))%do% lambda.f(i))
  ###########################
  Gam.f<-function(j){
    fit<-function(Gam){
      s<-matrix(0,nrow=n,ncol=1)
      for(i in 1:n){
        x<-Y[which(Y$unit==i),]
        tau<-apply(lambda.t(t,Gam),2,diff)
        s[i]<-sum(lambda.s[j]*tau*(Elnmu[i]+log(lambda.s[j])-log(alpha.s[j]))-log(gamma(lambda.s[j]*tau))+(lambda.s[j]*tau-1)*log(x[,j+2]))
      }
      s<--sum(s)
      return(s)
    }
    Fit<-optim(par=Gam.s[j],fit,method='Brent',lower=0,upper=2)
    return(Fit$par)
  }
  library(doParallel)
  Gam.s<-unlist(foreach(i=1:length(Gam.s))%do% Gam.f(i))
  
  Delta<-max(c(abs(((Para[[L]]$alpha.s-alpha.s)/Para[[L]]$alpha.s)),abs(Para[[L]]$beta.s-beta.s)/beta.s,abs(Para[[L]]$alpha.s-alpha.s)/alpha.s,
               abs(Para[[L]]$lambda.s-lambda.s)/Para[[L]]$lambda.s,abs(Para[[L]]$Gam.s-Gam.s)/Para[[L]]$Gam.s))
  L<-L+1
}
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code执行时间：',round(exc_time,2),'mins'))
for(i in 1:n){
  x<-Y[which(Y$unit==i),3:4]
  for(j in 1:2){
    for(k in 1:m){
      x[k,j]<-lambda.s[j]*delta.t.s[k,j]*(Elnmu[i]+log(lambda.s[j])-log(alpha.s[j]))-log(gamma(lambda.s[j]*delta.t.s[k,j]))+(lambda.s[j]*delta.t.s[k,j]-1)*log(x[k,j])-(Emu[i]*lambda.s[j]*x[k,j])/alpha.s[j]
    }
  }
  Y[which(Y$unit==i),3:4]<-x
}
L<-beta.s*log(beta.s)*n-log(gamma(beta.s))*n+(beta.s-1)*sum(Elnmu)-beta.s*sum(Emu)+sum(Y[,3:4])
Lik<--2*L+2*7
