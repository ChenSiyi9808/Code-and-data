load(file='C:\\Users\\syephd\\Desktop\\Bivariate exponential model\\Bivariate exponential model\\Code\\SimulationIG.Rdata')
Data<-read.csv('C:\\Users\\syephd\\Desktop\\Bivariate exponential model\\Bivariate exponential model\\Code\\Crack.meeker.csv')
library(tweedie)
rm (list = ls ())
f<-function(S){
lambda.t<-function(t,Gam){
  tau<-matrix(NA,ncol=length(Gam),nrow=length(t))
  for(j in 1:length(Gam)){
    tau[,j]<-t^Gam[j] 
  }
  
  return(tau)
}
Simu.deg<-Simulation[[S]]
n<-50;m<-50
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
alpha.s<-c(1.5,1)
sigma.s<-0.12
lambda.s<-c(10,10)
Gam.s<-c(1,1.5)
Delta<-0.1;L<-1
Para<-list()
Emu<-matrix(0,nrow=n,ncol=1)
Emu2<-matrix(0,nrow=n,ncol=1)
library(parallel)
library(doParallel)
time_start<-Sys.time()
AIC<-0
L<-1
Sigma.s<-0
while(Delta>1e-5){
  Para[[L]]<-list(sigma.s=sigma.s,alpha.s=alpha.s,lambda.s=lambda.s,Gam.s=Gam.s,AIC=AIC)
  ##############E-step:
  delta.t.s<-apply(lambda.t(t,Gam.s),2,diff)
  t.s<-lambda.t(t,Gam.s)
  for(i in 1:n){
    x<-cumsum(Y[which(Y$unit==i),])#ith individual degradation
    y<-x[nrow(x),3:4]
    Sigma.s<-(sum((lambda.s* y)/alpha.s^2)+sigma.s^(-2))
    Emu[i]<-(sum((lambda.s* t.s[nrow(t.s),])/alpha.s)+sigma.s^(-2))/Sigma.s
    Emu2[i]<-1/Sigma.s+(Emu[i])^2
  }
  
  # plot(par[,3])
  #sigma.s<-sqrt(mean(apply((par[501:1000,]-1)^2,2,mean)))
  
  #############M-step
  sigma.s<-sqrt((sum(Emu2)-2*sum(Emu)+n)/n)
  ####################################alpha
  s1<-matrix(0,nrow=n,ncol=2)
  s2<-matrix(0,nrow=n,ncol=2)
  for(i in 1:n){
    x<-as.matrix(Simu.deg[which(Simu.deg$unit==i),3:4])
    s1[i,]<-(x[nrow(x),])*Emu2[i]
    s2[i,]<-t.s[nrow(t.s),]*Emu[i]
  }
  alpha.s<-c(sum(s1[,1])/sum(s2[,1]),sum(s1[,2])/sum(s2[,2]))
  ##########################
  Z<-Y[,3:4]
  for(i in 1:n){
    x<-as.matrix(Z[Y$unit==i,])
    z<-matrix(0,nrow=m,ncol=2)
    for(k in 1:m){
      z[k,]<-(x[k,]^2*Emu2[i]-2*x[k,]*delta.t.s[k,]*alpha.s*Emu[i]+alpha.s^2*delta.t.s[k,]^2)/((alpha.s^2*x[k,]))
    }
    Z[Y$unit==i,]<-z
  }
  lambda.s<-c((n*m)/(sum(Z[,1])),(n*m)/(sum(Z[,2])))
  ###########################
  Gam.f<-function(j){
    fit<-function(Gam){
      s<-matrix(0,nrow=n,ncol=1)
      for(i in 1:n){
        x<-Y[which(Y$unit==i),]
        tau<-apply(lambda.t(t,Gam),2,diff)
        s[i]<-sum(2*log(tau)-(x[,j+2]^2*Emu2[i]-2*alpha.s[j]*x[,j+2]*tau*Emu[i]+alpha.s[j]^2*tau^2)/((alpha.s[j]^2*x[,j+2])/lambda.s[j]))
      }
      s<--sum(s)
      return(s)
    }
    Fit<-optim(par=Gam.s[j],fit,method='Brent',lower=0,upper=3)
    return(Fit$par)
  }
  library(doParallel)
  Gam.s<-unlist(foreach(i=1:length(Gam.s))%do% Gam.f(i))
  
  Delta<-max(c(abs(((Para[[L]]$alpha.s-alpha.s)/Para[[L]]$alpha.s)),abs(Para[[L]]$sigma.s-sigma.s)/sigma.s,abs(Para[[L]]$alpha.s-alpha.s)/alpha.s,
               abs(Para[[L]]$lambda.s-lambda.s)/Para[[L]]$lambda.s,abs(Para[[L]]$Gam.s-Gam.s)/Para[[L]]$Gam.s))
  L<-L+1
}
return(Para[[L-1]])
}
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code执行时间：',round(exc_time,2),'mins'))
for(i in 1:n){
  x<-Y[which(Y$unit==i),3:4]
  for(j in 1:2){
    for(k in 1:m){
      x[k,j]<-0.5*log((lambda.s[j]*delta.t.s[k,j]^2)/(2*pi*x[k,j]^3))-(x[k,j]^2*Emu2[i]-2*x[k,j]*alpha.s[j]*delta.t.s[k,j]*Emu[i]+alpha.s[j]^2*delta.t.s[k,j]^2)/((2*alpha.s[j]^2*x[k,j])/lambda.s[j])
    }
  }
  Y[which(Y$unit==i),3:4]<-x
}
L<--0.5*log(2*pi*sigma.s)*n-(sum(Emu2)-2*sum(Emu)+n)/(2*sigma.s^2)+sum(Y[,3:4])
Lik<--2*L+2*7
