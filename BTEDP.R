rm (list = ls ())
#load Data

library(ggplot2)
library(tweedie)
lambda.t<-function(t,Gam){
  tau<-matrix(NA,ncol=length(Gam),nrow=length(t))
  for(j in 1:length(Gam)){
    tau[,j]<-t^Gam[j] 
  }
  return(tau)
}

S<-776
Simu.deg<-Data#Simulation[[S]]
ggplot(data=Simu.deg,aes(x=time,y=pc.1,fill=factor(unit)))+
  geom_point()+geom_line()+
  geom_line(aes(x=time,y=pc.2,color="red"))+
  #geom_line(aes(x=time,y=pc.3,color="blue"))+
  scale_fill_discrete(guide=FALSE)+guides(color=guide_legend(title=NULL))
n<-10;m<-9
Psedo_failure<-matrix(0,nrow=n,ncol=2)
Df<-c(0.7,0.4)
eta.s<-matrix(0,nrow=n,ncol=2)
Gam.s<-matrix(0,nrow=n,ncol=2)
for(i in 1:n){
  for(j in 1:2){
    Fit.init<-nls(y~a*t^b,data=data.frame(t=Simu.deg[which(Simu.deg$unit==i),2],y=Simu.deg[which(Simu.deg$unit==i),j+2]),start=list(a=10,b=1))
    eta.s[i,j]<-coefficients(Fit.init)[1]
    Gam.s[i,j]<-coefficients(Fit.init)[2]
    Psedo_failure[i,j]<-(Df[j]/eta.s[i,j])^(1/Gam.s[j])
  }
}

apply(eta.s,2,mean)
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
alpha.s<-c(13.678996,7.274991)
sigma.s<-0.1251855
rho.s<-c(3.171255,1.058856)
lambda.s<-c(107199,511)
Gam.s<-c(1.346468,1.282084)

alpha.s<-c(1.5,1)
sigma.s<-0.1
rho.s<-c(2.7,2.2)
lambda.s<-c(10,10)
Gam.s<-c(1,1.5)
alpha.s<-c(0.03,0.02)
sigma.s<-0.15
rho.s<-c(2,2)
lambda.s<-c(1,2)
Gam.s<-c(1.5,1.5)
Delta<-0.1;L<-1
Para<-list()

library(parallel)
library(doParallel)
#time_start<-Sys.time()
AIC<-0

while(L<201){
  Para[[L]]<-list(sigma.s=sigma.s,alpha.s=alpha.s,lambda.s=lambda.s,Gam.s=Gam.s,rho.s=rho.s,AIC=AIC)
  ##############E-step:MCMC algorithm
  init<-matrix(1,nrow=1,ncol=n)#rnorm(n,mean=1,sd=0.1)#
  par<-matrix(0,nrow=200,ncol=n)
  delta.t.s<-apply(lambda.t(t,Gam.s),2,diff)
  for(b in 1:200){
    init_hat<-init
    for(i in 1:n){
      y<-Y[Y$unit==i,]
      init_hat[i]<-(rnorm(1,mean=init[i],sd=0.05))
      z<-runif(1,min=0,max=1)
      Alpha<-min(1,exp(h1(y,delta.t.s,init_hat[i],alpha.s,lambda.s,rho.s,sigma.s)-h1(y,delta.t.s,init[i],alpha.s,lambda.s,rho.s,sigma.s)))
      if(z>Alpha){
        init[i]<-init[i]
      }else{
        init[i]<-init_hat[i]
      }
      
    }
    par[b,]<-init
  }
  plot(par[,1])
  #sigma.s<-sqrt(mean(apply((par[501:1000,]-1)^2,2,mean)))
  #############M-step
  sigma.s<-sqrt(mean(apply((par[101:200,]-1)^2,2,mean)))
  ###########################
  ###################
  lambda.f<-function(j){
    return(lambda.F(alpha.s[j],rho.s[j],Gam.s[j],Y[,j+2],t,par))
  }
  lambda.s<-unlist(foreach(i=1:length(lambda.s))%do% lambda.f(i))
  ##########################
  rho.f<-function(j){
    fit<-function(rho){
      return(-para_M_step(alpha.s[j],lambda.s[j],rho,Gam.s[j],Y[,j+2],t,par))
    }
    Fit<-optim(par=rho.s[j],fit,method='Brent',lower=0,upper=10)
    return(Fit$par)
  }
  
  
  rho.s<-unlist(foreach(i=1:length(rho.s))%do% rho.f(i))
  ####################################alpha
  alpha.f<-function(j){
    fit<-function(alpha){
      return(-para_M_step(alpha,lambda.s[j],rho.s[j],Gam.s[j],Y[,j+2],t,par))
    }
    Fit<-optim(par=alpha.s[j],fit,method='Brent',lower=0,upper=20)
    return(Fit$par)
  }
  alpha.s<-unlist(foreach(i=1:length(alpha.s))%do% alpha.f(i))
  
  ###########################
  Gam.f<-function(j){
    fit<-function(Gam){
      return(-para_M_step(alpha.s[j],lambda.s[j],rho.s[j],Gam,Y[,j+2],t,par))
    }
    Fit<-optim(par=Gam.s[j],fit,method='Brent',lower=0,upper=2)
    return(Fit$par)
  }
  Gam.s<-unlist(foreach(i=1:length(Gam.s))%do% Gam.f(i))
  ##########################################
  
  AIC<--2*aic(alpha.s,lambda.s,rho.s,Gam.s,sigma.s,Y,delta.t.s,par[101:200,],n,m)+2*9
  Delta<-max(c(abs(((Para[[L]]$alpha.s-alpha.s)/Para[[L]]$alpha.s)),
               abs(Para[[L]]$lambda.s-lambda.s)/Para[[L]]$lambda.s,abs(Para[[L]]$Gam.s-Gam.s)/Para[[L]]$Gam.s))
  L<-L+1
}


alpha<-matrix(0,nrow=L,ncol=2)
for(i in 1:(L-1)){
  alpha[i,]<-Para[[i]]$AIC
  
}
plot(alpha[2:(L-1),1])
mean(rho.s[1000:2000,2])



alpha.s^rho.s/lambda.s
alpha^rho*lambda
Data<-read.csv('/home/zhanghu1/BEDM/Crack.meeker.csv')

eta.s1<-matrix(0,nrow=500,ncol=2)
lambda.s1<-matrix(0,nrow=500,ncol=2)
Gam.s1<-matrix(0,nrow=500,ncol=2)
rho.s1<-matrix(0,nrow=500,ncol=2)
sigma.s1<-matrix(0,nrow=500,ncol=1)
AIC1<-matrix(0,nrow=500,ncol=1)
for(i in 1:500){
  eta.s1[i,]<-Para[[i]]$alpha.s
  lambda.s1[i,]<-Para[[i]]$lambda.s
  Gam.s1[i,]<-Para[[i]]$Gam.s
  rho.s1[i,]<-Para[[i]]$rho.s
  sigma.s1[i]<-Para[[i]]$sigma.s
  AIC1[i]<-Para[[i]]$AIC
}
par<-list(sigma.s=mean(sigma.s1[251:500]),
          alpha.s=apply(eta.s1[251:500,],2,mean),
          lambda.s=apply(lambda.s1[251:500,],2,mean),
          Gam.s=apply(Gam.s1[251:500,],2,mean),
          rho.s=apply(rho.s1[251:500,],2,mean),
          AIC=mean(AIC1[251:500]))
plot(eta.s1[,1])
#############################################
load('/home/zhanghu1/BEDM/Simulationmeeker.Rdata')
#save(Para,file='/home/zhanghu1/BEDM/crackparatrain.Rdata')
#load('/home/zhanghu1/BEDM/simulationsample.Rdata')
Para<-list()
library(doParallel)
cl <- makeCluster(60)
registerDoParallel(cl)
time_start<-Sys.time()
Para<- foreach(r=1:1000) %dopar% f(r)
stopCluster(cl)
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code执行时间：',round(exc_time,2),'mins'))
save(Para, file='/home/zhanghu1/BEDM/Simulationparameeeker2000.Rdata')
################################################
load('F:\\degradation analysis\\Bivariate exponential model\\Code\\模拟结果\\Simulationpara1010.Rdata')
#save(Para,file='F:\\degradation analysis\\Bivariate exponential model\\Code\\paracracktrain.Rdata')

N<-1000
eta.s<-matrix(0,nrow=N,ncol=2)
lambda.s<-matrix(0,nrow=N,ncol=2)
Gam.s<-matrix(0,nrow=N,ncol=2)
rho.s<-matrix(0,nrow=N,ncol=2)
sigma.s<-matrix(0,nrow=N,ncol=1)
AIC<-matrix(0,nrow=N,ncol=1)
for(i in 1:N){
  eta.s[i,]<-Para[[i]]$alpha.s
  lambda.s[i,]<-Para[[i]]$lambda.s
  Gam.s[i,]<-Para[[i]]$Gam.s
  rho.s[i,]<-Para[[i]]$rho.s
  sigma.s[i]<-Para[[i]]$sigma.s
  AIC[i]<-Para[[i]]$AIC
}

a<-rho.s[,2]
plot(a)
CI<-matrix(NA,nrow=1000,ncol=3)
CI[,1]<-a-1.96*sd(a)
CI[,2]<-2.4
CI[,3]<-a+1.96*sd(a)
index<-function(y){
  if(y[1]<=y[2]&y[2]<=y[3]){
    x=1
  }else{x=0}
  return(x)
}
sum(apply(CI,1,index))/1000
quantile(a,c(0.025,0.975))[2]-quantile(a,c(0.025,0.975))[1]


(mean((Gam.s[,1])))
(mean((lambda.s[,2])))
(mean((rho.s[,2])))
(mean((sigma.s)))
plot(lambda.s[,1])

sqrt(mean((eta.s[,2]-1)^2))
sqrt(mean((Gam.s[,2]-1.5)^2))
sqrt(mean((lambda.s[,2]-10)^2))
sqrt(mean((rho.s[,2]-2.2)^2))
sqrt(mean((sigma.s-0.1)^2))

sd(lambda.s[,1])
out<-boxplot(rho.s[,1])$out
rho.s[which(rho.s[,1]%in%out ),1]<-NA
lambda.s[which(rho.s[,1]%in%out ) ,1]<-NA
out<-boxplot(lambda.s[,1])$out
lambda.s[which(lambda.s[,1]%in%out ),1]<-NA

rho.s[]
rho.s[which(rho.s[,2]<0.01 ),2]<-NA
boxplot(lambda.s[,1])
lambda.s<-na.omit(lambda.s)
quantile(na.omit(lambda.s[,1]),c(0.025,0.975))
length(na.omit(lambda.s[,1]))
help("boxplot")
plot(lambda.s[,2])
plot(rho.s[,2])
##############################################
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code执行时间：',round(exc_time,2),'mins'))

Data$fitted1<-alpha.s[1]*Data$time^Gam.s[1]
ggplot(data=Data,aes(x=time,y=pc.1,fill=factor(unit)))+
  geom_point()+geom_line()+
  #geom_line(aes(x=time,y=pc.2,color="red"))+
  geom_line(aes(x=time,y=fitted1,color="blue"))+
  scale_fill_discrete(guide=FALSE)+guides(color=guide_legend(title=NULL))
library(ggplot2)

for(i in 1:n){
  x<-Y[which(Y$unit==i),3:4]
  for(j in 1:2){
    Q<-0
    for(k in 1:m){
      Q<-d(x[k,j],delta.t.s[k],alpha.s[j]*par[101:200,i],rho.s[j])
      x[k,j]<-0.5*log(lambda.s[j]/(2*pi))-(1-rho.s[j])/2*log(delta.t.s[k,j])-rho.s[j]/2*log(x[k,j])-(lambda.s[j]*delta.t.s[k,j])/2*(mean(Q))
    }
  }
  Y[which(Y$unit==i),3:4]<-x
  P[i,1]<-mean(log(dnorm(par[101:200,i],mean=1,sd=sigma.s)))
}
-2*(sum(Y[,3:4])+sum(P))+2*9
save(R,file='F:\\degradation analysis\\Bivariate exponential model\\Code\\Reliabity.Rdata')
###############################################################
load('F:\\degradation analysis\\Bivariate exponential model\\Code\\Crack.meeker.csv')
#PC图像绘制
library(ggplot2)
library('viridis')
library(dplyr)
library(reshape2)
V<-melt(Data,id= c('unit','time'))
Z<-V%>%
  mutate(variable=recode(variable,'pc.1'='PC1','pc.2'='PC2'))#,'pc.3'='DC3'
ggplot(data=Z,aes(x=time,y=value,color=factor(unit)))+
  geom_point(size=1)+geom_line()+
  scale_color_viridis(discrete = TRUE, option = "D")+
  theme_bw()+
  guides(color=guide_legend(title=''))+
  xlab('Million of cycles')+ylab('Fatigue crack length')+
  
  facet_wrap(variable~.,scales='free_y')
######################################MCEM算法收敛图像
load('F:\\degradation analysis\\Bivariate exponential model\\Code\\Crack.meekerpara.Rdata')
eta.s<-matrix(0,nrow=2000,ncol=2)
lambda.s<-matrix(0,nrow=2000,ncol=2)
Gam.s<-matrix(0,nrow=2000,ncol=2)
rho.s<-matrix(0,nrow=2000,ncol=2)
sigma.s<-matrix(0,nrow=2000,ncol=1)
AIC<-matrix(0,nrow=2000,ncol=1)
for(i in 1:2000){
  eta.s[i,]<-Para[[i]]$alpha.s
  lambda.s[i,]<-Para[[i]]$lambda.s
  Gam.s[i,]<-Para[[i]]$Gam.s
  rho.s[i,]<-Para[[i]]$rho.s
  sigma.s[i]<-Para[[i]]$sigma.s
  AIC[i]<-Para[[i]]$AIC
}
mean(AIC[1001:2000])
apply(Gam.s,2,mean)
para<-data.frame(eta1=eta.s[,1])
para$eta2<-eta.s[,2]
para$lambda1<-lambda.s[,1]
para$lambda2<-lambda.s[,2]
para$gamma1<-Gam.s[,1]
para$gamma2<-Gam.s[,2]
para$rho1<-rho.s[,1]
para$rho2<-rho.s[,2]
para$sigma<-sigma.s
para$AIC<-c(NA,AIC[2:2000])
para$step<-c(1:2000)
V<-melt(para,id= c('step'))
Z<-V%>%
  mutate(variable=recode(variable,
                         'eta1'='eta[1]','lambda1'='lambda[1]',
                         'eta2'='eta[2]','lambda2'='lambda[2]',
                         'gamma1'='gamma[1]','gamma2'='gamma[2]',
                         'rho1'='rho[1]','rho2'='rho[2]'
  ))#,'pc.3'='DC3'
ggplot(data=Z,aes(x=step,y=value))+
  geom_point(size=0.5)+
  theme_bw()+
  guides(color=guide_legend(title=''))+
  xlab('Iteration')+ylab('Value ')+
  
  facet_wrap(variable~.,scales='free_y',labeller = label_parsed)
