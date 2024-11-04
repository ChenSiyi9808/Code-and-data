BootstrapRandomSample<-function(b,Para){
  alpha<-Para[[b]]$alpha.s
  Gam<-Para[[b]]$Gam.s
  lambda<-1/Para[[b]]$lambda.s
  rho<-Para[[b]]$rho.s
  sigma<-Para[[b]]$sigma.s
  library(tweedie)
  S<-1
  Simulation<-list()
  while(S<61){
    n<-10#######sample number
    m<-10;#measurement numbber
    eta<-1
    mu<-matrix(0,nrow=n,ncol=length(eta))
    for(i in 1:n){
      mu[i,]<-rnorm(1,mean=eta,sd=sigma)
    }
    t<-c(0:m)*(10/m)
    #############################time transformation function
    lambda.t<-function(t,Gam){
      tau<-matrix(NA,ncol=length(Gam),nrow=length(t))
      for(j in 1:length(Gam)){
        tau[,j]<-t^Gam[j] 
      }
      
      return(tau)
    }
    tau<-lambda.t(t,Gam)
    delta.tau<-apply(tau,2,diff)
    #generate the degradation path
    for(i in 1:n){
      x<-matrix(0,nrow=m,ncol=length(Gam))
      for(j in 1:length(Gam)){
        for(k in 1:m){
          x[k,j]<-rtweedie(n=1, power=rho[j], mu=alpha[j]*mu[i]*delta.tau[k,j], phi=lambda[j]*((delta.tau[k,j])^(1-rho[j])))
        }
      }
      x<-rbind(c(0,0),x)
      if(i==1){
        Simu.deg<-data.frame(unit=i,time=t,pc=apply(x,2,cumsum))
      }else{
        simu.deg<-data.frame(unit=i,time=t,pc=apply(x,2,cumsum))
        Simu.deg<-rbind(Simu.deg,simu.deg)
      }
    }
    Simulation[[S]]<-Simu.deg
    S<-S+1
    
  }
  return(Simulation)
}  
#EM for Bootstrap method
f<-function(S,Simulation){
  Simu.deg<-Simulation[[S]]
  source(file='/home/zhanghu1/BEDM/func.R')
  library(tweedie)
  lambda.t<-function(t,Gam){
    tau<-matrix(NA,ncol=length(Gam),nrow=length(t))
    for(j in 1:length(Gam)){
      tau[,j]<-t^Gam[j] 
    }
    
    return(tau)
  }
  
  Simu.deg<-Simulation[[S]]
  n<-10;m<-10
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
  sigma.s<-0.1
  rho.s<-c(2.7,2.2)
  lambda.s<-c(0.1,0.1)
  Gam.s<-c(1,1.5)
  Delta<-0.1;L<-1
  Para<-list()
  library(parallel)
  library(doParallel)
  #time_start<-Sys.time()
  AIC<-0
  while(L<201){
    Para[[L]]<-list(sigma.s=sigma.s,alpha.s=alpha.s,lambda.s=lambda.s,Gam.s=Gam.s,rho.s=rho.s,AIC=AIC)
    ##############E-step:MCMC algorithm
    init<-matrix(1,nrow=1,ncol=n)#rnorm(n,mean=1,sd=0.05)#
    par<-matrix(0,nrow=200,ncol=n)
    delta.t.s<-apply(lambda.t(t,Gam.s),2,diff)
    for(b in 1:200){
      init_hat<-init
      for(i in 1:n){
        y<-Y[Y$unit==i,]
        init_hat[i]<-abs(rnorm(1,mean=init[i],sd=0.05))
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
    #plot(par[,2])
    #sigma.s<-sqrt(mean(apply((par[101:200,]-1)^2,2,mean)))
    #############M-step
    sigma.s<-sqrt(mean(apply((par[101:200,]-1)^2,2,mean)))
    lambda.f<-function(j){
      return(lambda.F(alpha.s[j],rho.s[j],Gam.s[j],Y[,j+2],t,par))
    }
    
    lambda.s<-unlist(foreach(i=1:length(lambda.s))%do% lambda.f(i))
    
    
    ####################################alpha
    alpha.f<-function(j){
      fit<-function(alpha){
        return(-para_M_step(alpha,lambda.s[j],rho.s[j],Gam.s[j],Y[,j+2],t,par))
      }
      Fit<-optim(par=alpha.s[j],fit,method='Brent',lower=0,upper=20)
      return(Fit$par)
    }
    alpha.s<-unlist(foreach(i=1:length(alpha.s))%do% alpha.f(i))
    ##########################
    rho.f<-function(j){
      fit<-function(rho){
        return(-para_M_step(alpha.s[j],lambda.s[j],rho,Gam.s[j],Y[,j+2],t,par))
      }
      Fit<-optim(par=rho.s[j],fit,method='Brent',lower=0,upper=5)
      return(Fit$par)
    }
    rho.s<-unlist(foreach(i=1:length(rho.s))%do% rho.f(i))
    ###########################
    Gam.f<-function(j){
      fit<-function(Gam){
        return(-para_M_step(alpha.s[j],lambda.s[j],rho.s[j],Gam,Y[,j+2],t,par))
      }
      Fit<-optim(par=Gam.s[j],fit,method='Brent',lower=0,upper=2)
      return(Fit$par)
    }
    Gam.s<-unlist(foreach(i=1:length(Gam.s))%do% Gam.f(i))
    AIC<--2*aic(alpha.s,lambda.s,rho.s,Gam.s,sigma.s,Y,delta.t.s,par[101:200,],n,m)+2*9
    Delta<-max(c(abs(((Para[[L]]$alpha.s-alpha.s)/Para[[L]]$alpha.s)),
                 abs(Para[[L]]$lambda.s-lambda.s)/Para[[L]]$lambda.s,abs(Para[[L]]$Gam.s-Gam.s)/Para[[L]]$Gam.s))
    L<-L+1
  }
  eta.s1<-matrix(0,nrow=200,ncol=2)
  lambda.s1<-matrix(0,nrow=200,ncol=2)
  Gam.s1<-matrix(0,nrow=200,ncol=2)
  rho.s1<-matrix(0,nrow=200,ncol=2)
  sigma.s1<-matrix(0,nrow=200,ncol=1)
  AIC1<-matrix(0,nrow=200,ncol=1)
  for(i in 1:200){
    eta.s1[i,]<-Para[[i]]$alpha.s
    lambda.s1[i,]<-Para[[i]]$lambda.s
    Gam.s1[i,]<-Para[[i]]$Gam.s
    rho.s1[i,]<-Para[[i]]$rho.s
    sigma.s1[i]<-Para[[i]]$sigma.s
    AIC1[i]<-Para[[i]]$AIC
  }
  par<-list(sigma.s=mean(sigma.s1[101:200]),
            alpha.s=apply(eta.s1[101:200,],2,mean),
            lambda.s=apply(lambda.s1[101:200,],2,mean),
            Gam.s=apply(Gam.s1[101:200,],2,mean),
            rho.s=apply(rho.s1[101:200,],2,mean),
            AIC=mean(AIC1[101:200]))
  
  return(par)
}



load('/home/zhanghu1/BEDM/Simulationpara1010.Rdata')
hat.parameter2<-list()
library(doParallel)
cl <- makeCluster(60)
registerDoParallel(cl)
time_start<-Sys.time()
hat.parameter2<- foreach(h=1:100) %dopar% Boot(h,Para,BootstrapRandomSample)
stopCluster(cl)
exc_time<-difftime(Sys.time(),time_start,units = 'mins')
print(paste0('code执行时间：',round(exc_time,2),'mins'))
save(hat.parameter2,'/home/zhanghu1/BEDM/hat.parameterBoot.Rdata')
Boot<-function(b,Para,BootstrapRandomSample){
  ###########################################
  Simulation<-BootstrapRandomSample(b,Para)
  #############################EM算法求解上述模拟数据集
  lambda.t<-function(t,Gam){
    tau<-matrix(NA,ncol=length(Gam),nrow=length(t))
    for(j in 1:length(Gam)){
      tau[,j]<-t^Gam[j] 
    }
    
    return(tau)
  }
  library(parallel)
  source('/home/zhanghu1/BEDM/EMBootstrap.R')
  # 启用parallel作为foreach并行计算的后端
  initial<-Para[[b]]
  hat.parameter1<-list()
  library(doParallel)
  cl <- makeCluster(60)
  registerDoParallel(cl)
  hat.parameter1<- foreach(r=1:60) %dopar% f(r,Simulation)
  stopCluster(cl)
  ##########################################################
  L<-60
  eta.s1<-matrix(0,nrow=L,ncol=2)
  lambda.s1<-matrix(0,nrow=L,ncol=2)
  Gam.s1<-matrix(0,nrow=L,ncol=2)
  rho.s1<-matrix(0,nrow=L,ncol=2)
  sigma.s1<-matrix(0,nrow=L,ncol=1)
  for(i in 1:L){
    eta.s1[i,]<-hat.parameter1[[i]]$alpha.s
    lambda.s1[i,]<-hat.parameter1[[i]]$lambda.s
    Gam.s1[i,]<-hat.parameter1[[i]]$Gam.s
    rho.s1[i,]<-hat.parameter1[[i]]$rho.s
    sigma.s1[i]<-hat.parameter1[[i]]$sigma.s
  }
  parameter<-as.data.frame(cbind(eta.s1,Gam.s1,lambda.s1,rho.s1))
  parameter$cycle<-b
  return(parameter)
}

