###################generating bivariate ED process
for(S in 1:1000){
  Simu.deg<-Simulation[[S]]
  for(i in 1:n){
    x<-Simu.deg[which(Simu.deg$unit==i),]
    y<-Simu.deg[which(Simu.deg$unit==i),]
    x[2:nrow(x),3:(3+1)]<-apply(x[1:nrow(x),3:(3+1)],2,diff)
    x[1,3:(3+1)]<-NA
    x<-na.omit(x)
    if(length(which(x[,3:4]<1e-4))>0){#nrow(x[,3:4])<50
      print(S)
      Simu.deg<-Simulation[[S+1]]
    }
    
  }
  Simulation[[S]]<-Simu.deg
}

y[which(x[,4]<1e-6)+1,4]<-0.5*(y[which(x[,4]<1e-6)+2,4]+y[which(x[,4]<1e-6),4])
y[which(x[,3]<1e-6)+1,3]<-0.5*(y[which(x[,3]<1e-6)+2,3]+y[which(x[,3]<1e-6),3])
Simu.deg[which(Simu.deg$unit==i),]<-y
if(length(which(x[,3:4]==0))>0){
  y[which(x[,4]==0)+1,4]<-0.5*(y[which(x[,4]==0)+2,4]+y[which(x[,4]==0),4])
  y[which(x[,3]==0)+1,3]<-0.5*(y[which(x[,3]==0)+2,3]+y[which(x[,3]==0),3])
  Simu.deg[which(Simu.deg$unit==i),]<-y
}
library(tweedie)
library(ggplot2)
Simulation<-list()
S<-1
while(S<1001){
  n<-10#######sample number
  m<-50;#measurement numbber
  eta<-1
  alpha<-c(1.5,1)
  sigma<-0.1
  rho<-c(2.7,2.2)
  lambda<-c(0.1,0.1)
  Gam<-c(1,1.5)
  mu<-matrix(0,nrow=n,ncol=length(eta))
  for(i in 1:n){
    mu[i,]<-rnorm(1,mean=eta,sd=sigma)
  }
  t<-c(0:m)*(10/m)##c(0.00 ,0.02, 0.03, 0.04, 0.05 ,0.06, 0.07, 0.08 ,0.09)
    #
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
  #main='F:\\degradation analysis\\Bivariate exponential model\\Code\\'
  #Link<-paste(main,paste(S,'Rdata',sep='.'))
  #save(Simu.deg,file=Link)
  print(S)
  #if(any(Simu.deg<0)==TRUE){}else{
  Simulation[[S]]<-Simu.deg
  S<-S+1
  #}
}
Simulation[[S]]<-Simu.deg
ggplot(data=Simu.deg,aes(x=time,y=pc.1,fill=factor(unit)))+
  geom_point()+geom_line()+
  geom_line(aes(x=time,y=pc.2,color="red"))+
  #geom_line(aes(x=time,y=pc.3,color="blue"))+
  scale_fill_discrete(guide=FALSE)+guides(color=guide_legend(title=NULL))

load(file='F:\\degradation analysis\\Bivariate exponential model\\Code\\Simulationmeeker.Rdata')
for(S in 501:1000){
  main='F:\\degradation analysis\\Bivariate exponential model\\Code\\'
  Link<-paste(main,paste(S,'Rdata',sep='.'))
  load(Link)
  Simulation[[S]]<-Simu.deg
}

save(Simulation,file='F:\\degradation analysis\\Bivariate exponential model\\Code\\Simulationmeeker.Rdata')
source(file='/home/zhanghu1/BEDM/randomsample.R')
  Simulation<-list()
  library(doParallel)
  cl <- makeCluster(60)
  registerDoParallel(cl)
  time_start<-Sys.time()
  Simulation<- foreach(r=1:1000) %dopar% Randsample(r)
  stopCluster(cl)
  exc_time<-difftime(Sys.time(),time_start,units = 'mins')
  print(paste0('code执行时间：',round(exc_time,2),'mins'))
  save(Simulation, file='/home/zhanghu1/BEDM/Simulation1050.Rdata')
 