

#######################d function for ED distribution
d<-function(y,t,mu,rho){
  #y:degradation measure
  #t:measure interval
  #mu:mean
  #rho:model coefficient rho=0 for wiener process
  if(rho==0){
    z<-(y/t-mu)^2
  }else {
    if(rho==1){
      z<-2*(y/t*log(y/(t*mu))-(y/t-mu))
    } else{
      if(rho==2){
        z<-2*(log((mu*t)/y)+y/(mu*t)-1)
      }else{
        z<-2*(((y/t)^(2-rho))/((1-rho)*(2-rho))-(mu^(1-rho)*y)/((1-rho)*t)+(mu^(2-rho))/(2-rho))
      }
    }
  }
  return(z)
}

#########################################power time transformation
lambda.t<-function(t,Gam){
  tau<-matrix(NA,ncol=length(Gam),nrow=length(t))
  for(j in 1:length(Gam)){
    tau[,j]<-t^Gam[j] 
  }
  
  return(tau)
}

#############################


pdf<-function(y,t,mu,lambda,rho){
  #y:degradation measure
  #t:measure interval
  #mu:mean
  #lambda: dispersion
  #rho:model coefficient rho=0 for wiener process
  #likelihoodtype:lik>likelihood;loglik>log-likelihood
  d<-function(y,t,mu,rho){
    #y:degradation measure
    #t:measure interval
    #mu:mean
    #rho:model coefficient rho=0 for wiener process
    if(rho==0){
      z<-(y/t-mu)^2
    }else {
      if(rho==1){
        z<-2*(y/t*log(y/(t*mu))-(y/t-mu))
      } else{
        if(rho==2){
          z<-2*(log((mu*t)/y)+y/(mu*t)-1)
        }else{
          z<-c(2*(((y/t)^(2-rho))/((1-rho)*(2-rho))-(mu^(1-rho)*y)/((1-rho)*t)+(mu^(2-rho))/(2-rho)))
        }
      }
    }
    return(z)
  }
  z<-0
  z<-0.5*log(lambda/(2*pi*t^(1-rho)*y^rho))-(lambda*t)/2*d(y,t,mu,rho)
  return(z)
}


h1<-function(y,delta.t.s,mu,alpha,beta,lambda,rho,sigma){
  #y:degradation measure
  #t:measure interval
  #mu:mean
  #lambda: dispersion
  #rho:model coefficient rho=0 for wiener process
  for(j in 1:2){
    y[,(3+j)]<-pdf(y[,(3+j)],delta.t.s[,j],mu*alpha[j]*(y[,2])^beta[j],lambda[j],rho[j])
  }
  y1<-log(dnorm(mu,mean=1,sd=sigma))
  z<-0
  z<-(sum(y[,4:5])+sum(y1))
  return(z)
}

##################################likelihood for MLE
para_M_step<-function(alpha,beta,lambda,rho,Gam,y,S,t,par){
  n<-ncol(par)
  tau<-(apply(lambda.t(t,Gam),2,diff))
  m<-length(tau)
  y<-matrix(y,nrow=n,ncol=m,byrow=T)
  S<-matrix(S,nrow=n,ncol=m,byrow=T)
  Q<-matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    x<-0
    for(k in 1:m){
      x<-d(y[i,k],tau[k],par[,i]*alpha*(S[i,k])^beta,rho)
      Q[i,k]<--0.5*(1-rho)*log(tau[k])-0.5*rho*log(y[i,k])-0.5*(lambda*tau[k])*(mean(x[101:200]))#0.5*log(lambda/(2*pi))
    }
  }
  z<-sum(Q)
  return(z)
}

#################################
########################################

lambda.F<-function(alpha,beta,rho,Gam,y,S,t,par){
  n<-ncol(par)
  tau<-(apply(lambda.t(t,Gam),2,diff))
  m<-length(tau)
  y<-matrix(y,nrow=n,ncol=m,byrow=T)
  S<-matrix(S,nrow=n,ncol=m,byrow=T)
  Q<-matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    x<-0
    for(k in 1:m){
      x<-d(y[i,k],tau[k],alpha*(S[i,k])^beta*par[,i],rho)
      Q[i,k]<-tau[k]*(mean(x[101:200]))
    }
  }
  z<-(n*m)/sum(Q)
  return(z)
}

##################################AIC value 
aic<-function(alpha.s,beta.s, lambda.s,rho.s,Gam.s,sigma.s,Y,delta.t.s,par,n,m){
  P<-matrix(NA,nrow=n,ncol=1)
  
  for(i in 1:n){
    x<-Y[which(Y$unit==i),4:5]
    for(j in 1:2){
      Q<-0
      for(k in 1:m-1){
        Q<-d(x[k,j],delta.t.s[k,j],par[,i]*alpha.s[j]*(y[,2])^beta.s[j],rho.s[j])
        x[k,j]<-0.5*log(lambda.s[j]/(2*pi))-(1-rho.s[j])/2*log(delta.t.s[k,j])-rho.s[j]/2*log(x[k,j])-(lambda.s[j]*delta.t.s[k,j])/2*(mean(Q))
      }
    }
    Y[which(Y$unit==i),4:5]<-x
    P[i,1]<-mean(log(dnorm(par[,i],mean=1,sd=sigma.s)))
  }
  
  return(sum(Y[,4:5])+sum(P))
}
############################################ NORMAL CONDITION OF UV: 5%

Reliability<-function(alpha.s,beta.s,lambda.s,rho.s,Gam.s,sigma.s,t){
  Df<-c(0.4,0.27)
  if(rho.s[1]==2){
    par<-1/rgamma(100000,shape=b.s,scale=1/b.s)
  }else if(rho.s[1]==3){
    par<-1/rnorm(100000,mean=1,sd=sigma.s)
  }else{
    par<-rnorm(100000,mean=1,sd=sigma.s)
  }
  
  delta.t.s<-lambda.t(t,Gam.s)
  x1<-sqrt(lambda.s[1]/(par^(rho.s[1])*(alpha.s[1]*(5)^beta.s[1])^rho.s[1]))*(par*(alpha.s[1]*(5)^beta.s[1])*sqrt(delta.t.s[1])-Df[1]/(sqrt(delta.t.s[1])))
  x2<-sqrt(lambda.s[2]/(par^(rho.s[2])*(alpha.s[2]*(5)^beta.s[2])^rho.s[2]))*(par*(alpha.s[2]*(5)^beta.s[2])*sqrt(delta.t.s[2])-Df[2]/(sqrt(delta.t.s[2])))
  F1<-mean(pnorm(x1))
  F2<-mean(pnorm(x2))
  F12<-mean(pnorm(x1)*pnorm(x2))
  R<-1-F1-F2+F12
  return(R)
}
