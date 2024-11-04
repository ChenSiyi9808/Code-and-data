rm(list = ls())

# Load required R packages
library(openxlsx)
library(ggplot2)
library(tweedie)
library(gridExtra)
library(dplyr)
library(viridis)
library(parallel)
library(doParallel)
library(reshape2)


#Load func_aS^b.R 

# Load multideg.3d.dat

# Filter data where rh is 0, temp is 25, and types is 1250
filtered_data1 <- subset(data, rh == 0 & temp == 25 & types == 1250)
# Filter data where rh is 0, temp is 25, and types is 1510
filtered_data2 <- subset(data, rh == 0 & temp == 25 & types == 1510)
# Filter data where rh is 0, temp is 25, and types is 2925
filtered_data3 <- subset(data, rh == 0 & temp == 25 & types == 2925)

# Generate a mapping from unit to id
id_to_unit <- unique(filtered_data1$id)
unit_map <- data.frame(id = id_to_unit, unit = 1:length(id_to_unit))

# Convert id to unit
data1 <- filtered_data1 %>%
  left_join(unit_map, by = "id") %>%
  select(unit, den, days, resp) %>%
  rename(time1 = days, pc.1 = resp)

# Generate a mapping from unit to id
id_to_unit <- unique(filtered_data2$id)
unit_map <- data.frame(id = id_to_unit, unit = 1:length(id_to_unit))

# Convert id to unit
data2 <- filtered_data2 %>%
  left_join(unit_map, by = "id") %>%
  select(unit, den, days, resp) %>%
  rename(time2 = days, pc.2 = resp)

# Generate a mapping from unit to id
id_to_unit <- unique(filtered_data3$id)
unit_map <- data.frame(id = id_to_unit, unit = 1:length(id_to_unit))

# Convert id to unit
data3 <- filtered_data3 %>%
  left_join(unit_map, by = "id") %>%
  select(unit, den, days, resp) %>%
  rename(time3 = days, pc.3 = resp)

################################################################################


# Combine data1 and data3 by columns
Data <- cbind(data1, data3[, c("time3", "pc.3")])
# Rename columns to change pc.3 and time3 to pc.2 and time2
colnames(Data)[5:6] <- c("time2", "pc.2")

Simu.deg  <- Data 
# Set initial values
n <- length(unique(Data$unit))
m <- 63
p <- 12 ###number of parameters
S <- c(10, 40, 60, 100)#####consider how to map to n one-to-one
num_S <- 4 
Psedo_failure <- matrix(0, nrow = n, ncol = 2)
Df <- c(0.5, 0.5) ####unknown
a.s <- matrix(0, nrow = n, ncol = 2)
b.s <- matrix(0, nrow = n, ncol = 2)
Gam.s <- matrix(0, nrow = n, ncol = 2)

## Data processing
for (i in 1:n) {
  # Extract data for the current unit
  x <- Simu.deg[which(Simu.deg$unit == i), ]
  
  # Check for consecutive identical values in the 4th column, and adjust if found
  for (j in 2:nrow(x)) {
    if (x[j, 4] == x[j - 1, 4]) {
      x[j - 1, 4] <- (x[j - 2, 4] + x[j, 4]) / 2
    }
  }
  
  # Check for consecutive identical values in the 6th column, and adjust if found
  for (j in 2:nrow(x)) {
    if (x[j, 6] == x[j - 1, 6]) {
      x[j - 1, 6] <- (x[j - 2, 6] + x[j, 6]) / 2
    }
  }
  
  # Calculate the difference
  x[2:nrow(x), c(4,6)] <- apply(x[1:nrow(x), c(4,6)], 2, diff)
  
  # Set the first row difference to NA
  x[1, c(4,6)] <- NA
  
  # Merge data
  if (i == 1) {
    Y <- x
  } else {
    Y <- rbind(Y, x)
  }
}

# Remove NA values and take absolute values
Y <- abs(na.omit(Y))
# Drop the time2 column
Y <- Y[, !names(Y) %in% "time2"]

# Rename the time1 column to time
names(Y)[names(Y) == "time1"] <- "time"

# View the result
#print(Y)


alpha.s<-c(6e-09,5e-09)#######the slope of pc2 should be smaller than pc1
sigma.s<-0.005
beta.s<-c(0.0008363282,0.0006363282)
rho.s<-c(2.25,2.29)
lambda.s<-c(0.01,0.01)
Gam.s<-c(0.9,0.8)
Delta<-0.1;L<-1

Para<-list()
#############################
#for(j in 1:2){
#     y[,(3+j)]<-pdf(y[,(3+j)],delta.t.s[,j],alpha.s[j]*exp(-beta.s[j]/y[,2]),lambda.s[j],rho.s[j])
# }
library(parallel)
library(doParallel)
# Initialize AIC
AIC <- 0

while(L<2001){
  Para[[L]]<-list(sigma.s=sigma.s,alpha.s=alpha.s,beta.s=beta.s,lambda.s=lambda.s,Gam.s=Gam.s,rho.s=rho.s,AIC=AIC)
  ##############E-step:MCMC algorithm
  init<-matrix(1,nrow=1,ncol=n)#rnorm(n,mean=1,sd=0.1)#
  par<-matrix(0,nrow=2000,ncol=n)
  
  for(b in 1:2000){
    init_hat<-init
    for(i in 1:n){
      y<-abs(Y[Y$unit==i,])
      t<-Simu.deg[Simu.deg$unit==i,3]
      delta.t.s<-apply(lambda.t(t,Gam.s),2,diff)
      init_hat[i]<-abs(rnorm(1,mean=init[i],sd=sigma.s))
      z<-runif(1,min=0,max=1)
      Alpha<-min(1,exp(h1(y,delta.t.s,init_hat[i],alpha.s,beta.s,lambda.s,rho.s,sigma.s)-h1(y,delta.t.s,init[i],alpha.s,beta.s,lambda.s,rho.s,sigma.s)))
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
  sigma.s<-sqrt(mean(apply((par[1001:2000,]-1)^2,2,mean)))
  ###########################
  S<-Y[,2]
  ###################
  lambda.f<-function(j){
    return(lambda.F(alpha.s[j],beta.s[j],rho.s[j],Gam.s[j],Y[,j+3],S,t,par))
  }
  lambda.s<-unlist(foreach(i=1:length(lambda.s))%do% lambda.f(i))
  ##########################
  rho.f<-function(j){
    fit<-function(rho){
      return(-para_M_step(alpha.s[j],beta.s[j],lambda.s[j],rho,Gam.s[j],Y[,j+3],S,t,par))
    }
    Fit<-optim(par=rho.s[j],fit,method='Brent',lower=0,upper=5)
    return(Fit$par)
  }
  rho.s<-unlist(foreach(j=1:length(rho.s))%do% rho.f(j))
  ####################################alpha
  alpha.f<-function(j){
    fit<-function(alpha){
      return(-para_M_step(alpha,beta.s[j],lambda.s[j],rho.s[j],Gam.s[j],Y[,j+3],S,t,par))
    }
    Fit<-optim(par=alpha.s[j],fit,method='Brent',lower=0,upper=1)#lower=1,1e+09)
    return(Fit$par)
  }
  alpha.s<-unlist(foreach(i=1:length(alpha.s))%do% alpha.f(i))
    ##########################beta
    beta.f<-function(j){
      fit<-function(beta){
        return(-para_M_step(alpha.s[j],beta,lambda.s[j],rho.s[j],Gam.s[j],Y[,j+3],S,t,par))
      }
      Fit<-optim(par=beta.s[j],fit,method='Brent',lower=0,upper=1)#10000)
      return(Fit$par)
    }
    beta.s<-unlist(foreach(i=1:length(beta.s))%do% beta.f(i))
    ###########################
    Gam.f<-function(j){
      fit<-function(Gam){
        return(-para_M_step(alpha.s[j],beta.s[j],lambda.s[j],rho.s[j],Gam,Y[,j+3],S,t,par))
      }
      Fit<-optim(par=Gam.s[j],fit,method='Brent',lower=0,upper=0.999)
      return(Fit$par)
    }
    Gam.s<-unlist(foreach(i=1:length(Gam.s))%do% Gam.f(i))
    ##########################################
    AIC <- -2 * aic(alpha.s, beta.s, lambda.s, rho.s, Gam.s, sigma.s, Y, delta.t.s, do.call(rbind, lapply(as.data.frame(par), function(x) x[1001:2000]))
                    , n, m) + 2 * 11 #101:200
    
    
    Delta<-max(c(abs(((Para[[L]]$alpha.s-alpha.s)/Para[[L]]$alpha.s)),
                 abs(Para[[L]]$lambda.s-lambda.s)/Para[[L]]$lambda.s,abs(Para[[L]]$Gam.s-Gam.s)/Para[[L]]$Gam.s))
    L<-L+1
  }
  
  alpha <- matrix(0, nrow=L, ncol=2)
  for (i in 1:(L-1)) {
    alpha[i, ] <- Para[[i]]$AIC
  }
  plot(alpha[2:(L-1), 1])
  
  
  ######################################
  N <- 2000
  a.s1 <- matrix(0, nrow=N, ncol=2)
  b.s1 <- matrix(0, nrow=N, ncol=2)
  lambda.s1 <- matrix(0, nrow=N, ncol=2)
  Gam.s1 <- matrix(0, nrow=N, ncol=2)
  rho.s1 <- matrix(0, nrow=N, ncol=2)
  sigma.s1 <- matrix(0, nrow=N, ncol=1)
  AIC1 <- matrix(0, nrow=N, ncol=1)
  

  for(i in 1:N){
    a.s1[i,] <- Para[[i]]$alpha.s
    b.s1[i,] <- Para[[i]]$beta.s
    lambda.s1[i,] <- Para[[i]]$lambda.s
    Gam.s1[i,] <- Para[[i]]$Gam.s
    rho.s1[i,] <- Para[[i]]$rho.s
    sigma.s1[i] <- Para[[i]]$sigma.s
    AIC1[i] <- Para[[i]]$AIC
  }
  
  

  para <- data.frame(
    alpha1 = c(NA, a.s1[, 1][2:N]),#a.s1[, 1],
    alpha2 = c(NA, a.s1[, 2][2:N]),#a.s1[, 2],
    beta1 = c(NA, b.s1[, 1][2:N]),#b.s1[, 1],
    beta2 = c(NA, b.s1[, 2][2:N]),#b.s1[, 2],
    lambda1 = lambda.s1[, 1],
    lambda2 = lambda.s1[, 2],
    gamma1 = Gam.s1[, 1],
    gamma2 = Gam.s1[, 2],
    rho1 = rho.s1[, 1],
    rho2 = rho.s1[, 2],
    sigma = sigma.s1,
    AIC = c(NA,NA, AIC1[3:N]),
    step = 1:N
  )
  

  V <- melt(para, id = c('step'))
  

  Z <- V %>%
    mutate(variable = recode(variable,
                             'alpha1' = 'alpha[1]', 'beta1' = 'beta[1]', 'lambda1' = 'lambda[1]',
                             'alpha2' = 'alpha[2]', 'beta2' = 'beta[2]', 'lambda2' = 'lambda[2]',
                             'gamma1' = 'gamma[1]', 'gamma2' = 'gamma[2]',
                             'rho1' = 'rho[1]', 'rho2' = 'rho[2]'
    ))
  
  

  ggplot(data = Z, aes(x = step, y = value)) +
    geom_point(size = 0.9) +
    theme_bw() +
    guides(color = guide_legend(title = '')) +
    xlab('Iteration') + ylab('Value ') +
    facet_wrap(variable ~ ., scales = 'free_y', labeller = label_parsed)
  
    