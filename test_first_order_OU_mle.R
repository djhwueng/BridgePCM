rm(list=ls())
library(ape)
library(sde)
library(geiger)
library(mvMORPH)
library(phangorn)
library(phytools)
library(phyclust)

first.order.mu<-function(alpha,path=path){
  mu_hat <- (path[length(path)] - path[1]*exp(-alpha))/(1-exp(-alpha)) + sum(path[2:(length(path)-1)])
  mu_hat <- mu_hat*(1+exp(-alpha))/(length(path)*(1-exp(-alpha)))
  return(mu_hat)
  }
first.order.sigma.sq<-function(alpha,path=path){
  mu_hat <- first.order.mu(alpha,path=path)
  sigma.sq_hat<-0
  for(pathIndex in 2:length(path)){
    sigma.sq_hat <- sigma.sq_hat + (path[pathIndex]-mu_hat-(path[pathIndex-1]-mu_hat)*exp(-alpha))^2
    }
  sigma.sq_hat<-sigma.sq_hat*(2*alpha/(length(path)*(1-exp(-2*alpha))))
  return(sigma.sq_hat)
  }
first.order.ounegloglike.onepath<-function(alpha,path=path){
  mu_hat<-first.order.mu(alpha,path=path)
  sigma.sq_hat<-first.order.sigma.sq(alpha,path=path)
  negloglike<- length(path)/2*log(sigma.sq_hat/(2*alpha)) + length(path)/2*log(1-exp(-2*alpha))
  for(pathIndex in 2:length(path)){
    negloglike<-negloglike + alpha/(1-exp(-2*alpha))/sigma.sq_hat*(path[pathIndex]-mu_hat-(path[pathIndex-1] - mu_hat)*exp(-alpha))^2
    }
    return(negloglike)
  }



alpha<-50
d <-expression(50*(20-x))
s <-expression(3)
N<-100
sims<-100
mu.array<-array(0,c(sims))
sigma.array<-array(0,c(sims))
plot(sde.sim(X0=0,drift=d,sigma=s,sigma.x = 0,N=N))
for(Index in 1:sims){
  print(Index)
  path<-sde.sim(X0=0,drift=d,sigma=s,sigma.x = 0,N=N)
  mu.array[Index]<-first.order.mu(alpha,path=path)
  sigma.array[Index]<-sqrt(first.order.sigma.sq(alpha,path=path))
  }
print(mean(mu.array))
print(mean(sigma.array))
#parameter estimation in OU process is issue. We 

sims<-10
alp.array<-array(0,c(sims))
mu.array<-array(0,c(sims))
sigma.sq.array<-array(0,c(sims))
for(simIndex in 1:sims){
  print(simIndex)
  path<-sde.sim(X0=0,drift=d,sigma=s,sigma.x = 0,N=100)
  temp.result<-optimize(first.order.ounegloglike.onepath,c(0,20),path=path)
  alp.array[simIndex]<-temp.result$minimum
  mu.array[simIndex]<-first.order.mu(alp.array[simIndex],path=path)
  sigma.sq.array[simIndex]<-first.order.sigma.sq(alp.array[simIndex],path=path)
  }

output<-cbind(alp.array,mu.array,sigma.sq.array)
colnames(output)<-c("alp","mu","sig.sq")
head(output)
apply(output,2,mean)
#plot(path)
#alpha <- 4; mu<- 2; sigma <-0.5
#print(first.order.ounegloglike.onepath(alpha,path=path))
