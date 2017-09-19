rm(list=ls())
library(ape)
library(sde)
library(geiger)
library(mvMORPH)
library(phangorn)
library(phytools)
library(phyclust)

first.order.mu<-function(alpha,path=path, T=T, N=N){
  dt<-T/N
  A<-0
  mu.hat<-0
  for(Index in 2:(N+1)) {
    A <-  A  + (1-exp(-alpha*dt))/(1+exp(-alpha*dt))
    mu.hat <- mu.hat + (path[Index] - path[Index-1]*exp(-alpha*dt))/(1+exp(-alpha*dt))
    }
  A <- 1/A
  return(mu.hat*A)
  }

first.order.sigma.sq <- function(alpha, path=path, T=T, N=N){
  dt<-T/N
  mu.hat <- first.order.mu(alpha, path=path, T=T, N=N)
  sigma.sq.hat<-0
  for(Index in 2:(N+1)){
    sigma.sq.hat <- sigma.sq.hat + (path[Index] - mu.hat - ( path[Index-1] - mu.hat )*exp(-alpha*dt))^2/(1-exp(-2*alpha*dt))
    }
  sigma.sq.hat<-sigma.sq.hat*2*alpha/N
  return(sigma.sq.hat)
  }

first.order.ounegloglike.onepath<-function(alpha,path=path,T=T,N=N){
  dt<-T/N
  mu.hat<-first.order.mu(alpha,path=path,T=T,N=N)
  sigma.sq.hat<-first.order.sigma.sq(alpha,path=path,T=T,N=N)
  negloglike<- N/2*log(sigma.sq.hat/(2*alpha))
  for(Index in 2:(N+1)) {
    negloglike <- negloglike + 1/2*log(1-exp(-2*alpha*dt))
    negloglike <- negloglike + alpha/sigma.sq.hat*(path[Index] - mu.hat - (path[Index-1] - mu.hat)*exp(-alpha*dt))^2/(1-exp(-2*alpha*dt))
    }
    return(negloglike)
  }

ou.path.sim <- function(model.params, T=T, N=N, x0=x0){
  mu <- model.params[1]
  alpha <- model.params[2]
  sigma <- model.params[3]
  dw <- rnorm(N, 0, sqrt(T/N))
  dt <- T/N
  path<-c(x0)
  for(Index in 2:(N+1)){
    path[Index] <- path[Index-1] + alpha*(mu-path[Index-1])*dt + sigma*dw[Index-1]
    }
  return(ts(path, start=x0, deltat=dt))
  }

T<-1
N<-1000
x0<-0
sims<-1000
true.alpha<-1.2
true.mu<-16
true.sigma<-4
true.params<-c(true.mu, true.alpha, true.sigma)

mu.array.true.alpha<-array(0,c(sims))
sigma.array.true.alpha<-array(0,c(sims))

mu.array.mle.alpha<-array(0,c(sims))
alpha.array.mle.alpha<-array(0,c(sims))
sigma.array.mle.alpha<-array(0,c(sims))

for(Index in 1:sims){
  print(Index)
  path<-ou.path.sim(true.params, T=T, N=N ,x0=x0)
  mu.array.true.alpha[Index]<-first.order.mu(true.alpha,path=path, T=T, N=N)
  sigma.array.true.alpha[Index]<-sqrt(first.order.sigma.sq(true.alpha,path=path, T=T, N=N))
  alpha.array.mle.alpha[Index]<-optimize(first.order.ounegloglike.onepath,c(1,10),path=path,T=T,N=N)$minimum
  mu.array.mle.alpha[Index]<-first.order.mu(alpha.array.mle.alpha[Index], path=path, T=T, N=N)
  sigma.array.mle.alpha[Index]<-sqrt(first.order.sigma.sq(alpha.array.mle.alpha[Index], path=path, T=T, N=N))
  }
print(mean(mu.array))
print(mean(sigma.array))


output<-cbind(mu.array.true.alpha,sigma.array.true.alpha,alpha.array.mle.alpha,mu.array.mle.alpha,sigma.array.mle.alpha)
colnames(output)<-c("mu.true.alpha","sigma.true.alpha","alp.mle","mu.mle.alpha","sigma.mle.alpha")
head(output)
apply(output,2,mean)
#plot(path)
#alpha <- 4; mu<- 2; sigma <-0.5
#print(first.order.ounegloglike.onepath(alpha,path=path))
