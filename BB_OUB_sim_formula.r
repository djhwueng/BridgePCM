rm(list=ls())
library(sde)

T<-1
N=100
t0<-0
dt<-1
t<-seq(t0,T,length=N+1)

start.value<-2
end.value<-30

sigma<-2
X<-c(0,sigma*cumsum(rnorm(N)*sqrt(dt)))
BB<-start.value+X-(t-t0)/(T-t0)*(X[N+1]-end.value + start.value)
BB<-ts(BB,start=t0,deltat=dt)

?sde.sim
alp<-0.0001
d<-expression(-0.0001*x)
s<-expression(20)
X<-sde.sim(X0=start.value,drift=d,sigma=s,N=N)
#plot(X)
Lambda <- exp(-alp*(T-t)) - exp(-alp*(T+t))
Lambda <- Lambda/(1-exp(-2*alp*T))
Lambda
OUB <- X - Lambda*(X[N+1] - end.value) 
OUB <- as.numeric(OUB)
#length(BB)
#length(OUB)
#par(mfrow=c(1,2))
plot(NA,xlim=c(0,length(BB)), ylim=c(min(BB,OUB,as.numeric(X)),max(BB,OUB,as.numeric(X))),ylab="trait value",main="BB and OUB")
lines(1:(N+1),BB,type="l",col="red")
lines(1:(N+1),OUB,type="l",col="blue")
lines(1:(N+1),as.numeric(X),type="l",col="purple")

#plot(BB,type="l",main="BB")
#plot(OUB,type="l",xlab="Time",main="OUB")