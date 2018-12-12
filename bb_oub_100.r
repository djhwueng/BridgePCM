setwd("/Users/djhwueng/Dropbox/CollabJhwueng_pcmBridge/Rcode/")
rm(list=ls())
library(sde)

############################BB
T<-1
N=1000
t0<-0
dt<-1
t<-seq(t0,T,length=N+1)

start.value=2
end.value=10

sigma<-2

sim<-20
BBarray<-array(0,c(sim,N+1))

for (simIndex in 1:sim){
  X<-c(0,sigma*cumsum(rnorm(N)*sqrt(dt)))
  BB<-start.value+X-(t-t0)/(T-t0)*(X[N+1]-end.value+start.value)
  BB<-ts(BB,start=t0,deltat = dt)
  BBarray[simIndex,]<-BB
}

BBarray[1,]



####################OUB

T<-1
N=1000
t0<-0
dt<-1
t<-seq(t0,T,length=N+1)

start.value=2
end.value=10

alp<-1
d<-expression(-1*x)
s<-expression(20)

sims<-20
OUBarray<-array(0,c(sims,N+1))
for(simIndex in 1:sims){
  X<-sde.sim(X0=start.value,drift=d,sigma=s,N=N)
  Lambda<-exp(-alp*(T-t))-exp(-alp*(T+t))
  Lambda<-Lambda/(1-exp(-2*alp*T))
  OUB<-X-Lambda*(X[N+1]-end.value)
  OUB<-as.numeric(OUB)
  OUBarray[simIndex,]<-OUB
}

####################PLOTS

png("bb100.png")
plot(NA,xlim=c(0,length(BB)), ylim=c(min(BBarray),max(BBarray)),xlab="time step",ylab="trait value", main="Brownian Bridge")
for(simIndex in 1:sims){
  lines(1:(N+1),BBarray[simIndex,],type = "l")
  }
dev.off()

png("oub100.png")
plot(NA,xlim=c(0,length(OUB)), ylim=c(min(BBarray),max(BBarray)),xlab="time step",ylab="trait value", main="Ornstein-Uhlenbeck Bridge")
for(simIndex in 1:sims){
  lines(1:(N+1),OUBarray[simIndex,],type = "l")
  }
dev.off()

