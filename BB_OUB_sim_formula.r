rm(list=ls())
#library(sde)

# T<-1
# N=1000
# t0<-0
# dt<-1
# t<-seq(t0,T,length=N+1)
#
# start.value<-2
# end.value<-30
#
# sigma<-2
# X<-c(0,sigma*cumsum(rnorm(N)*sqrt(dt)))
# BB<-start.value+X-(t-t0)/(T-t0)*(X[N+1]-end.value + start.value)
# BB<-ts(BB,start=t0,deltat=dt)
#
# #?sde.sim
# alp<-1
# d<-expression(1*(20-x))
# s<-expression(2)
# X<-sde.sim(X0=start.value,drift=d,sigma=s,N=N)
# plot(X)
# Lambda <- exp(-alp*(T-t)) - exp(-alp*(T+t))
# Lambda <- Lambda/(1-exp(-2*alp*T))
# Lambda
# OUB <- X - Lambda*(X[N+1] - end.value)
# OUB <- as.numeric(OUB)
# #length(BB)
# #length(OUB)
# #par(mfrow=c(1,2))

sim.bb<-function (sigma, T = T, N = N, start.value=start.value, end.value=end.value){
#    if (T <= t0)
#    stop("wrong times")
    dt <- T/N #dt will be difference given difference brnlen, but overall we want 100 step on each brnlen
    t <- seq(0, T, length = N+1)
    path <- c(0, sigma*cumsum(rnorm(N) * sqrt(dt)))
    #Cite Archambeau
    BB <- start.value + path - t/T * (path[N + 1] - end.value + start.value)
    path <- ts(BB, start = 0, deltat = dt)
    return(invisible(path))
    }

sim.ou.one.path <- function(model.params, T=T, N=N, start.value=start.value){
		mu <- model.params[1]
		alpha <- model.params[2]
		sigma <- model.params[3]
		dw <- rnorm(N, 0, sqrt(T/N))
		dt <- T/N
		path<-c(start.value)
		for(Index in 2:(N+1)){
		  path[Index] <- path[Index-1] + alpha*(mu-path[Index-1])*dt + sigma*dw[Index-1]
		  }
		return(ts(path, start=start.value, deltat=dt))
		}

sim.ou.bridge <- function(model.params, T=T, N=N, start.value=start.value, end.value=end.value){
  #mu<-model.params[1]
	alpha<-model.params[2]
	#sigma<-model.params[3]
	#dt<-T/N
	#OUpath<-sde.sim(X0=start.value,drift=drift, sigma=diffusion, sigma.x=diffusion.x,N=N)
  OUpath<-sim.ou.one.path(model.params,T=T,N=N,start.value=start.value)
	t<-seq(0,T,length=N+1)
  #From Alison Etheridge. 1995. Stochastic partial differential equations. Cambridge University Press.
	Lambda <- exp(-alpha*(T-t)) - exp(-alpha*(T+t))
	Lambda <- Lambda/(1-exp(-2*alpha*T))
	OUB <- OUpath - Lambda*(OUpath[N+1] - end.value)
	return(as.numeric(OUB))
	}

mu<-10
alpha<-10
sigma<-4
model.params<-c(mu,alpha,sigma)
names(model.params)<-c("mu","alpha","sigma")

N=2000
T=1
start.value=2
end.value=10

BB<-sim.bb(sigma,T=T,N=N,start.value=start.value, end.value=end.value)
OUB<-sim.ou.bridge(model.params,T=T, N=N, start.value=start.value, end.value=end.value)


plot(NA,xlim=c(0,length(BB)), ylim=c(min(BB,OUB),max(BB,OUB)),ylab="Trait value",main="BB and OUB on a branch",xlab="Time step")
lines(1:(N+1),BB,type="l",lwd=2,col="blue")
lines(1:(N+1),OUB,type="l",lwd=2,col="red")
#lines(1:(N+1),as.numeric(X),type="l",col="purple",lwd=2)
#plot(BB,type="l",main="BB")
#plot(OUB,type="l",xlab="Time",main="OUB")
