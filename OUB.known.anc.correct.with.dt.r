rm(list=ls())
library(ape)
library(sde)
library(geiger)
library(mvMORPH)
library(phangorn)
library(phytools)
library(phyclust)

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

sim.ou.bridge <- function(model.params,start.value=start.value, end.value=end.value, T=T, N=N){
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

condOU <- function(model.params, phy=phy, tip.states=tip.states) {
	mu<-model.params[1]
  alpha<-model.params[2]
  sigma<-model.params[3]
	ntip <- length(phy$tip.label)
	N <- dim(phy$edge)[1] + 1
	covmatrix = sigma^2/2/alpha*exp(-alpha*dist.nodes(phy))

	cov22 = covmatrix[1:ntip, 1:ntip]
	cov11 = covmatrix[(ntip+1):N, (ntip+1):N]
	cov12 = covmatrix[1:ntip, (ntip+1):N]

	inv_cov22 = solve(cov22)

	condmean = mu + t(cov12)%*%inv_cov22%*%(tip.states - rep(mu,ntip))
	condvar =  cov11 - t(cov12)%*%inv_cov22%*%cov12
	return(list(condmean = condmean, condvar = condvar))
	}

sim.ou.tree.path<-function(model.params, phy=phy, tip.states=tip.states, N=N, SE=NULL){
	ntips<-length(phy$tip.label)
	edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
	anc<- phy$edge[,1]
	des<- phy$edge[,2]
	anc.states<-condOU(model.params,phy=phy,tip.states=tip.states)$condmean
  full.node.data<-c(tip.states,anc.states)
	path.objects<- list()
	for(edgeIndex in edge.number:1){
      brlen<-edge.length[edgeIndex]
      end.value<-full.node.data[des[edgeIndex]]
			start.value<-full.node.data[anc[edgeIndex]]
#        assign(paste("path",edgeIndex,sep=""), modified.DBridge(x = start.state, y = end.state, t0 = 0, N=ceiling(brnlen), drift=drift, sigma=sigma, sigma.x= 0, N.give.up=10))
      assign(paste("path",edgeIndex,sep=""), sim.ou.bridge(model.params, start.value=start.value, end.value=end.value, T=brlen, N=N))
      path.objects<-c(path.objects, list(get(paste("path",edgeIndex,sep=""))))
			}#loop of edgeIndex
		return(path.objects)
		}#end of function

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
#so far onepath ounegloglike optimization works well using first order method.
#when move to tree case, we shall be able to use first order method.
first.order.ounegloglike.tree<-function(alpha,phy=phy,path.data=path.data, N=N){
	  badval<-(0.5)*.Machine$double.xmax
 		edge.number<-dim(phy$edge)[1]
    edge.length<-phy$edge.length
		negloglike.array<-array(0,c(edge.number))
 		for(edgeIndex in edge.number:1){
			brlen<-edge.length[edgeIndex]
   		path<-unlist(path.data[edgeIndex])
   		negloglike.array[edgeIndex]<- first.order.ounegloglike.onepath(alpha,path=path,T=brlen,N=N)
   		}
 		return(sum(negloglike.array))
 		}#end of loglike function

N=5
model.params<-c(10,0.1,4)
names(model.params)<-c("mu","alpha","sigma")
tree.size<-3
min.time.step<-5
phy<-rcoal(tree.size)
#while(min(phy$edge.length)<min.time.step){
#  phy$edge.length<-1.001*phy$edge.length
#  }
#print(phy$edge.length)
tip.states<-rnorm(tree.size,sd=1)
path.data<-sim.ou.tree.path(model.params, phy=phy, tip.states=tip.states, N=N, SE=NULL)
#source("~/GitHub/BridgePCM/plot_history.r")
#par(mfrow=c(1,2))
plot(phy)
#plot.history(phy=phy,path.data=path.data,main="OUB")
print(first.order.ounegloglike.tree(0.1,phy=phy,path.data=path.data,N=N))
path.data
phy$edge.length

?optimize
optimize(first.order.ounegloglike.tree,c(0,1000000),phy=phy,path.data=path.data,N=N)
#no good for alphas, we do not have 

edge.number<-dim(phy$edge)[1]
edge.length<-phy$edge.length
negloglike.array<-array(0,c(edge.number))
alpha<-0.1
for(edgeIndex in edge.number:1){
  path<-unlist(path.data[edgeIndex])
  brlen<-edge.length[edgeIndex]
  mu.hat<-first.order.mu(alpha,path=path,T=brlen,N=N)
  sigma.sq.hat<-first.order.sigma.sq(alpha,path=path,T=brlen,N=N)
  print(c(mu.hat,sigma.sq.hat))
}
#single path case estimate well, but tree case not.......yet