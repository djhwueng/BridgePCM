rm(list=ls())
library(ape)
library(sde)
library(geiger)
library(mvMORPH)
library(phangorn)
library(phytools)
library(phyclust)

DriftDiffusion<-function(model.params){
		mu.offset <<- model.params[1]
		alpha.offset <<- model.params[2]
		sigma.offset <<- model.params[3]
		drift <- expression(alpha.offset*(mu.offset-x))
		diffusion <- expression(sigma.offset)
		diffusion.x <- expression(0)

		D1<-function(t,x,theta){theta[2]*(theta[1]-x)}
		S1<-function(t,x,theta){theta[3]}
		return(list(expression.drift=drift,expression.diffusion=diffusion,expression.diffusion.x=diffusion.x,function.drift=D1,function.diffusion=S1))
		}

OUBridge <- function(model.params,start.value=start.value, end.value=end.value, t0=0, T=T, N=N, drift=drift,diffusion=diffusion,diffusion.x=diffusion.x){
  #mu<-model.params[1]
	alpha<-model.params[2]
	#sigma<-model.params[3]
	OUpath<-sde.sim(X0=start.value,drift=drift, sigma=diffusion, sigma.x=diffusion.x,N=N)
  t<-seq(t0,T,length=N+1)
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

sim.ou.path<-function(model.params, phy=phy, tip.states=tip.states, drift=drift, diffusion=diffusion, diffusion.x=diffusion.x, SE=NULL){
		ntips<-length(phy$tip.label)
	  edge.number<-dim(phy$edge)[1]
    edge.length<-phy$edge.length
		anc<- phy$edge[,1]
		des<- phy$edge[,2]
		anc.states<-condOU(model.params,phy=phy,tip.states=tip.states)$condmean
    full.node.data<-c(tip.states,anc.states)
		path.objects<- list()
		for(edgeIndex in edge.number:1){
        brnlen<-edge.length[edgeIndex]
        end.state<-full.node.data[des[edgeIndex]]
				start.state<-full.node.data[anc[edgeIndex]]
#        assign(paste("path",edgeIndex,sep=""), modified.DBridge(x = start.state, y = end.state, t0 = 0, N=ceiling(brnlen), drift=drift, sigma=sigma, sigma.x= 0, N.give.up=10))
        assign(paste("path",edgeIndex,sep=""), OUBridge(model.params, start.value=start.state, end.value=end.state, t0=0, T=1, N=ceiling(brnlen), drift=drift,diffusion=diffusion, diffusion.x=diffusion.x))
        path.objects<-c(path.objects, list(get(paste("path",edgeIndex,sep=""))))
			}#loop of edgeIndex
		return(path.objects)
		}#end of function


ounegloglike_onepath<-function(model.params,path=path){
			mu<-model.params[1]
			alpha<-model.params[2]
			sigma<-model.params[3]
	 		negloglike<- (length(path)+1)/2*log(pi*sigma^2)+path[1]^2/sigma^2
	 		for(pathIndex in 2:length(path)){
				#ref Bridge_NossanEvans2015.pdf formula (27), (38) for CIR but unknown contruction for bridge
				negloglike<-negloglike+ 1/2*log(1-exp(-2*alpha))+(path[pathIndex]-mu-(path[pathIndex-1]-mu)*exp(-alpha))^2/(sigma^2*(1-exp(-2*alpha)))
	 			}
	 		return(negloglike)
	 		}

ounegloglike<-function(model.params,phy=phy,path.data=path.data){
	mu<-model.params[1]
	alpha<-model.params[2]
	sigma<-model.params[3]
	#see Valdivieso.pdf for likelihood of ou
	 edge.number<-dim(phy$edge)[1]
	 negloglike.array<-array(0,c(edge.number))
	 for(edgeIndex in 1:edge.number){
	   path<-unlist(path.data[edgeIndex])
     negloglike.array[edgeIndex]<- ounegloglike_onepath(model.params,path=path)
     }
	 return(sum(negloglike.array))
	 }#end of loglike function

model.params<-c(10,1,2)
names(model.params)<-c("mu","alpha","sigma")
size<-10
min.time.step<-5
phy<-rcoal(size)
while(min(phy$edge.length)<min.time.step){
  phy$edge.length<-1.001*phy$edge.length
  }
print(phy$edge.length)
drift.diffusion<-DriftDiffusion(model.params)
expression.drift<-drift.diffusion$expression.drift
expression.diffusion<-drift.diffusion$expression.diffusion
expression.diffusion.x<-drift.diffusion$expression.diffusion.x
par(mfrow=c(1,2))
plot(phy)
tip.states<-rnorm(size,sd=1)

loop.size<-100
onepath.like.array<-array(0,c(loop.size))
for (i in seq(loop.size)) {
if(i %% 100 ==0){print(i)}
path.data<-sim.ou.path(model.params, phy=phy, tip.states=tip.states, SE=NULL, drift=expression.drift, diffusion=expression.diffusion, diffusion.x=expression.diffusion.x)
#print(ounegloglike_onepath(model.params,path=unlist(path.data[1])))
onepath.like.array[i] <- ounegloglike_onepath(model.params,path=unlist(path.data[1]))
}
print(mean(onepath.like.array))



source("~/GitHub/BridgePCM/plot_history.r")
plot.history(phy=phy,path.data=path.data,main="OUB")
print(ounegloglike(model.params,phy=phy,path.data=path.data))
true.tip.states<-tip.states
true.path.data<-path.data
optim(model.params,ounegloglike_onepath, method="L-BFGS-B", lower=c(-10,1e-8,1e-8),upper=Inf,path=unlist(path.data[1]))
optim(model.params,ounegloglike,method="L-BFGS-B",lower=c(-10,1e-8,1e-8),upper=Inf, phy=phy, path=path.data)
