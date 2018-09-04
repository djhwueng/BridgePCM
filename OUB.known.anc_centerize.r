rm(list=ls())
library(ape)
library(sde)
library(geiger)
library(mvMORPH)
library(phangorn)
library(phytools)
library(phyclust)

DriftDiffusion<-function(model.params){
		alpha.offset <<- model.params[1]
		sigma.offset <<- model.params[2]
		drift <- expression(-alpha.offset*x)
		diffusion <- expression(sigma.offset)
		diffusion.x <- expression(0)

		D1<-function(t,x,theta){-theta[1]*x}
		S1<-function(t,x,theta){theta[2]}
		return(list(expression.drift=drift,expression.diffusion=diffusion,expression.diffusion.x=diffusion.x,function.drift=D1,function.diffusion=S1))
		}

OUBridge <- function(model.params, start.value=start.value, end.value=end.value, t0=0, T=T, N=N, drift=drift, diffusion=diffusion, diffusion.x=diffusion.x){
	alpha<-model.params[1]
	OUpath<-sde.sim(X0 = start.value, drift=drift, sigma=diffusion, sigma.x = diffusion.x, N=N)
  t <- seq(t0,T,length=N+1)
	Lambda <- exp(-alpha*(T-t)) - exp(-alpha*(T+t))
	Lambda <- Lambda/(1-exp(-2*alpha*T))
	OUB <- OUpath - Lambda*(OUpath[N+1] - end.value)
	return(as.numeric(OUB))
	}

# model.params<-c(1,1,1)
# drift.diffusion<-DriftDiffusion(model.params)
# start.value<-122
# end.value<-23
# N<-100
# drift<-drift.diffusion$expression.drift
# diffusion<-drift.diffusion$expression.diffusion
# diffusion.x<-drift.diffusion$expression.diffusion.x
# OUBridge(model.params,start.value=start.value, end.value=end.value, t0=0, T=T, N=N, drift=drift,diffusion=diffusion,diffusion.x=diffusion.x)

condOU <- function(model.params, phy=phy, tip.states=tip.states) {
	alpha<-model.params[1]
  sigma<-model.params[2]
	ntip <- length(phy$tip.label)
	N <- dim(phy$edge)[1] + 1
	covmatrix = sigma^2/2/alpha*exp(-alpha*dist.nodes(phy))

	cov22 = covmatrix[1:ntip, 1:ntip]
	cov11 = covmatrix[(ntip+1):N, (ntip+1):N]
	cov12 = covmatrix[1:ntip, (ntip+1):N]

	inv_cov22 = solve(cov22)

	condmean = t(cov12)%*%inv_cov22%*%(tip.states)
	condvar =  cov11 - t(cov12)%*%inv_cov22%*%cov12
	return(list(condmean = condmean, condvar = condvar))
	}

sim.ou.path<-function(model.params, phy=phy, tip.states=tip.states, drift=drift, diffusion=diffusion, diffusion.x=diffusion.x){
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
        assign(paste("path",edgeIndex,sep=""), OUBridge(model.params, start.value=start.state, end.value=end.state, t0=0, T=ceiling(brnlen), N=100, drift=drift,diffusion=diffusion, diffusion.x=diffusion.x))
        path.objects<-c(path.objects, list(get(paste("path",edgeIndex,sep=""))))
			}#loop of edgeIndex
		return(path.objects)
		}#end of function

k.alpha<-function(alpha,t=t){
  return(exp(2*alpha*t-1)/(2*alpha))
  }

ousigmasqhat<-function(alpha,phy=phy,path.data=path.data){
  #alpha<-0.5
  edge.number<-dim(phy$edge)[1]
  sigmasqhat<-0
  for(edgeIndex in 1:edge.number){
    #edgeIndex<-1
    path<-unlist(path.data[edgeIndex])
    for(pathIndex in 2:length(path)){
      ti<-pathIndex/length(path)
      timinus1<-(pathIndex-1)/length(path)
      T<-length(path)/length(path)
      # pathIndex<-2
      if(ti!=T){
        A<-(path[pathIndex]-exp(alpha*(ti-timinus1))*path[pathIndex-1])^2/k.alpha(alpha,t=ti-timinus1)
        B<-exp(2*alpha*(T-ti))*path[pathIndex]^2/k.alpha(alpha,t=(T-ti))
        C<-exp(2*alpha*(T-timinus1))*path[pathIndex]^2/k.alpha(alpha,t= (T-timinus1))
        #if( is.nan(C)){print(c(edgeIndex,A,B,C,k.alpha(alpha,t= (length(path)-(pathIndex-1)))))}
        sigmasqhat <- sigmasqhat + A + B - C
      }
    }
  }
  sigmasqhat
  return(sigmasqhat)
}

ounegloglike<-function(alpha,phy=phy,path.data=path.data){
  edge.number<-dim(phy$edge)[1]
  sigmasq<-ousigmasqhat(alpha,phy=phy,path.data=path.data)
  negloglike<-0
  for(edgeIndex in 1:edge.number){
    #edgeIndex<-1
    path<-unlist(path.data[edgeIndex])
    for(pathIndex in 2:length(path)){
      ti<-pathIndex/length(path)
      timinus1<-(pathIndex-1)/length(path)
      T<-length(path)/length(path)
      
      #pathIndex<-1
      if(ti!=T){   
      D<- 1/2* log(k.alpha(alpha,t=(T-timinus1))/(2*pi*sigmasq*k.alpha(alpha,t=timinus1)*k.alpha(alpha,t=T-ti))) 
      E1<--(path[pathIndex]-exp(alpha*(ti-timinus1))*path[pathIndex-1])^2/(2*sigmasq*k.alpha(alpha,t=ti-timinus1))
      E2<--exp(2*alpha*(T- ti))*path[pathIndex]^2/(2*sigmasq*k.alpha(alpha,t=(T-ti)))
      E3<-exp(2*alpha*(T-timinus1))*path[pathIndex]^2/(2*sigmasq*k.alpha(alpha,t= T-timinus1))
      negloglike <- negloglike + D + E1 + E2 + E3
      }
    }
  }
  return(negloglike)
}

plot.history<-function(phy=phy,path.data=path.data,main=main){
    ntips<-length(phy$tip.label)
 	  edge.number<-dim(phy$edge)[1]
    x.start<-array(0,c(edge.number,1))
    x.end<-array(0,c(edge.number,1))
    step<-array(0,c(edge.number,1))
 	  anc.des.start.end<-cbind(phy$edge,step,x.start, x.end)
 	  colnames(anc.des.start.end)<-c("anc","des","step","x.start","x.end")
    anc.des.start.end<-apply(anc.des.start.end,2,rev)
 	  anc.des.start.end<-data.frame(anc.des.start.end)
    anc<- anc.des.start.end$anc
    des<- anc.des.start.end$des
 	  for(edgeIndex in 1:edge.number){
       path<-unlist(path.data[edgeIndex])
       anc.des.start.end$step[edgeIndex]<-length(path)
 	 		 }

    for(edgeIndex in 1:edge.number){

		  if(anc[edgeIndex]== Ntip(phy)+1){
         anc.des.start.end$x.start[edgeIndex]<- 1+ ceiling(nodeheight(phy,node=anc.des.start.end$anc[edgeIndex]))
         }else{
         			anc.des.start.end$x.start[edgeIndex]<- ceiling(nodeheight(phy,node=anc.des.start.end$anc[edgeIndex]))
			 				}
				 anc.des.start.end$x.end[edgeIndex]<- 1 + ceiling(nodeheight(phy,node=anc.des.start.end$des[edgeIndex]))
        }
    plot( NA,type="l",xlim=c(1,ceiling(get.rooted.tree.height(phy))+1 ),ylim=c(min(unlist(path.data)), max(unlist(path.data))),ylab="Trait value", xlab="Time Steps",main=main)
    for(edgeIndex in 1:edge.number){
      path<-unlist(path.data[edgeIndex])
      starting.x<-ceiling(nodeheight(phy,node=anc[edgeIndex]))
        points((starting.x+1): (starting.x+length(path)),path,type="l")
        }
#    abline(v=ceiling(get.rooted.tree.height(phy))+1)
	  }#end of plot history

model.params<-c(0.1,3)
names(model.params)<-c("alpha","sigma")
size<-10
tip.states<-rnorm(size,sd=1)
tip.states
phy<-rcoal(size)

if(min(phy$edge.length)<1){
  if(min(phy$edge.length)==0){ phy$edge.length[which(phy$edge.length==0)]<-1e-5 }
  phy$edge.length <-  phy$edge.length *(1/min(phy$edge.length)+1)
}

while(min(phy$edge.length)<2){
  phy$edge.length<-1.005*phy$edge.length
  }
phy$edge.length
print(phy$edge.length)
drift.diffusion<-DriftDiffusion(model.params)
expression.drift<-drift.diffusion$expression.drift
expression.diffusion<-drift.diffusion$expression.diffusion
expression.diffusion.x<-drift.diffusion$expression.diffusion.x
# function.drift<-drift.diffusion$function.drift
# function.diffusion<-drift.diffusion$function.diffusion
path.data<-sim.ou.path(model.params, phy=phy, tip.states=tip.states, drift=expression.drift, diffusion=expression.diffusion, diffusion.x=expression.diffusion.x)
#print(path.data)
par(mfrow=c(1,2))
plot(phy)
plot.history(phy=phy,path.data=path.data,main="OUB")
alpha<-0.2
print(ounegloglike(alpha,phy=phy,path.data=path.data))

testou<-optimize(ounegloglike,c(0,10),phy=phy,path.data=path.data)
print(testou)
print(ousigmasqhat(0.1,phy=phy,path.data=path.data)) 
print(model.params)