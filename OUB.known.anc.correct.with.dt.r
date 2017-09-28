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

#we estimate ancestor first, so we have nodes value for tree. Then we do bridge likelihood for parameter estimation.
sim.ou.tree.path<-function(model.params, phy=phy, tip.states=tip.states, N=N, SE=NULL){
	ntips<-length(phy$tip.label)
	edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
	anc<- phy$edge[,1]
	des<- phy$edge[,2]
	anc.states<-condOU(model.params,phy=phy,tip.states=tip.states)$condmean
  full.node.data<-c(tip.states,anc.states)
	path.objects<- list()
	for(edgeIndex in 1:edge.number){
      brlen<-edge.length[edgeIndex]
      end.value<-full.node.data[des[edgeIndex]]
			start.value<-full.node.data[anc[edgeIndex]]
#        assign(paste("path",edgeIndex,sep=""), modified.DBridge(x = start.state, y = end.state, t0 = 0, N=ceiling(brnlen), drift=drift, sigma=sigma, sigma.x= 0, N.give.up=10))
      assign(paste("path",edgeIndex,sep=""), sim.ou.bridge(model.params, start.value=start.value, end.value=end.value, T=brlen, N=N))
      path.objects<-c(path.objects, list(get(paste("path",edgeIndex,sep=""))))
			}#loop of edgeIndex
		return(path.objects)
		}#end of function


#
# first order work very well for each branch as it is the 1D case.
# but likelihood given paths of tree yet to determine.
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
#when move to tree case, we shall be able to use first order method, however, yet to find a method to do.
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


ounegloglike.onepath<-function(model.params,path=path,T=T,N=N){
  badval<-(0.5)*.Machine$double.xmax
	mu<-model.params[1]
	alpha<-model.params[2]
	sigma<-model.params[3]
	if(alpha<0 || alpha>100 || sigma<0){
		return(badval)
		}
	dt<-T/N
	negloglike<- N/2*log(sigma^2/(2*alpha))
	for(Index in 2:(N+1)) {
		negloglike <- negloglike + 1/2*log(1-exp(-2*alpha*dt))
		negloglike <- negloglike + alpha/sigma^2*(path[Index] - mu - (path[Index-1] - mu)*exp(-alpha*dt))^2/(1-exp(-2*alpha*dt))
		}
	if (!is.finite(negloglike)){
		return(badval)
		}
		return(negloglike)
	}

ounegloglike.tree<-function(model.params,phy=phy,path.data=path.data,N=N){
	badval<-(0.5)*.Machine$double.xmax
	mu<-model.params[1]
	alpha<-model.params[2]
	sigma<-model.params[3]
	edge.number<-dim(phy$edge)[1]
	edge.length<-phy$edge.length
	negloglike<-0
	for(edgeIndex in edge.number:1){
		brlen<-edge.length[edgeIndex]
		path<-unlist(path.data[edgeIndex])
		negloglike<-negloglike+ounegloglike.onepath(model.params,path=path,T=brlen,N=N)
		}
  if(alpha<0 || alpha > 100 || sigma<0 || !is.finite(negloglike) ){
		return(badval)
		}
  return(negloglike)
	}


plot.history.dt<-function(phy=phy,path.data=path.data,main=main){
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
    plot( NA,type="l",xlim=c(1,ceiling(get.rooted.tree.height(phy))+1 ),ylim=c(min(unlist(path.data)), max(unlist(path.data))),ylab="Trait value", xlab="Time Steps",main=main)
    abline(a=1e-8,b=0,lty=2)

    for(edgeIndex in 1:edge.number){
      if(anc[edgeIndex]== Ntip(phy)+1){
         anc.des.start.end$x.start[edgeIndex]<- 1+ round(nodeheight(phy,node=anc.des.start.end$anc[edgeIndex]))
         }else{
              anc.des.start.end$x.start[edgeIndex]<- 1+ceiling(nodeheight(phy,node=anc.des.start.end$anc[edgeIndex]))
              }
              anc.des.start.end$x.end[edgeIndex]<- 1 + ceiling(nodeheight(phy,node=anc.des.start.end$des[edgeIndex]))
          }
    anc.des.start.end
    for(edgeIndex in 1:edge.number){
      path<-unlist(path.data[edgeIndex])
      #starting.x<-ceiling(nodeheight(phy,node=anc[edgeIndex]))
      x.start <-  anc.des.start.end$x.start[edgeIndex]
      x.end   <-  anc.des.start.end$x.end[edgeIndex]
      gap <- round(anc.des.start.end$step[edgeIndex]/( x.end-x.start))
      #points((starting.x+1): (starting.x+length(path)),path[seq(1, anc.des.start.end$step[edgeIndex], by = gap)],type="l")
      point.to.use<-(anc.des.start.end$x.start[edgeIndex]: anc.des.start.end$x.end[edgeIndex])
      sample.path<-path[seq(1,  anc.des.start.end$step[edgeIndex],by=gap)]
#        print("before")
#        print( c(length(point.to.use), length(sample.path) ) )
      if(length(point.to.use)!=length(sample.path)){
          if(length(point.to.use) > length(sample.path)){
            point.to.use<-point.to.use[1:length(sample.path)]
            #sample.path<-c(sample.path,sample.path[length(sample.path)])
          }else{
            sample.path<-sample.path[1:length(point.to.use)]
            #point.to.use<-c(point.to.use, point.to.use[length(point.to.use)]+1 )
            }
        }
#        print("After")
#        print( c(length(point.to.use), length(sample.path) ) )
      points(point.to.use , sample.path, type="l",lwd=1.5)
      #path[seq(1, anc.des.start.end$step[edgeIndex], by = gap)  # we need to find a better way to sample this and get right number with x - axis
       #path[sample(    )    ]  LOOK UP SAMPLE FOR GOOD
        }
    }#end of plot history


N=2000
model.params<-c(10,0.1,4)
names(model.params)<-c("mu","alpha","sigma")
tree.size<-3
phy<-rcoal(tree.size)
phy<-reorder(phy,"postorder")
min.length<-N/50
while(min(phy$edge.length)<min.length){
  phy$edge.length<-1.005*phy$edge.length
  }
tip.states<-rnorm(tree.size,mean=10,sd=10)
path.data<-sim.ou.tree.path(model.params, phy=phy, tip.states=tip.states, N=N, SE=NULL)
#source("~/GitHub/BridgePCM/plot_history_dt.r")
par(mfrow=c(1,2))
plot(phy)
plot.history.dt(phy=phy,path.data=path.data,main="OUB")
#path.data
print(first.order.ounegloglike.tree(0.1,phy=phy,path.data=path.data,N=N))
print(ounegloglike.tree(model.params,phy=phy,path.data=path.data,N=N))

print(optim(par=model.params,ounegloglike.tree, method="Nelder-Mead",phy=phy,path.data=path.data,N=N))
#path.data
#phy$edge.length

# #?optimize
#optimize(first.order.ounegloglike.tree,c(0,1000000),phy=phy,path.data=path.data,N=N)
# #no good for alphas, we do not have
#
# edge.number<-dim(phy$edge)[1]
# edge.length<-phy$edge.length
# negloglike.array<-array(0,c(edge.number))
# alpha<-0.1
# for(edgeIndex in edge.number:1){
#   path<-unlist(path.data[edgeIndex])
#   brlen<-edge.length[edgeIndex]
#   mu.hat<-first.order.mu(alpha,path=path,T=brlen,N=N)
#   sigma.sq.hat<-first.order.sigma.sq(alpha,path=path,T=brlen,N=N)
#   print(c(mu.hat,sigma.sq.hat))
# }
#single path case estimate well, but tree case not.......yet
