#use bb_oub_100.r, CIRbridge.R, BB.known.anc.r (ace), OUB.known.anc.correct.with.dt.r (condOU), plot.history_dt.r to generate this plot

rm(list=ls())
library(ape)
library(sde)
library(phytools)
library(phyclust)

#plot tree using path space
PlotTreeHistory<-function(phy=phy,pathdata=pathdata,main=main){
  ntips<-length(phy$tip.label)
  edge.number<-dim(phy$edge)[1]
  start<-rep(0,length(edge.number))
  end<-rep(0,length(edge.number))
  
  mtx<-(cbind(phy$edge,ceiling(phy$edge.length)+1,start,end))
  mtx<-mtx[dim(mtx)[1]:1,]
  colnames(mtx)<-c("anc","des","step","start","end")
  mtx<-as.data.frame(mtx)
  mtx$start[1:2]<-1
  mtx$end[1:2]<-mtx$step[1:2]
  for(edgeIndex in 3:dim(mtx)[1]){
    mtx$start[edgeIndex] <-  mtx$end[mtx$anc[edgeIndex]==mtx$des]
    mtx$end[edgeIndex] <-    mtx$end[mtx$anc[edgeIndex]==mtx$des]+mtx$step[edgeIndex]-1
  }
  xlim<-c(1, max(mtx$end) )
  unlistpath<-unlist(pathdata)
  ylim<-c(min(unlistpath),max(unlistpath))
  plot(NA,xlim=xlim,ylim=ylim,ylab="Trait value", xlab="Times Steps",main=main)
  points(1, pathdata[[1]][1] ,lwd=5,col= 1)
  for(edgeIndex in 1:edge.number){
    path<-unlist(pathdata[edgeIndex]) 
    points(mtx$start[edgeIndex]:mtx$end[edgeIndex], path,type="l",lwd=2)
    points(mtx$end[edgeIndex]:mtx$end[edgeIndex], path[length(path)],lwd=5,col= edgeIndex+1)
  }    
}#end of plot history


#Brownian bridge
bb<-function(sigma,a=0,b=0,t0=0,T=1,N=100){
  if(T<=t0){stop("wrong times")}
  dt <- (T-t0)/N
  t<- seq(t0,T,length=N+1)
  X<-c(0,sigma*cumsum(rnorm(N)*sqrt(dt)))
  bbpath<-a+X-(t-t0)/(T-t0)*(X[N+1]-b+a)
  return(bbpath)
  }
  
bbtree<-function(sigma,phy=phy,tips=tips){
  anclist<-ace(x=tips,phy=phy)
  fullnodedata<-c(tips,anclist$ace)
  anc<-phy$edge[,1]
  des<-phy$edge[,2]
  ntips<-length(phy$tip.label)
  edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
  pathdata<-list()
  for(edgeIndex in edge.number:1){
    brlen<-edge.length[edgeIndex]
    N<-ceiling(brlen)
    b<-fullnodedata[des[edgeIndex]]
    a<-fullnodedata[anc[edgeIndex]]
    assign(paste("path",edgeIndex,sep=""),bb(sigma,a=a,b=b,t0=0,T=brlen,N=N))
    pathdata<-c(pathdata, list(get(paste("path",edgeIndex,sep=""))) )
  }
  return(pathdata)
}


#Ornstein Uhlenbeck bridge
oub<-function(Theta, a=0,b=0,t0=0,T=1,N=100){
  alpha<-Theta[1]
  theta<-Theta[2]
  sigma<-Theta[3]
  if(T<=t0){stop("wrong times")}
  dt<-(T-t0)/N
  t<-seq(t0,T,length=N+1)
  d<-expression(-1*x)#alp
  s<-expression(2)#sigma=2
  X<-sde.sim(X0=a,drift=d,sigma=s,N=N)
  Lambda<-exp(-alpha*(T-t))-exp(-alpha*(T+t))
  Lambda<-Lambda/(1-exp(-2*alpha*T))
  oubpath<-X-Lambda*(X[N+1]-b)
  oubpath<-as.numeric(oubpath)
  return(oubpath)
}

condOU<-function(Theta,phy=phy,tips=tips){
  alpha<-Theta[1]
  theta<-Theta[2]
  sigma<-Theta[3]
  n<-length(phy$tip.label)
  m<-dim(phy$edge)[1]+1
  covmtx<-sigma^2/2/alpha*exp(-alpha*dist.nodes(phy))
  
  cov22 <- covmtx[1:n,1:n]
  cov11<-covmtx[(n+1):m,(n+1):m]
  cov12<-covmtx[1:n,(n+1):m]
  
  inv_cov22 <-solve(cov22)
  
  condmean <- theta + t(cov12)%*%inv_cov22%*%(tips-rep(theta,n))
  condvar <- cov11 - t(cov12)%*%inv_cov22%*%cov12
  return(list(condmean=condmean,condvar=condvar))
}



oubtree<-function(Theta,phy=phy,tips=tips){
  #anclist<-condOU(Theta=Theta,phy=phy,tips=tips)
  #fullnodedata<-c(tips,anclist$condmean)
  anclist<-ace(x=tips,phy=phy)
  fullnodedata<-c(tips,anclist$ace)
  anc<-phy$edge[,1]
  des<-phy$edge[,2]
  ntips<-length(phy$tip.label)
  edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
  pathdata<-list()
  for(edgeIndex in edge.number:1){
    brlen<-edge.length[edgeIndex]
    N<-ceiling(brlen)
    b<-fullnodedata[des[edgeIndex]]
    a<-fullnodedata[anc[edgeIndex]]
    assign(paste("path",edgeIndex,sep=""),oub(Theta,a=a,b=b,t0=0,T=brlen,N=N))
    pathdata<-c(pathdata, list(get(paste("path",edgeIndex,sep=""))) )
  }
  return(pathdata)
}


#CIR bridge
A<-function(x,t,alpha=alpha,theta=theta,b=b,T=T){
  jx<- 4*alpha*sqrt(x*b*exp(-alpha*(T-t)))/(sigma^2*(1-exp(-alpha*(T-t))))
  J1<-besselI(x=jx,nu=2*alpha*theta/sigma^2)
  J2<-besselI(x=jx, nu=2*alpha*theta/sigma^2-1)
  J1J2<-J1/J2
  if(!is.finite(J1) || !is.finite(J2)){J1J2 <-1} 
  Ax<- alpha*(theta-x) + 2*alpha*x/(1-exp(-alpha*(T-t)))*(sqrt(b/x*exp(-alpha*(T-t)))*J1J2-exp(-alpha*(T-t)))
  return(Ax)
}


B<-function(x,sigma=sigma,p=p){
  return(sigma^2*x^p)
}

cirb<-function(Theta,a=0,b=0,t0=0,T=1,N=100){
  alpha<-Theta[1]
  theta<-Theta[2]
  sigma<-Theta[3]

  noise<-rnorm(n=N-1,mean=0,sd=sqrt(T/N))
  cirbpath<-array(0,c(N+1))
  cirbpath[1]<-a
  
  for(Index in 2:N){
    cirbpath[Index] <-abs(cirbpath[Index-1] + A(cirbpath[Index-1],(Index-2)*T/N,alpha=alpha,theta=theta,b=b,T=T)*T/N + sqrt(B(cirbpath[Index-1],sigma=sigma,p=0.6))*noise[Index-1]) 
    }
  cirbpath[N+1]<-b
  return(cirbpath)
  }

cirbtree<-function(Theta,phy=phy,tips=tips){
  #anclist<-condOU(Theta=Theta,phy=phy,tips=tips)
  #fullnodedata<-c(tips,anclist$condmean)
  anclist<-ace(x=tips,phy=phy)
  fullnodedata<-c(tips,anclist$ace)
  anc<-phy$edge[,1]
  des<-phy$edge[,2]
  ntips<-length(phy$tip.label)
  edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
  pathdata<-list()
  for(edgeIndex in edge.number:1){
    brlen<-edge.length[edgeIndex]
    N<-ceiling(brlen)
    b<-fullnodedata[des[edgeIndex]]
    a<-fullnodedata[anc[edgeIndex]]
    assign(paste("path",edgeIndex,sep=""),cirb(Theta,a=a,b=b,t0=0,T=brlen,N=N))
    pathdata<-c(pathdata, list(get(paste("path",edgeIndex,sep=""))) )
  }
  return(pathdata)
}
#Main program

T<-1
N=1000
t0<-0
dt<-1
t<-seq(t0,T,length=N+1)

a=2
b=10

#Brownian bridge
sigma=2
plot(bb(sigma,a=a,b=b,t0=t0,T=T,N=N),type="l")

#OU bridge
Theta=c(1,1,2) #(alp,theta,sigma)
plot(oub(Theta=Theta,a=a,b=b,t0=t0,T=T,N=N),type="l")

#CIR bridge
Theta=c(1,1,2)
plot(cirb(Theta=Theta,a=a,b=b,t0=t0,T=T,N=N),type="l")


#tree
size<-4
tips<-c(0.5,1,2,4)
phy<-rcoal(size)
phy<-reorder(phy,"postorder")
phy$edge.length
min.length<-10
while(min(phy$edge.length)<min.length){
  phy$edge.length<-1.005*phy$edge.length
}
phy$tip.label<-paste("Spe.",1:length(phy$tip.label),sep="")
phy$edge.length
plot(phy)

#BB tree
sigma=2
bbpathdata<-bbtree(sigma,phy=phy,tips=tips)
par(mfrow=c(2,2))
plot(phy,edge.width=3,cex=1)
#tiplabels(pch=21,col="black",adj=1,bg="black",cex=2)
PlotTreeHistory(phy=phy,pathdata=bbpathdata,main="Brownian Bridge")

## OUB tree
Theta=c(0.01,mean(tips),1) #(alp,theta,sigma)
oubpathdata<-oubtree(Theta,phy=phy,tips=tips)
PlotTreeHistory(phy=phy,pathdata=oubpathdata,main="OU Bridge")

## CIRB tree
Theta=c(0.1,mean(tips),1) #(alp,theta,sigma)
cirbpathdata<-cirbtree(Theta,phy=phy,tips=tips)
PlotTreeHistory(phy=phy,pathdata=cirbpathdata,main="CIR Bridge")

