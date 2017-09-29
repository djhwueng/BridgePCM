rm(list=ls())
library(ape)
library(sde)
library(geiger)
library(mvMORPH)
library(phangorn)
library(phytools)
library(phyclust)

#Brownian motion with drift
#SDE dX_t = mu dt + sigma dW_t
sim.bb.drift.one.path<-function(model.params,T=T,N=N,start.value=start.value){
  mu<-model.params[1]
  sigma<-model.params[2]
  dw<-rnorm(N,0,sqrt(T/N))
  dt<-T/N
  path<-c(start.value)
  for(Index in 2:(N+1)){
    path[Index]<-path[Index-1] + mu*dt + sigma*dw[Index-1]
    }
  return(ts(path,start=start.value,deltat=dt))
  }

# mu<- -5
# sigma<-2
# T=0.345
# N=100
# start.value=10
# path=sim.bb.drift.one.path(c(mu,sigma),T,N,start.value)
# print(path)
# plot(path)

sim.bb.drift.bridge<-function(model.params,start.value=start.value,end.value=end.value,T=T,N=N){
  BB.drift.path<-sim.bb.drift.one.path(model.params,T=T,N=N,start.value=start.value)
  t<-seq(0,T,length=N+1)
  Lambda<-t/T
  BB.drift<-BB.drift.path - Lambda*(BB.drift.path[N+1]-end.value)
  return(as.numeric(BB.drift))
}

sim.bb.drift.tree.path<-function(model.params,phy=phy,tip.states=tip.states,N=N){
  anc.list<-ace(x=tip.states,phy=phy)
  full.node.data<-c(tip.states,anc.list$ace)
  anc<-phy$edge[,1]
  des<-phy$edge[,2]
  ntips<-length(phy$tip.label)
  edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
  path.data<-list()
  for(edgeIndex in edge.number:1){
    brlen<-edge.length[edgeIndex]
    end.value<-full.node.data[des[edgeIndex]]
    start.value<-full.node.data[anc[edgeIndex]]
    assign(paste("path",edgeIndex,sep=""),sim.bb.drift.bridge(model.params,T=brlen,N=N,start.value=start.value, end.value=end.value))
    path.data<-c(path.data,list(get(paste("path",edgeIndex,sep=""))))
    }
  return(path.data)
  }


bbdriftnegloglike.onepath<-function(model.params,path=path,T=T,N=N){
  badval<-(0.5)*.Machine$double.xmax
  mu<-model.params[1]
  sigma<-model.params[2]
  if(sigma<0){
    return(badval)
    }
  dt<-T/N
  negloglike<-N/2*log(dt*sigma^2)
  for(pathIndex in 2:length(path)){
    negloglike<-negloglike+1/2/dt/sigma^2*(path[pathIndex]-path[pathIndex-1]-mu*dt)^2 #see reference problem 3 ,4
    }
  if(!is.finite(negloglike)){return(badval)}
  return(negloglike)
  }

bbdriftnegloglike.tree<-function(model.params,phy=phy,path.data=path.data,N=N){
  badval<-(0.5)*.Machine$double.xmax
  mu<-model.params[1]
  sigma<-model.params[2]
  edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
  if(sigma<0){return(badval)}
  negloglike<-0
  for(edgeIndex in edge.number:1){
    brlen<-edge.length[edgeIndex]
    path<-unlist(path.data[edgeIndex])
    negloglike<-negloglike+bbdriftnegloglike.onepath(model.params,path=path,T=brlen,N=N)
    }
  if(!is.finite(negloglike)){return(badval)}
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

N=10000
model.params<-c(0,1)
names(model.params)<-c("mu","sigma")
tree.size<-3
phy<-rcoal(tree.size)
phy<-reorder(phy,"postorder")
min.length<-N/50
while(min(phy$edge.length)<min.length){
  phy$edge.length<-1.0005*phy$edge.length
}
tip.states<-rnorm(tree.size,mean=10,sd=10)
path.data<-sim.bb.drift.tree.path(model.params,phy=phy,tip.states=tip.states,N=N)
par(mfrow=c(1,2))
plot(phy)
plot.history.dt(phy=phy,path.data=path.data,main="BB drift")
print(bbdriftnegloglike.tree(model.params,phy=phy,path.data=path.data,N=N))
print(optim(par=model.params, bbdriftnegloglike.tree, method="Nelder-Mead",phy=phy,path.data=path.data,N=N))
