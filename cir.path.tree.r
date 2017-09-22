rm(list=ls())
library(ape)
library(sde)

sim.cir.path <- function(model.params,T=T,N=N,x0=x0){
    mu <- model.params[1]
    alpha<-model.params[2]
    sigma<-model.params[3]
    dw<- rnorm(N,0,sqrt(T/N))
    dt<-T/N
    path<-c(x0)
    for(Index in 2:(N+1)){
      path[Index]<-path[Index-1]+alpha*(mu-path[Index-1])*dt +sigma*sqrt(path[Index-1])*dw[Index-1]
      }
    if(2*alpha*mu-sigma^2>=0){
    return(ts(path,start=x0,deltat=dt))}else{return(NULL)}
    }

sim.cir.path.tree<-function(model.params,phy=phy,N=N,root=root){
  sim.node.data<-integer(length(phy$tip.label)+phy$Nnode)
  edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
  ntips<-length(phy$tip.label)
  ROOT<-ntips+1
  anc<- phy$edge[,1]
  des<- phy$edge[,2]
  path.objects<-list()
  start.state<-root
  for(edgeIndex in edge.number:1){
    brnlen<-edge.length[edgeIndex]
#    print(anc[edgeIndex])
#    print(start.state)
    start.state<-sim.node.data[anc[edgeIndex]]
    assign(paste("path",edgeIndex,sep=""), sim.cir.path(model.params,T=brnlen,N=N,x0=start.state))
    temp.path<-get(paste("path",edgeIndex,sep=""))
    sim.node.data[des[edgeIndex]]<-temp.path[length(temp.path)]
    print(sim.node.data[des[edgeIndex]])
  #  print(sim.node.data)
    path.objects<-c(path.objects,list(get(paste("path",edgeIndex,sep=""))))
    }
  return(list(path.objects=path.objects,sim.node.data=sim.node.data))
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
         path<-unlist(path.data$path.objects[edgeIndex])
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
      main="cir"
      plot( NA,type="l",xlim=c(1,ceiling(get.rooted.tree.height(phy))+1 ),ylim=c(min(unlist(path.data)), max(unlist(path.data))),ylab="Trait value", xlab="Time Steps",main=main)
      anc.des.start.end
      for(edgeIndex in 1:edge.number){
        edgeIndex <- 1
        path<-unlist(path.data$path.objects[edgeIndex])
        #starting.x<-ceiling(nodeheight(phy,node=anc[edgeIndex]))
        x.start <-  anc.des.start.end$x.start[edgeIndex]
        x.end   <-  anc.des.start.end$x.end[edgeIndex]
        gap <- round(anc.des.start.end$step[edgeIndex]/( x.end-x.start))
        #points((starting.x+1): (starting.x+length(path)),path[seq(1, anc.des.start.end$step[edgeIndex], by = gap)],type="l")
        points(  (anc.des.start.end$x.start[edgeIndex]: anc.des.start.end$x.end[edgeIndex]),


        path[seq(1, anc.des.start.end$step[edgeIndex], by = gap)  # we need to find a better way to sample this and get right number with x - axis



        ],type="l")
          }
  #    abline(v=ceiling(get.rooted.tree.height(phy))+1)
  	  }#end of plot history


T<-1
N<-3000
x0<-0.1
sims<-1000
true.alpha<-0.01
true.mu<-1
true.sigma<-0.1
model.params<-c(true.mu,true.alpha,true.sigma)
#plot(sim.cir.path(model.params, T=T, N=N,x0=0.1))
root<-0
# we can do sims to different regimes
size<-3
phy<-rcoal(size)
phy<-reorder(phy,"postorder")
min.length<-15
while(min(phy$edge.length)<min.length){
  phy$edge.length<-1.005*phy$edge.length
  }
phy$edge.length
plot(phy)
path.data<-sim.cir.path.tree(model.params,phy=phy,N=N,root=root)
#print(path.data)
plot.history(phy=phy,path.data=path.data,main="CIR")
