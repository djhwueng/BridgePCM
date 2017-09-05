rm(list=ls())
library(ape)
library(sde)
library(geiger)
library(mvMORPH)
library(phytools)
library(phangorn)
library(phyclust)
modified.BBridge<-function (model.params,x = 0, y = 0, t0 = 0, T = 1, N = 100){
    if (T <= t0)
        stop("wrong times")
    dt <- 1#(T - t0)/N #dt will be difference given difference brnlen, but overall we want 100 step on each brnlen
    t <- seq(t0, T, length = N+1)
    X <- c(0, model.params*cumsum(rnorm(N) * sqrt(dt)))
    BB <- x + X - (t - t0)/(T - t0) * (X[N + 1] - y + x)
    X <- ts(BB, start = t0, deltat = dt)
    return(invisible(X))
    }

sim.bb.path<-function(model.params,tip.states=tip.states,phy=phy,SE=NULL){ #SE is for tips
    anc.list<-ace(x=tip.states,phy=phy)
    full.node.data<- c(tip.states, anc.list$ace) #trait data and estimated ancestral values
    anc<-phy$edge[,1]
    des<-phy$edge[,2]
    ntips<-length(phy$tip.label)
    edge.number <-dim(phy$edge)[1]
    edge.length<-phy$edge.length
    path.data <- list()

    for(edgeIndex in edge.number:1){
      brnlen<- edge.length[edgeIndex]    #here could be amplified for checking issues
      end.state<- full.node.data[des[edgeIndex]] #can do sampling here
      start.state<-full.node.data[anc[edgeIndex]] #can do sampling here
      if(anc[edgeIndex]!=Ntip(phy)+1){
      assign(paste("path",edgeIndex,sep=""),modified.BBridge(model.params,x=start.state,y=end.state,t0=0,T=ceiling(brnlen)-1,N=ceiling(brnlen)-1))}else{
        assign(paste("path",edgeIndex,sep=""),modified.BBridge(model.params,x=start.state,y=end.state,t0=0,T=ceiling(brnlen),N=ceiling(brnlen)))
      }
      path.data<-c(path.data,list(get(paste("path",edgeIndex,sep=""))))
      }
    return(path.data)
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
         if( anc[edgeIndex]== Ntip(phy)+1){
             anc.des.start.end$x.start[edgeIndex]<- 1+ ceiling(nodeheight(phy,node=anc.des.start.end$anc[edgeIndex]))
               }else{
             anc.des.start.end$x.start[edgeIndex]<- ceiling(nodeheight(phy,node=anc.des.start.end$anc[edgeIndex]))
                 }
             anc.des.start.end$x.end[edgeIndex]<-   1+ceiling(nodeheight(phy,node=anc.des.start.end$des[edgeIndex]))
             }
       anc.des.start.end

       main="BB"
       plot( NA,type="l",xlim=c(1,ceiling(get.rooted.tree.height(phy))+1 ),ylim=c(min(unlist(path.data)), max(unlist(path.data))),ylab="Trait value", xlab="Time Steps",main=main)
      for(edgeIndex in 1:edge.number){
         path<-unlist(path.data[edgeIndex])
         length(path)
         starting.x<-ceiling(nodeheight(phy,node=anc[edgeIndex]))
#         starting.x
#         edgeIndex<-1
         #if(anc[edgeIndex]  == Ntip(phy)+1  ){
           points((starting.x+1): (starting.x+length(path)),path,type="l")
 	        # }else{
 	         #  points((starting.x): (starting.x+length(path)-1),path,type="l")
 	        #  }
           }

       abline(v=ceiling(get.rooted.tree.height(phy))+1)

      }#end of plot history

bmnegloglike<-function(model.params,phy=phy,path.data=path.data){
      sigma2<-model.params
      edge.number<-dim(phy$edge)[1]
      negloglike<-0
      for(edgeIndex in 1:edge.number){
      	 path<-unlist(path.data[edgeIndex])
      	 for(pathIndex in 2:length(path)){
             negloglike<-negloglike+ 1/2*log(2*pi)+1/2*log(sigma2)+1/(2*sigma2)*(path[pathIndex]-path[pathIndex-1])^2
             }
      	 	}
      return(negloglike)
      }#end of loglike function



model.params <- 1
names(model.params)<-"sigma2"
size<-10
phy<-rcoal(size)
while(min(phy$edge.length)<10){
  phy$edge.length<-1.005*phy$edge.length
  }
phy$edge.length
plot(phy)
tip.states<-rnorm(size,sd=10)
path.data<-sim.bb.path(model.params,tip.states=tip.states,phy=phy)
plot.history(phy=phy,path.data=path.data,main="Test BB") #finish plotting
print(bmnegloglike(model.params,phy=phy,path.data=path.data))
optimize(bmnegloglike,c(0,10) ,phy=phy,path.data=path.data)
