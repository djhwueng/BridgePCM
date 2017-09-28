
rm(list=ls())
library(ape)
library(geiger)
library(mvMORPH)
library(phytools)
library(phangorn)
library(phyclust)

#Euler scheule
sim.bb.one.path<-function (sigma, T = T, N = N, start.value=start.value, end.value=end.value){
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

sim.bb.tree.path<-function(sigma, phy=phy, tip.states=tip.states, N=N){
  anc.list<-ace(x=tip.states,phy=phy)
  full.node.data<-c(tip.states,anc.list$ace)
  anc<-phy$edge[,1]
  des<-phy$edge[,2]
  ntips<-length(phy$tip.label)
  edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
  path.objects<-list()
  for(edgeIndex in edge.number:1){
#    edgeIndex<-1
    brlen<-edge.length[edgeIndex]
    end.value<-full.node.data[des[edgeIndex]]
    start.value<-full.node.data[anc[edgeIndex]]
    assign(paste("path",edgeIndex,sep=""),sim.bb.one.path(sigma,T=brlen,N=N,start.value=start.value, end.value=end.value))
    path.objects<-c(path.objects,list(get(paste("path",edgeIndex,sep=""))))
    }
  return(path.objects)
  }

bbnegloglike.onepath<-function(sigma,path=path,T=T,N=N){
  badval<-(0.5)*.Machine$double.xmax
  if(sigma<0){return(badval)}
  dt<-T/N
  negloglike<- 1/2*log(2*pi*dt)+1/2*log(sigma^2)
  for(pathIndex in 2:length(path)){
    negloglike <- negloglike+ 1/(2*dt*sigma^2)*(path[pathIndex]-path[pathIndex-1])^2
    }
  return(negloglike)
  }

bbnegloglike.tree<-function(sigma,phy=phy,path.data=path.data,N=N){
  badval<-(0.5)*.Machine$double.xmax
  edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
  negloglike<-0
  if(sigma<0 || !is.finite(negloglike)){return(badval)}
  for(edgeIndex in edge.number:1){
    brlen<-edge.length[edgeIndex]
    path<-unlist(path.data[edgeIndex])
    negloglike<-negloglike+bbnegloglike.onepath(sigma,path=path,T=brlen,N=N)
    }
  return(negloglike)
  }

bbnegloglike.tree.v2<-function(sigma,phy=phy,path.data=path.data,N=N){
  badval<-(0.5)*.Machine$double.xmax
  edge.number<-dim(phy$edge)[1]
  edge.length<-phy$edge.length
  negloglike<-0
  if(sigma<0 || !is.finite(negloglike)){return(badval)}
  for(edgeIndex in edge.number:1){
    brlen<-edge.length[edgeIndex]
    path<-unlist(path.data[edgeIndex])
    for(pathIndex in 2:length(path)){
    negloglike<-negloglike+1/2*log(2*pi)+1/2*log(sigma^2)+1/(2*sigma^2)*(path[pathIndex]-path[pathIndex-1])^2
     }
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

N<-2000
sigma<-10
tree.size<-3
phy<-rcoal(tree.size)
phy<-reorder(phy,"postorder")
min.length<-N/50
while(min(phy$edge.length)<min.length){
  phy$edge.length<-1.005*phy$edge.length
  }

tip.states<-rnorm(tree.size,mean=10,sd=10)
path.data<-sim.bb.tree.path(sigma, phy=phy, tip.states=tip.states, N=N)
par(mfrow=c(1,2))
plot(phy)
tiplabels(pch=21, col="black", adj=1, bg="black", cex=2)
plot.history.dt(phy=phy,path.data=path.data,main="BB")
print(bbnegloglike.tree(sigma,phy=phy,path.data=path.data,N=N))
print(optimize(bbnegloglike.tree,c(0,100),phy=phy,path.data=path.data,N=N))

print(optimize(bbnegloglike.tree.v2,c(0,100),phy=phy,path.data=path.data,N=N))
