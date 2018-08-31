rm(list=ls())

library(ape)
library(sde)
library(geiger)
library(mvMORPH)
library(phytools)
library(phangorn)
library(phyclust)

modified.BBridge <- function(sigma, x = 0, y = 0, t0 = 0, T = 1, N = 100){
  if (T <= t0)
        stop("wrong times")
    dt <- 1# (T - t0)/N #dt will be difference given difference brnlen, but overall we want 100 step on each brnlen
    # 
    # x=start.state
    # y=end.state
    # t0=0
    # T=brnlen
    # N=100
    # sigma<-2
    #print(dt)
    t <- seq(t0, T, length = N+1)
    X <- c(0, cumsum(rnorm(N) * sqrt(dt)))
    BB <- x + X - (t - t0)/(T - t0) * (X[N + 1] - y + x)
#    print("here is sigma")
#    print(sigma)
#    print(BB)
#    print(BB)
#    BB<-c(sigma*BB)
#    X <- ts(BB, start = t0, deltat = dt)
#    return(invisible(X))
     return(sigma*BB)
     }

# sigma<-12
# single.path<-modified.BBridge(sigma, x = 0, y = 0, t0 = 0, T = 1, N = 100)
# #sd(single.path)
# plot(single.path)
# 
# testloglike<-function(sigma){
#   negloglike<-0
#   path<-unlist(single.path)
#   for(pathIndex in 2:length(path)){
#       negloglike<-negloglike+ 1/2*log(2*pi)+1/2*2*log(sigma)+1/(2*sigma^2)*(path[pathIndex]-path[pathIndex-1])^2
#       }
#   return(negloglike)
#   }
# 
# optimize(testloglike,c(0,100))

sim.bb.path<-function(sigma,tip.states=tip.states,phy=phy,SE=NULL){ #SE is for tips
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
      assign(paste("path",edgeIndex,sep=""),modified.BBridge(sigma,x=start.state,y=end.state,t0=0,T=brnlen,N=100))}else{
      assign(paste("path",edgeIndex,sep=""),modified.BBridge(sigma,x=start.state,y=end.state,t0=0,T=brnlen,N=100))
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

bmnegloglike<-function(sigma,phy=phy,path.data=path.data){
      edge.number<-dim(phy$edge)[1]
      negloglike<-0
      for(edgeIndex in 1:edge.number){
        path<-unlist(path.data[edgeIndex])
      	 for(pathIndex in 2:length(path)){
             negloglike<-negloglike+ 1/2*log(2*pi)+1/2*2*log(sigma)+1/(2*sigma^2)*(path[pathIndex]-path[pathIndex-1])^2
             }
#        print(negloglike)
        }
      return(negloglike)
      }#end of loglike function

size<-30
phy<-rcoal(size)
while(min(phy$edge.length)<10){
   phy$edge.length<-1.005*phy$edge.length
   }
phy$edge.length
plot(phy)
tip.states<-rnorm(size,sd=10)

if(min(phy$edge.length)<1){
  if(min(phy$edge.length)==0){ phy$edge.length[which(phy$edge.length==0)]<-1e-5 }
  phy$edge.length <-  phy$edge.length *(1/min(phy$edge.length)+1)
  }

sigma<-200
tip.states<-rnorm(size,sd=10)
path.data<-sim.bb.path(sigma,tip.states=tip.states,phy=phy)

plot.history(phy=phy,path.data=path.data,main="Test BB") #finish plotting
print(bmnegloglike(sigma,phy=phy,path.data=path.data))
optimize(bmnegloglike,c(0,1000),phy=phy,path.data=path.data)


####Data 1 good, paramter estimate is closed to geiger
phy <- read.tree("~/Dropbox/CollabHoJhwueng/data/Molina_Borja.et.al_2003/tree.phy")
plot(phy)

# if(min(phy$edge.length)<1){
#   if(min(phy$edge.length)==0){ phy$edge.length[which(phy$edge.length==0)]<-1e-5 }
#   phy$edge.length <-  phy$edge.length *(1/min(phy$edge.length)+1)
#   }
# 
# while(min(phy$edge.length)<10){
#   phy$edge.length<-1.005*phy$edge.length
#   }

plot(phy)
traits<-read.csv("~/Dropbox/CollabHoJhwueng/data/Molina_Borja.et.al_2003/trait.csv")
head(traits)
tip.states<-traits$HatchingMass
names(tip.states)<-phy$tip.label
geiger.results<-fitContinuous(phy, tip.states, model="BM")
geiger.results$lik
sigma<-geiger.results$opt$sigsq
sigma
path.data<-sim.bb.path(sigma,tip.states=tip.states,phy=phy)
result <-optimize(bmnegloglike,c(0,100),phy=phy,path.data=path.data)
result   
 
####Data 2 good, paramter estimate is closed to geiger
phy <- read.tree("~/Dropbox/CollabHoJhwueng/data/Niewiarowski.et.al_2004/tree.phy")
plot(phy)

# if(min(phy$edge.length)<1){
#   if(min(phy$edge.length)==0){ phy$edge.length[which(phy$edge.length==0)]<-1e-5 }
#   phy$edge.length <-  phy$edge.length *(1/min(phy$edge.length)+1)
# }
# 
# while(min(phy$edge.length)<10){
#   phy$edge.length<-1.005*phy$edge.length
# }

plot(phy)
traits<-read.csv("/Users/djhwueng/Dropbox/CollabHoJhwueng/data/Niewiarowski.et.al_2004/trait.csv")
head(traits)
tip.states<-traits$egg.mass
names(tip.states)<-phy$tip.label
geiger.results<-fitContinuous(phy, tip.states, model="BM")
sigma<-geiger.results$opt$sigsq
sigma
path.data<-sim.bb.path(sigma,tip.states=tip.states,phy=phy)
result <-optimize(bmnegloglike,c(0,100),phy=phy,path.data=path.data)
result   

####Data 3 good, paramter estimate is closed to geiger
phy <- read.tree("~/Dropbox/CollabHoJhwueng/data_we_will_use/webster.purvis_length/intree")
plot(phy)

#if(min(phy$edge.length)<1){
#  if(min(phy$edge.length)==0){ phy$edge.length[which(phy$edge.length==0)]<-1e-5 }
#  phy$edge.length <-  phy$edge.length *(1/min(phy$edge.length)+1)
#}

#while(min(phy$edge.length)<10){
#  phy$edge.length<-1.005*phy$edge.length
#}

plot(phy)
traits<-read.csv("~/Dropbox/CollabHoJhwueng/data_we_will_use/webster.purvis_length/purvis.length.csv")
head(traits)
tip.states<-traits$length.mm.
names(tip.states)<-phy$tip.label
geiger.results<-fitContinuous(phy, tip.states, model="BM")
sigma<-geiger.results$opt$sigsq
sigma
path.data<-sim.bb.path(sigma,tip.states=tip.states,phy=phy)
result <-optimize(bmnegloglike,c(0,100),phy=phy,path.data=path.data)
result   

