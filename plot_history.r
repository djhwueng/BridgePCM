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
    abline(v=ceiling(get.rooted.tree.height(phy))+1)
	  }#end of plot history
