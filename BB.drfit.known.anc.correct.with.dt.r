#Brownian motion with drift
#SDE dX_t = mu dt + sigma dW_t
sim.bbd.one.path<-function(model.params,T=T,N=N,start.value=start.value){
  mu<-model.params
}


T<-1
N<-100
sigma<-1
W<-rep(NA,N)
W[1]<-0
mu<-0.5
for(i in 2:(N+1)){
  W[i]<-
}
