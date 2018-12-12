rm(list=ls())
A<-function(x,t,alpha=alpha,theta=theta,a=a,T=T){
  jx<- 4*alpha*sqrt(x*a*T*exp(-alpha*(T-t)))/(sigma^2*(1-exp(-alpha*(T-t))))
  J1<-besselI(x=jx,nu=2*alpha*theta/sigma^2)
  #besselI(x=760,nu=2*alpha*theta/sigma^2)
  #Inf
  J2<-besselI(x=jx, nu=2*alpha*theta/sigma^2-1)
  J1J2<-J1/J2
  if(!is.finite(J1) || !is.finite(J2)){J1J2 <-1} 
  Ax<- alpha*(theta-x) + 2*alpha*x/(1-exp(-alpha*(T-t)))*(sqrt(a*T/x*exp(-alpha*(T-t)))*J1J2-exp(-alpha*(T-t)))
  return(Ax)
}


B<-function(x,sigma=sigma,p=p){
  return(sigma^2*x^p)
}

alpha = 10
theta = 5
sigma = 2
x0 = 2
a = 10
T = 1
q = 2*alpha*theta/sigma^2 - 1

L = 1000
delta = T/L

# condition for the approximation to have strong convergence
check <- sigma^2*q^2/8  >  max(3*alpha,2*sigma)
check
noise <- rnorm(n=L-1,mean=0,sd=sqrt(delta))

y<-array(0,c(L+1))
y[1]<-x0

for (Index in 2:L){
  #    Index<-L
  #    print(c(Index,A(x=y[Index-1],t=(Index-2)*delta,alpha=alpha,theta=theta,a=a,T=T)))
  y[Index] <- abs(y[Index-1] + A(y[Index-1],(Index-2)*delta,alpha=alpha,theta=theta,a=a,T=T)*delta + sqrt(B(y[Index-1],sigma=sigma,p=0.6))*noise[Index-1])
}
y[L+1]<-a*T
tail(y)
steps<-seq(0*delta,L*delta,by=delta)

plot(steps,y,type="l",main="CIR Bridge", xlab="Times", ylab="Trait Value")


##########################
sims<-20
y.array<-array(0,c(sims,L+1))

for(simIndex in 1:sims){
  noise <- rnorm(n=L-1,mean=0,sd=sqrt(delta))
  y<-array(0,c(L+1))
  y[1]<-x0
  for (Index in 2:L){
  #    Index<-L
  #    print(c(Index,A(x=y[Index-1],t=(Index-2)*delta,alpha=alpha,theta=theta,a=a,T=T)))
    y[Index] <- abs(y[Index-1] + A(y[Index-1],(Index-2)*delta,alpha=alpha,theta=theta,a=a,T=T)*delta + sqrt(B(y[Index-1],sigma=sigma,p=0.6))*noise[Index-1])
    }
  y[L+1]<-a*T
  y.array[simIndex,]<-y
  }
steps<-seq(0*delta,L*delta,by=delta)
plot(steps,y.array[1,],type="l",main="CIR Bridge", xlab="Times", ylab="Trait Value")

plot(NA,xlim=c(0,length(y)), ylim=c(min(y.array),max(y.array)),xlab="time step",ylab="trait value", main="CIR Bridge")
for(simIndex in 1:sims){
  lines(1:(L+1),y.array[simIndex,],type = "l")
}
