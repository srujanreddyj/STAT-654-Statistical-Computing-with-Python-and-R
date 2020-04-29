X <- cbind(c(0.57,0.10,0.26,0.85,0.50,0.30),c(0.91,0.37,0.14,0.04,0.30,0.70))
Y <- matrix(c(0,0,0,0.9266,0.6103,0.3631,
              0,0,0.2970,0.8524,0,0,
              0,0,0,0,0.2995,0.5656,
              0,0,0,0,0,0,
              0,0,0,0,0,0,
              0,0,0,0,0,0),nrow=6,byrow=TRUE)
Y <- Y + t(Y)
Y
## create the W matrix
W <- 1*(Y!=0)
W
Xk <- X[5:6,]


## this is L(x_1,x_2,x_3,x_4) above Equation 12 in paper
## note there is a typo in the paper, first exp should be
## to the W_ij power. corrected in this code
likelihood <- function(theta,Xk,Y){
  thetaX <- rbind(matrix(theta,nrow=length(theta)/2,byrow=TRUE),Xk) # columns filled first, Xk appended
  Xd <- as.matrix(dist(thetaX))
  temp <- -Xd^2 / (2*.3^2) # Log prob. observe distance
  out <- W*(-(Y - Xd)^2 / (2*0.02^2) + temp) + (1-W)*log(1-exp(temp))
  return(sum(out[lower.tri(out)]))
}


## we run the optimizer N times from
## random starting locations
N <- 500
soln <- array(0,dim=c(N,4,2))

maxit <- 1000
for(ii in 1:N){
  init <- runif(8,min=-.5,max=1.5)
  a <- optim(init,likelihood,method="BFGS",Xk=Xk,Y=Y,control=list(maxit=maxit,fnscale=-1))
  soln[ii,,] <- matrix(a$par,ncol=2,byrow=TRUE)
}

##Plot the results
plot(0,0,xlim=c(-.5,1.5),ylim=c(-.6,1.5),col=0,xlab="",ylab="")
pch_est <- c(0,1,2,5)
## make colors 1:4 transparent
cols <- rgb(t(col2rgb(1:4)),alpha=70,maxColorValue=255)
for(jj in 1:4){
  points(soln[,jj,1],soln[,jj,2],xlab="x_11",ylab="x_12",col=cols[jj],pch=pch_est[jj])
  points(X[jj,1],X[jj,2],col=jj,cex=2,pch=as.character(jj),font=2)
}
legend("topright",paste0("Sensor ",1:4),col=1:4,pch=pch_est)
