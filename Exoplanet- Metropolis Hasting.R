rm(list=ls())

home <- "D:/Texas A&M/Semester 2/STAT654/Project/Exoplanet/Exoplanet/" ## please modify this data home
data_home <- paste(home,"data/",sep="") ## please modify this data home
code_home <- paste(home,"R model code/",sep="") ## please modify this code home

# Get posterior (including data)
setwd(code_home)
source("log_post.R")
source("cov_fun.R")
source("cov_make.R")
source("planet_model.R") # Added in sqrt to phi

setwd(data_home)
dataset_num <- 1
data_now <- read.table(paste("rvs_000",dataset_num,".txt",sep=""))
colnames(data_now) <- c("time","rv","sd")
head(data_now)

prior_bounds <- read.table(paste("prior_bounds_000",dataset_num,".txt",sep=""),sep=",")
colnames(prior_bounds) <- c("parameter","planet","min","max")
prior_bounds
setwd(code_home)

source('restricted_priors.R') ## restricted priors
# Given values
tau <- 20 # days (stellar rotation period)
alpha <- sqrt(3) # m/s
lambda_e <- 50 # days
lambda_p <- 0.5 # unitless

Dim <- 7
num_planets <- 1

post_den <- function(x,log=F){
  m <- length(x)/Dim
  x <- matrix(x,m,Dim)
  value <- numeric(m)
  for (i in 1:m){
    planet_paras <- list()
    planet_paras[[1]] <- list(tau=x[i,2],K=x[i,3],e=x[i,4],w=x[i,5],M0=x[i,6],gamma=0)
    for_post <- list(t=data_now[,1],y=data_now[,2],sd=data_now[,3],C=x[i,1],sigmaJ=x[i,7],
                     cov_paras=c(log(tau),log(alpha),log(lambda_e),log(lambda_p)),planet_paras=planet_paras)
    value[i] <- log_post(for_post)
  }
  if (log==F){
    value <- exp(value)
  }
  return(value)
}

##Metropolis Hasting
Niter = 100
chain = matrix(0,nrow = Niter,ncol=7)
chain[1,] = c(-1.9,40,2.4,0.12,4.2,0.37,1.2)


metropolis <- function(current.location, loc.number,scale ){
  x.c <- current.location 
  x.c.den <- post_den(x.c)
  x.p1 <- x.c 
  x.p1[loc.number] <- x.p1[loc.number] + rnorm(1, 0, scale)
  x.p1.den <- post_den(x.p1)
  if ((x.p1.den / x.c.den) > runif(1)){
    x.c <- x.p1
  }
  return(x.c)
}


MHwG.metropolis <- function(initial.loc, jump.scale, 
                     n.sample = 10, n.burn = 10) {
  
  print(Sys.time())
  n.total <- n.sample + n.burn
  out <- matrix(NA, nrow = n.total, ncol = 7)
  loc.t <- initial.loc
  
  for (i in 1 : n.total) {
    for (j in 1 : 7) {
      TEMP <- metropolis(loc.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 7]
    }
    out[i, ] <- loc.t
  }
  print(Sys.time())
  list(x = out[-c(1 : n.burn), ])
}

j.scale <- rep(2,7)
a<- c(-1.9,40,2.4,0.12,4.2,0.37,1.2)

out <- MHwG.metropolis(a,j.scale,n.sample = 10000,n.burn=2000)

###Trace plot
plot(1:10000,out$x[,6],type='l',xlab = 'Iterations',ylab='Density',main = 'Mean Anomaly')