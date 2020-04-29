rm(list=ls())

home <- "D:/Texas A&M/Semester 2/STAT654/Project/Exoplanet/Exoplanet/"
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
 
##RAM Function 

ram.kernel <- function(current.location, current.aux, loc.number, scale) {
  
  #eps <- 10^(-308)
  eps <- 10^(-308)
  accept <- 0 
  x.c <- current.location 
  x.c.den <- post_den(x.c)
  log.x.c.den <- log(x.c.den)
  z.c <- current.aux
  z.c.den <- post_den(z.c) 
  log.z.c.den <- log(z.c.den)
  
  # downhill 
  x.p1 <- x.c 
  x.p1[loc.number] <- x.p1[loc.number] + rnorm(1, 0, scale)
  x.p1.den <- post_den(x.p1)
  log.x.p1.den <- log(x.p1.den)
  N.d <- 1
  while (-rexp(1) > log(x.c.den + eps) - log(x.p1.den + eps)) {
    x.p1 <- x.c 
    x.p1[loc.number] <- x.p1[loc.number] + rnorm(1, 0, scale)
    x.p1.den <- post_den(x.p1)
    log.x.p1.den <- log(x.p1.den)
    N.d <- N.d + 1
  }
  
  # uphill
  x.p2 <- x.p1
  x.p2[loc.number] <- x.p2[loc.number] + rnorm(1, 0, scale)
  x.p2.den <- post_den(x.p2)
  log.x.p2.den <- log(x.p2.den)
  N.u <- 1
  while (-rexp(1) > log(x.p2.den + eps) - log(x.p1.den + eps)) {
    x.p2 <- x.p1
    x.p2[loc.number] <- x.p2[loc.number] + rnorm(1, 0, scale)
    x.p2.den <- post_den(x.p2)
    log.x.p2.den <- log(x.p2.den)
    N.u <- N.u + 1
  }
  
  # downhill for N.d
  N.dz <- 1     # number of total downhill trials for estimate
  z <- x.p2
  z[loc.number] <- z[loc.number] + rnorm(1, 0, scale)
  z.den <- post_den(z)
  log.z.den <- log(z.den)
  while (-rexp(1) > log(x.p2.den + eps) - log(z.den + eps)) {
    z <- x.p2
    z[loc.number] <- z[loc.number] + rnorm(1, 0, scale)
    z.den <- post_den(z)
    log.z.den <- log(z.den)
    N.dz <- N.dz + 1
  }
  
  # accept or reject the proposal
  min.nu <- min(1, (x.c.den + eps) / (z.c.den + eps))
  min.de <- min(1, (x.p2.den + eps) / (z.den + eps))
  l.mh <- log.x.p2.den - log.x.c.den + log(min.nu) - log(min.de)
  
  if (l.mh > -rexp(1)) {
    x.c <- x.p2
    z.c <- z
    accept <- 1
  }
  
  c(x.c, z.c, N.d, N.u, N.dz, accept)
}


MHwG.RAM <- function(initial.loc, initial.aux, jump.scale, 
                    n.sample = 10, n.burn = 10) {
  
  print(Sys.time())
  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 7)
  out <- matrix(NA, nrow = n.total, ncol = 7)
  loc.t <- initial.loc
  aux.t <- initial.aux
  Nd <- matrix(NA, nrow = n.total, ncol = 7)
  Nu <- matrix(NA, nrow = n.total, ncol = 7)
  Nz <- matrix(NA, nrow = n.total, ncol = 7)
  
  for (i in 1 : n.total) {
    for (j in 1 : 7) {
      TEMP <- ram.kernel(loc.t, aux.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 7]
      aux.t <- TEMP[8 : 14]
      Nd[i, j] <- TEMP[15]
      Nu[i, j] <- TEMP[16]
      Nz[i, j] <- TEMP[17]
      accept[i, j] <- TEMP[18]
    }
    out[i, ] <- loc.t
  }
  print(Sys.time())
  list(x = out[-c(1 : n.burn), ], 
       accept = accept[-c(1 : n.burn), ],
       N.d = Nd[-c(1 : n.burn), ],
       N.u = Nu[-c(1 : n.burn), ],
       N.z = Nz[-c(1 : n.burn), ])
  
}

j.scale <- c(10,rep(2,6))   #SD of proposal density for the 7 parameters
a<- c(-1.9,40,2.4,0.12,4.2,0.37,1.2)    #Initial/Starting values of the parameters

#We run RAM for 10000 iterations with burin-in of 2000

planet.ram <- MHwG.RAM(a, a, jump.scale = j.scale, n.sample = 10000, n.burn = 2000)
planet.ram <- planet.ram$x

#Trace plots
par(mfrow=c(2,4))
plot(1:10000,planet.ram[,1],type='l',ylim=c(-1000,1000),xlab = 'Iterations',ylab='Density',main = 'RV Velocity offset')
plot(1:10000,planet.ram[,2],type='l',ylim=c(39,45),xlab = 'Iterations',ylab='Density',main = 'Planets orbital period')
plot(1:10000,planet.ram[,3],type='l',ylim=c(0,5),xlab = 'Iterations',ylab='Density',main = 'RV Semi Amplitude')
plot(1:10000,planet.ram[,4],type='l',ylim=c(0,1),xlab = 'Iterations',ylab='Density',main = 'Planets eccentricity')
plot(1:10000,planet.ram[,5],type='l',xlab = 'Iterations',ylab='Density',main = 'Argument of Epicenter')
plot(1:10000,planet.ram[,6],type='l',xlab = 'Iterations',ylab='Density',main = 'Mean Anomaly')
plot(1:10000,planet.ram[,7],type='l',ylim=c(1,3.5),xlab = 'Iterations',ylab='Density',main = 'White Noise')

#Density Plots
par(mfrow = c(2,4))
hist(planet.ram[,1],30,xlab = 'C', main = 'RV Velocity offset',freq = FALSE)
lines(density(planet.ram[,1]),col='Red')
hist(planet.ram[,2],30,xlab = 'P', main = 'Planets orbital period',freq = FALSE)
lines(density(planet.ram[,2]),col='Red')
hist(planet.ram[,3],30,xlab = 'K', main = 'RV Semi Amplitude',freq = FALSE)
lines(density(planet.ram[,3]),col='Red')
hist(planet.ram[,4],30,xlab = 'e', main = 'Planets eccentricity',freq = FALSE)
lines(density(planet.ram[,4]),col='Red')
hist(planet.ram[,5],30,xlab = 'W', main = 'Argument of Epicenter',freq = FALSE)
lines(density(planet.ram[,5]),col='Red')
hist(planet.ram[,6],30,xlab = 'M', main = 'Mean Anomaly',freq = FALSE)
lines(density(planet.ram[,6]),col='Red')
hist(planet.ram[,7],30,xlab = 'SigmaJ', main = 'White Noise',freq = FALSE)
lines(density(planet.ram[,7]),col='Red')

