#Generating simulated points

set.seed(6)
X<- matrix(round(runif(18),2),ncol = 2)
X
par(mfrow=c(1,1))
plot(X[,1],X[,2],pch=as.character(1:9),col=c(rep(1,6),rep(2,3)), xlim = c(0.15, 1.1), 
     ylim = c(-0.2, 1.1),xlab="x-coordinate",ylab="y-coordinate",main='Sensor network- 9 nodes')
dist(X)

distance<- matrix(0,nrow = 9,ncol = 9)
distance

norm2 <- function(loca, locb) {
  sqrt(sum((loca - locb)^2))
}

for (i in 1:9){
  for (j in 1:9){
    distance[i,j]=norm2(X[i,],X[j,])
  }
}

#Add random error
obs_distance = distance + matrix(rnorm(81,0,0.02),nrow = 9)

for (i in 1:9){
  obs_distance[i:9,i]=0
}
obs_distance

#Now we randomly remove some distance by randomly selecting 0 and 1
set.seed(1)
O <- matrix(rbinom(81,1,0.5),nrow=9)
O
for (i in 1:9){
  O[i:9,i]=0
}
O
O <- O + t(O)
O
obs_distance
Y = O*obs_distance

Y <- Y + t(Y)
Y

##Distance (binary) of 7th,8th and 9th to 1-6th 
Ob<-O[1:6,7:9]
Ob

Os<-O[1:6,1:6]
Os
#The true locations of 7th-9th sensor
Xb <- X[7:9,]
Xb
#The location of unknown sensors (1-6)
Xs <- X[1:6,]
Xs

#Distance of 7,8,9 to 1st 6 points and within 1st 6 points
Yb <- Y[1:6,7:9]
Yb
Ys <- Y[1:6,1:6]
Ys

######## Target joint posterior density
l.target <- function(loc, R = 0.3, sigma = 0.02, Ob, Os, Xb, Xs, Yb, Ys) {
  
  First.term <- NULL
  for (i in 1 : 3) {
    TEMP <- sapply(1 : 6, function(j) {
      exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Ob[j, i]) *
        (1 - exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Ob[j, i]) 
    })
    First.term <- c(First.term, TEMP)
  }  
  
  Second.term <- NULL
  for (i in 1 : 5) {
    TEMP <- sapply((i + 1) : 6, function(j) {
      exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                 loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Os[i, j]) *
        (1 - exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                        loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Os[i, j])  
    })
    Second.term <- c(Second.term, TEMP)
  }
  
  First.obs.term <- NULL
  for (i in 1 : 3) {
    TEMP <- sapply(1 : 6, function(j) {
      dnorm(Yb[j, i], mean = norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)]), 
            sd = sigma)^Ob[j, i]
    })
    First.obs.term <- c(First.obs.term, TEMP)
  }
  
  Second.obs.term <- NULL
  for (i in 1 : 5) {
    TEMP <- sapply((i + 1) : 6, function(j) {
      dnorm(Ys[i, j], mean = norm2(loc[(2 * i -1) : (2 * i)], 
                                   loc[(2 * j -1) : (2 * j)]), 
            sd = sigma)^Os[i, j]
    })
    Second.obs.term <- c(Second.obs.term, TEMP)
  }  
  
  log.lik <- sum(log(c(First.term, Second.term, First.obs.term, Second.obs.term)))
  post <- log.lik + sum(dnorm(loc, mean = rep(0, 12), sd = rep(10, 12), log = TRUE))
  post
  
}


######## RAM

ram.kernel <- function(current.location, current.aux, loc.number, scale) {
  
  eps <- 10^(-308)
  accept <- 0 
  x.c <- current.location 
  log.x.c.den <- l.target(x.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.c.den <- exp(log.x.c.den)
  z.c <- current.aux
  log.z.c.den <- l.target(z.c, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.c.den <- exp(log.z.c.den)
  
  # downhill 
  x.p1 <- x.c
  x.p1[(2 * loc.number - 1) : (2 * loc.number)] <- x.p1[(2 * loc.number - 1) : (2 * loc.number)] + 
    rnorm(2, 0, scale)
  log.x.p1.den <- l.target(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p1.den <- exp(log.x.p1.den)
  N.d <- 1
  while (-rexp(1) > log(x.c.den + eps) - log(x.p1.den + eps)) {
    x.p1 <- x.c 
    x.p1[(2 * loc.number - 1) : (2 * loc.number)] <- x.p1[(2 * loc.number - 1) : (2 * loc.number)] + 
      rnorm(2, 0, scale)
    log.x.p1.den <- l.target(x.p1, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p1.den <- exp(log.x.p1.den)
    N.d <- N.d + 1
  }
  
  # uphill
  x.p2 <- x.p1
  x.p2[(2 * loc.number - 1) : (2 * loc.number)] <- x.p2[(2 * loc.number - 1) : (2 * loc.number)] + 
    rnorm(2, 0, scale)
  log.x.p2.den <- l.target(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  x.p2.den <- exp(log.x.p2.den)
  N.u <- 1
  while (-rexp(1) > log(x.p2.den + eps) - log(x.p1.den + eps)) {
    x.p2 <- x.p1
    x.p2[(2 * loc.number - 1) : (2 * loc.number)] <- x.p2[(2 * loc.number - 1) : (2 * loc.number)] + 
      rnorm(2, 0, scale)
    log.x.p2.den <- l.target(x.p2, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    x.p2.den <- exp(log.x.p2.den)
    N.u <- N.u + 1
  }
  
  # downhill for N.d
  N.dz <- 1     # number of total downhill trials for estimate
  z <- x.p2
  z[(2 * loc.number - 1) : (2 * loc.number)] <- z[(2 * loc.number - 1) : (2 * loc.number)] + 
    rnorm(2, 0, scale)
  log.z.den <- l.target(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
  z.den <- exp(log.z.den)
  while (-rexp(1) > log(x.p2.den + eps) - log(z.den + eps)) {
    z <- x.p2
    z[(2 * loc.number - 1) : (2 * loc.number)] <- z[(2 * loc.number - 1) : (2 * loc.number)] + 
      rnorm(2, 0, scale)
    log.z.den <- l.target(z, 0.3, 0.02, Ob, Os, Xb, Xs, Yb, Ys)
    z.den <- exp(log.z.den)
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
                     Ob, Os, Xb, Xs, Yb, Ys, n.sample = 10, n.burn = 10) {
  
  print(Sys.time())
  n.total <- n.sample + n.burn
  accept <- matrix(0, nrow = n.total, ncol = 6)
  out <- matrix(NA, nrow = n.total, ncol = 12)
  loc.t <- initial.loc
  aux.t <- initial.aux
  Nd <- matrix(NA, nrow = n.total, ncol = 6)
  Nu <- matrix(NA, nrow = n.total, ncol = 6)
  Nz <- matrix(NA, nrow = n.total, ncol = 6)
  
  for (i in 1 : n.total) {
    for (j in 1 : 6) {
      TEMP <- ram.kernel(loc.t, aux.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 12]
      aux.t <- TEMP[13 : 24]
      Nd[i, j] <- TEMP[25]
      Nu[i, j] <- TEMP[26]
      Nz[i, j] <- TEMP[27]
      accept[i, j] <- TEMP[28]
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

j.scale <- rep(1.08, 6)  #SD for proposal density

#We run RAM for 250000 iterations with burin-in of 50000

res.ram <- MHwG.RAM(runif(12), runif(12), jump.scale = j.scale, 
                    Ob, Os, Xb, Xs, Yb, Ys, 
                    n.sample = 250000, n.burn = 50000)

res.ram_9sensor <- res.ram$x

########Trace plots
par(mfrow=c(2,3))
plot(1:250000,res.ram_9sensor[,1],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[11])))
plot(1:250000,res.ram_9sensor[,2],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[12])))
plot(1:250000,res.ram_9sensor[,3],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[21])))
plot(1:250000,res.ram_9sensor[,4],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[22])))
plot(1:250000,res.ram_9sensor[,5],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[31])))
plot(1:250000,res.ram_9sensor[,6],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[32])))
mtext("Trace Plot", outer = TRUE, cex = 1.5, line=-2)

par(mfrow=c(2,3))
plot(1:250000,res.ram_9sensor[,7],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[41])))
plot(1:250000,res.ram_9sensor[,8],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[42])))
plot(1:250000,res.ram_9sensor[,9],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[51])))
plot(1:250000,res.ram_9sensor[,10],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[52])))
plot(1:250000,res.ram_9sensor[,11],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[61])))
plot(1:250000,res.ram_9sensor[,12],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[62])))
mtext("Trace Plot", outer = TRUE, cex = 1.5, line=-2)

# Scatterplots

par(mfrow = c(2,3))

plot(res.ram_9sensor[, c(1, 2)], pch = 46, xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[1]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[12])), line = 1.9, cex = 1.2)
abline(v = Xs[1,1], lty = 2, lwd = 1)
abline(h = Xs[1,2], lty = 2, lwd = 1)

plot(res.ram_9sensor[, c(3, 4)], pch = 46, xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[2]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[22])), line = 1.9, cex = 1.2)
abline(v = Xs[2,1], lty = 2, lwd = 1)
abline(h = Xs[2,2], lty = 2, lwd = 1)

plot(res.ram_9sensor[, c(5, 6)], pch = 46, xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[3]))))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[32])), line = 1.9, cex = 1.2)
abline(v = Xs[3,1], lty = 2, lwd = 1)
abline(h = Xs[3,2], lty = 2, lwd = 1)

plot(res.ram_9sensor[, c(7, 8)], pch = 46, xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[4]))))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[42])), line = 1.9, cex = 1.2)
abline(v = Xs[4,1], lty = 2, lwd = 1)
abline(h = Xs[4,2], lty = 2, lwd = 1) 

plot(res.ram_9sensor[, c(9, 10)], pch = 46, xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[5]))))
mtext(side = 1, text = expression(bold(x[51])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[52])), line = 1.9, cex = 1.2)
abline(v = Xs[5,1], lty = 2, lwd = 1)
abline(h = Xs[5,2], lty = 2, lwd = 1) 

plot(res.ram_9sensor[, c(11, 12)], pch = 46, xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[6]))))
mtext(side = 1, text = expression(bold(x[61])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[62])), line = 1.9, cex = 1.2)
abline(v = Xs[6,1], lty = 2, lwd = 1)
abline(h = Xs[6,2], lty = 2, lwd = 1) 

#Histogram & Density plots

par(mfrow = c(3,4))
hist(res.ram_9sensor[, 1],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[11]))))
lines(density(res.ram_9sensor[, 1]), lwd = 1.5)
abline(v = Xs[1,1], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 2],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[12])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[12]))))
lines(density(res.ram_9sensor[, 2]), lwd = 1.5)
abline(v = Xs[1,2], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 3],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[21]))))
lines(density(res.ram_9sensor[, 3]), lwd = 1.5)
abline(v = Xs[2,1], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 4],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[22])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[22]))))
lines(density(res.ram_9sensor[, 4]), lwd = 1.5)
abline(v = Xs[2,2], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 5],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[31]))))
lines(density(res.ram_9sensor[, 5]), lwd = 1.5)
abline(v = Xs[3,1], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 6],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[32])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[32]))))
lines(density(res.ram_9sensor[, 6]), lwd = 1.5)
abline(v = Xs[3,2], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 7],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[41]))))
lines(density(res.ram_9sensor[, 7]), lwd = 1.5)
abline(v = Xs[4,1], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 8],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[42])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[42]))))
lines(density(res.ram_9sensor[, 8]), lwd = 1.5)
abline(v = Xs[4,2], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 9],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[51])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[51]))))
lines(density(res.ram_9sensor[, 9]), lwd = 1.5)
abline(v = Xs[5,1], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 10],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[52])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[52]))))
lines(density(res.ram_9sensor[, 10]), lwd = 1.5)
abline(v = Xs[5,2], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 11],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[61])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[61]))))
lines(density(res.ram_9sensor[, 11]), lwd = 1.5)
abline(v = Xs[6,1], lty = 2, lwd = 2)

hist(res.ram_9sensor[, 12],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[62])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[62]))))
lines(density(res.ram_9sensor[, 12]), lwd = 1.5)
abline(v = Xs[6,2], lty = 2, lwd = 2)
