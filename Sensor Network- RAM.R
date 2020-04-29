######## Data
 
# Observation indicators from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Ob <- matrix(c(1, 0, 1, 0, 1, 0, 1, 0), ncol = 2)

# Observation indicators among the first four sensors. 
Os <- matrix(c(0, 0, 0, 1,
               0, 0, 1, 1,
               0, 1, 0, 0,
               1, 1, 0, 0), ncol = 4)

# Each row indicates the location of the known sensors (5th and 6th).
Xb <- matrix(c(0.5, 0.3, 0.3, 0.7), ncol = 2)


# Each row indicates the location of the unknown sensors (1st, 2nd, 3rd, and 4th).
Xs <- matrix(c(0.5748, 0.0991, 0.2578, 0.8546, 
               0.9069, 0.3651, 0.1350, 0.0392), ncol = 2)

# The observed distances from the fifth sensor (1st column) to the first four sensors
# and those from the sixth sensor (2nd column) to the first four sensors.
Yb <- matrix(c(0.6103, 0, 0.2995, 0, 
               0.3631, 0, 0.5656, 0), ncol = 2)

# Observed distances among the first four sensors.
Ys <- matrix(c(0, 0, 0, 0.9266,
               0, 0, 0.2970, 0.8524,
               0, 0.2970, 0, 0,
               0.9266, 0.8524, 0, 0), ncol = 4)

######## Target joint posterior density
 
norm2 <- function(loca, locb) {
  sqrt(sum((loca - locb)^2))
}

l.target <- function(loc, R = 0.3, sigma = 0.02, Ob, Os, Xb, Xs, Yb, Ys) {
  
  First.term <- NULL
  for (i in 1 : 2) {
    TEMP <- sapply(1 : 4, function(j) {
      exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Ob[j, i]) *
        (1 - exp(-norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Ob[j, i]) 
    })
    First.term <- c(First.term, TEMP)
  }  
   
  Second.term <- NULL
  for (i in 1 : 3) {
    TEMP <- sapply((i + 1) : 4, function(j) {
      exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                 loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2 * Os[i, j]) *
        (1 - exp(-norm2(loc[(2 * i -1) : (2 * i)], 
                        loc[(2 * j -1) : (2 * j)])^2 / 2 / R^2))^(1 - Os[i, j])  
    })
    Second.term <- c(Second.term, TEMP)
  }
   
  First.obs.term <- NULL
  for (i in 1 : 2) {
    TEMP <- sapply(1 : 4, function(j) {
      dnorm(Yb[j, i], mean = norm2(Xb[i, ], loc[(2 * j -1) : (2 * j)]), 
            sd = sigma)^Ob[j, i]
    })
    First.obs.term <- c(First.obs.term, TEMP)
  }
  
  Second.obs.term <- NULL
  for (i in 1 : 3) {
    TEMP <- sapply((i + 1) : 4, function(j) {
      dnorm(Ys[i, j], mean = norm2(loc[(2 * i -1) : (2 * i)], 
                                   loc[(2 * j -1) : (2 * j)]), 
            sd = sigma)^Os[i, j]
    })
    Second.obs.term <- c(Second.obs.term, TEMP)
  }  
  
  log.lik <- sum(log(c(First.term, Second.term, First.obs.term, Second.obs.term)))
  post <- log.lik + sum(dnorm(loc, mean = rep(0, 8), sd = rep(10, 8), log = TRUE))
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
  accept <- matrix(0, nrow = n.total, ncol = 4)
  out <- matrix(NA, nrow = n.total, ncol = 8)
  loc.t <- initial.loc
  aux.t <- initial.aux
  Nd <- matrix(NA, nrow = n.total, ncol = 4)
  Nu <- matrix(NA, nrow = n.total, ncol = 4)
  Nz <- matrix(NA, nrow = n.total, ncol = 4)
  
  for (i in 1 : n.total) {
    for (j in 1 : 4) {
      TEMP <- ram.kernel(loc.t, aux.t, j, jump.scale[j])
      loc.t <- TEMP[1 : 8]
      aux.t <- TEMP[9 : 16]
      Nd[i, j] <- TEMP[17]
      Nu[i, j] <- TEMP[18]
      Nz[i, j] <- TEMP[19]
      accept[i, j] <- TEMP[20]
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

j.scale <- rep(1.08, 4)   #SD of the proposal density

#We run RAM for 500000 iterations with burn-in period 200000

res.ram <- MHwG.RAM(runif(8), runif(8), jump.scale = j.scale, 
                                Ob, Os, Xb, Xs, Yb, Ys, 
                                n.sample = 500000, n.burn = 200000)

res.ram_6sensor <- res.ram$x

########Summary- No. of proposals at each step

colMeans(res.ram$accept)
colMeans(res.ram$N.d)
colMeans(res.ram$N.u)
colMeans(res.ram$N.z)

########Trace plots
par(mfrow=c(2,2))
plot(1:500000,res.ram_6sensor[,1],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[11])))
plot(1:500000,res.ram_6sensor[,2],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[12])))
plot(1:500000,res.ram_6sensor[,3],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[21])))
plot(1:500000,res.ram_6sensor[,4],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[22])))
mtext("Trace Plot", outer = TRUE, cex = 1.5, line=-2)

plot(1:500000,res.ram_6sensor[,5],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[31])))
plot(1:500000,res.ram_6sensor[,6],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[32])))
plot(1:500000,res.ram_6sensor[,7],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[41])))
plot(1:500000,res.ram_6sensor[,8],type='l',xlab = 'Iterations',ylab='Density',main = expression(paste(X[42])))
mtext("Trace Plot", outer = TRUE, cex = 1.5, line=-2)

# Scatterplots & comparison with BFGS

par(mfrow = c(2,2))

plot(res.ram_6sensor[, c(1, 2)], pch = 46, xlim = c(-0.2, 0.8), ylim = c(0.3, 1.1),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[1]))))
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[12])), line = 1.9, cex = 1.2)
abline(v = Xs[1,1], lty = 2, lwd = 1)
abline(h = Xs[1,2], lty = 2, lwd = 1)
#points(soln[,1,1],soln[,1,2],col='Red',pch=8) #BFGS Solution

plot(res.ram_6sensor[, c(3, 4)], pch = 46, xlim = c(-0.3, 1.3), ylim = c(-0.5, 1.1),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[2]))))
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[22])), line = 1.9, cex = 1.2)
abline(v = Xs[2,1], lty = 2, lwd = 1)
abline(h = Xs[2,2], lty = 2, lwd = 1)
#points(soln[,2,1],soln[,2,2],col='Red',pch=8)  #BFGS Solution

plot(res.ram_6sensor[, c(5, 6)], pch = 46, xlim = c(0, 1), ylim = c(-0.1, 0.65),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[3]))))
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[32])), line = 1.9, cex = 1.2)
abline(v = Xs[3,1], lty = 2, lwd = 1)
abline(h = Xs[3,2], lty = 2, lwd = 1)
#points(soln[,3,1],soln[,3,2],col='Red',pch=8)  #BFGS Solution

plot(res.ram_6sensor[, c(7, 8)], pch = 46, xlim = c(-1, 2), ylim = c(-0.5, 2),
     xlab = "", ylab = "", main = "")
title(expression(bold(paste("RAM: ", x[4]))))
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = expression(bold(x[42])), line = 1.9, cex = 1.2)
abline(v = Xs[4,1], lty = 2, lwd = 1)
abline(h = Xs[4,2], lty = 2, lwd = 1) 
#points(soln[,4,1],soln[,4,2],col='Red',pch=8)  #BFGS Solution

# Density plots
par(mfrow = c(2,4))
hist(res.ram_6sensor[, 1],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[11])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[11]))))
lines(density(res.ram_6sensor[, 1]), lwd = 1.5)
abline(v = Xs[1,1], lty = 2, lwd = 2,col ='Red')

hist(res.ram_6sensor[, 2],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[12])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[12]))))
lines(density(res.ram_6sensor[, 2]), lwd = 1.5)
abline(v = Xs[1,2], lty = 2, lwd = 2,col ='Red')

hist(res.ram_6sensor[, 3],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[21])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[21]))))
lines(density(res.ram_6sensor[, 3]), lwd = 1.5)
abline(v = Xs[2,1], lty = 2, lwd = 2,col ='Red')

hist(res.ram_6sensor[, 4],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[22])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[22]))))
lines(density(res.ram_6sensor[, 4]), lwd = 1.5)
abline(v = Xs[2,2], lty = 2, lwd = 2,col ='Red')

hist(res.ram_6sensor[, 5],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[31])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[31]))))
lines(density(res.ram_6sensor[, 5]), lwd = 1.5)
abline(v = Xs[3,1], lty = 2, lwd = 2,col ='Red')

hist(res.ram_6sensor[, 6],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[32])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[32]))))
lines(density(res.ram_6sensor[, 6]), lwd = 1.5)
abline(v = Xs[3,2], lty = 2, lwd = 2,col ='Red')

hist(res.ram_6sensor[, 7],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[41])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[41]))))
lines(density(res.ram_6sensor[, 7]), lwd = 1.5)
abline(v = Xs[4,1], lty = 2, lwd = 2,col ='Red')

hist(res.ram_6sensor[, 8],30,freq = FALSE, main = "", xlab = "", ylab = "")
mtext(side = 1, text = expression(bold(x[42])), line = 1.6, cex = 1.2)
mtext(side = 2, text = "Density", line = 1.7, cex = 1.2, las = 0)
title(expression(bold(paste("RAM: ", x[42]))))
lines(density(res.ram_6sensor[, 8]), lwd = 1.5)
abline(v = Xs[4,2], lty = 2, lwd = 2,col ='Red')
