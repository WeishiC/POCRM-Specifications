#######################################################################
#######################################################################
##                                                                   ##
##  Simulation under the surface-free design (Mozgunov et al, 2020)  ##
##  Last update: 10 Jun 2026                                         ##
##                                                                   ##
#######################################################################
#######################################################################
library("rjags")

# Function to Compute Mean Prior Point Estimate of Connections Using the Monotherapy Data
compute.prior.means.SFD <- function(p1, p2) {
  t.prior <- mat.or.vec(length(p1)+length(p2)-1, 1)
  t.prior[1] <- 1  -p1[1] - p2[1] + p1[1]*p2[1]
  for (i in 2:(length(p1))) {
    t.prior[i] <- (1-p1[i])/(1-p1[i-1])
  }
  for (i in 1:(length(p2)-1)) {
    t.prior[length(p2)+i] <- (1-p2[i+1])/(1-p2[i])
  }
  return(t.prior)
}

# implementation function for SFD
# 3x3 grid
SFD.imp <- function(datas, datan, nextdose, a.prior, b.prior, target, iterations=1000) {

  # fit the model using jags
  model1.string <-"
    model {
      s[1] ~ dbin(p[1], n[1])
      p[1] <- 1-theta[1]
      s[2] ~ dbin(p[2], n[2])
      p[2] <- 1-theta[1]*theta[2]
      s[3] ~ dbin(p[3], n[3])
      p[3] <- 1-theta[1]*theta[2]*theta[3]
      s[4] ~ dbin(p[4], n[4])
      p[4] <- 1-theta[1]*theta[4]
      s[5] ~ dbin(p[5], n[5])
      p[5] <- 1-theta[1]*theta[2]*theta[4]
      s[6] ~ dbin(p[6], n[6])
      p[6] <- 1-theta[1]*theta[2]*theta[3]*theta[4]
      s[7] ~ dbin(p[7], n[7])
      p[7] <- 1-theta[1]*theta[4]*theta[5]
      s[8] ~ dbin(p[8], n[8])
      p[8] <- 1-theta[1]*theta[2]*theta[4]*theta[5]
      s[9] ~ dbin(p[9], n[9])
      p[9] <- 1-theta[1]*theta[2]*theta[3]*theta[4]*theta[5]

      theta[1] ~ dbeta(a[1],b[1])T(0,0.999999)
      theta[2] ~ dbeta(a[2],b[2])T(0,0.999999)
      theta[3] ~ dbeta(a[3],b[3])T(0,0.999999)
      theta[4] ~ dbeta(a[4],b[4])T(0,0.999999)
      theta[5] ~ dbeta(a[5],b[5])T(0,0.999999)
    }
  "
  model1.spec <- textConnection(model1.string)

  mydata <- list(s=datas, n=datan, a=a.prior, b=b.prior)
  jags <- jags.model(model1.spec, data =mydata, n.chains=2, n.adapt=iterations, quiet=TRUE)
  update(jags, iterations,progress.bar="none")
  tt <- jags.samples(jags, c('theta'), iterations, progress.bar="none")
  t <- cbind(c(tt$theta[1,,]), c(tt$theta[2,,]), c(tt$theta[3,,]), c(tt$theta[4,,]), c(tt$theta[5,,]))
  t.mean<-colMeans(t)
  # compute toxicity probabilities
  p1 <- 1 - t.mean[1]
  p2 <- 1 - t.mean[1]*t.mean[2]
  p3 <- 1 - t.mean[1]*t.mean[2]*t.mean[3]
  p4 <- 1 - t.mean[1]*t.mean[4]
  p5 <- 1 - t.mean[1]*t.mean[2]*t.mean[4]
  p6 <- 1 - t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]
  p7 <- 1 - t.mean[1]*t.mean[4]*t.mean[5]
  p8 <- 1 - t.mean[1]*t.mean[2]*t.mean[4]*t.mean[5]
  p9 <- 1 - t.mean[1]*t.mean[2]*t.mean[3]*t.mean[4]*t.mean[5]
  # remove not allowed escalations
  if(nextdose==1) p3 <- p5 <- p6 <- p7 <- p8 <- p9 <- 1
  if(nextdose==2) p6 <- p7 <- p8 <- p9 <- 1
  if(nextdose==3) p7 <- p8 <- p9 <- 1
  if(nextdose==4) p6 <- p8 <- p9 <- 1
  if(nextdose==5) p9 <- 1
  if(nextdose==7) p9 <- 1

  p <- c(p1,p2,p3,p4,p5,p6,p7,p8,p9)

  # escalation
  nextdose <- which( abs(p-target) == min(abs(p-target)) )

  return(list(p = p, nextdose=nextdose))
}

# simulate 1 trial
SFD.sim1 <- function(true, a.prior, b.prior, target, n.cohort, cohort=3, start.dose=1, iterations=1000) {
  datad <- mat.or.vec(1,length(true))
  doses.exp <- mat.or.vec(1,length(true))
  doses.tox <- mat.or.vec(1,length(true))
  
  current.dose <- start.dose
  
  for (i in 1:n.cohort) {
    # simulate pseudo data
    outcome <- rbinom(1, cohort, true[current.dose])
    doses.exp[current.dose] <- doses.exp[current.dose]+cohort
    doses.tox[current.dose] <- doses.tox[current.dose]+outcome
    datan <- c(doses.exp)
    datas <- c(doses.tox)
    # nextdose <- start.dose
    
    # implement the SFD
    imp <- SFD.imp(datas=datas, datan=datan, nextdose=current.dose, a.prior=a.prior, b.prior=b.prior, 
                   target=target, iterations=iterations)
    nextdose <- imp$nextdose
    p <- imp$p
    
    # escalation
    current.dose <- nextdose
  }
  d.sel <- rep(0, length(true))
  d.sel[nextdose] <- 1
  
  return(list(d.sel=d.sel, n=datan, n.tox=sum(datas)))
}

# simulate multiple trials
SFD.sim <- function(n.sim, true, a.prior, b.prior, target, n.cohort, cohort=3, start.dose=1, iterations=1000) {
  sim <- replicate(n.sim, SFD.sim1(true=true, a.prior=a.prior, b.prior=b.prior, target=target,
                                   n.cohort=n.cohort, cohort=cohort, start.dose=start.dose, 
                                   iterations=iterations))
  comb.select_mat <- matrix(nrow=n.sim, ncol=length(as.vector(true)))
  n.tox_vec <- c()
  for (r in 1:n.sim) {
    comb.select_mat[r,] <- sim[,r]$d.sel
    n.tox_vec[r] <- sim[,r]$n.tox
  }
  return(list(comb.select = colMeans(comb.select_mat),
              comb.select.full = comb.select_mat,
              n.tox = mean(n.tox_vec)))
}


# Computing the Mean Prior Point Estimates
pA <- c(0.05,0.10,0.20)
pB <- c(0.10,0.20,0.30)
t.prior <- compute.prior.means.SFD(p1=pB,p2=pA)

# Defining the Strength of Prior
c.prior <- rep(4, 5)   # stregth of Beta prior distributions on the connections
a.prior <- t.prior*c.prior # Finding the first parameters of Beta distribution
b.prior <- (1-t.prior)*c.prior # Finding the second parameters of Beta distribution

target<-0.30         # Target Toxicity Level
cohort<-3            # Cohort Size
start.dose<-1        # Starting Dose (corresponds to the combination (1,1))
n<-20                # Number of cohorts in the trial
nsims<-10          # Number of simulations  
iterations<-1000     # Number of samples in MCMC


scen <- matrix(c(0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70,
                 0.25, 0.30, 0.50, 0.35, 0.45, 0.55, 0.40, 0.60, 0.65,
                 0.10, 0.20, 0.30, 0.15, 0.25, 0.35, 0.40, 0.45, 0.50,
                 0.20, 0.25, 0.35, 0.30, 0.40, 0.50, 0.45, 0.55, 0.60,
                 0.15, 0.20, 0.35, 0.25, 0.30, 0.40, 0.45, 0.50, 0.55,
                 0.01, 0.03, 0.05, 0.10, 0.15, 0.30, 0.20, 0.25, 0.35,
                 0.01, 0.10, 0.15, 0.05, 0.20, 0.25, 0.30, 0.35, 0.40,
                 0.05, 0.20, 0.35, 0.10, 0.25, 0.40, 0.15, 0.30, 0.45,
                 0.01, 0.03, 0.07, 0.05, 0.15, 0.20, 0.10, 0.25, 0.30,
                 0.2, 0.3, 0.5, 0.3, 0.45, 0.6, 0.4, 0.55, 0.7,
                 0.1, 0.3, 0.45, 0.2, 0.4, 0.55, 0.3, 0.5, 0.6,
                 0.1, 0.2, 0.3, 0.3, 0.4, 0.55, 0.45, 0.5, 0.6,
                 0.05, 0.2, 0.3, 0.1, 0.3, 0.4, 0.45, 0.5, 0.6,
                 0.05, 0.1, 0.3, 0.15, 0.2, 0.5, 0.3, 0.4, 0.6,
                 0.05, 0.15, 0.3, 0.1, 0.2, 0.4, 0.25, 0.3, 0.5,
                 0.05, 0.2, 0.25, 0.1, 0.3, 0.4, 0.3, 0.5, 0.6,
                 0.05, 0.1, 0.25, 0.15, 0.2, 0.3, 0.3, 0.4, 0.6,
                 0.01, 0.05, 0.25, 0.1, 0.2, 0.3, 0.15, 0.3, 0.5,
                 0.1, 0.2, 0.3, 0.15, 0.3, 0.5, 0.3, 0.4, 0.6), ncol=9, byrow=TRUE)

sim <- SFD.sim(n.sim=1000, true=scen[1,], a.prior, b.prior, target=0.3, n.cohort=20, cohort=3, start.dose = 1)
