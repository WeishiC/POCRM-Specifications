###############################################################
###############################################################
##                                                           ##
##  Simulations under the BLRM design (Riviere et al, 2014)  ##
##  Last update: 10 Jun 2026                                 ##
##                                                           ##
###############################################################
###############################################################

library(rjags); library(coda)

#
#  Implement 
#
blrm.imp <- function(y.data=NULL, n.data=NULL, comb.curr, TTL, a, b, cohortsize, ncohort, 
                     c_e, c_d, n.adapt=1e+3){
  which.is.max <- function (x) 
  {
    y <- seq_along(x)[x == max(x)]
    if (length(y) > 1L) 
      sample(y, 1L)
    else y
  }
  #n.combo <- length(y.data)
  J <- length(a)
  K <- length(b)
  
  jags.script <- "
  model {
    # Likelihood
    for(j in 1:J) {
      for(k in 1:K) {
        p[j,k] <- 1/(1 + exp(-beta[1] - beta[2]*a[j] - beta[3]*b[k] - beta[4]*a[j]*b[k]))
        y[j,k] ~ dbin(p[j,k], n[j,k])
      }
    }

    # Prior
    beta[1] ~ dnorm(0, 1/10)
    beta[2] ~ dexp(1)
    beta[3] ~ dexp(1)
    beta[4] ~ dnorm(0, 1/10)
    
    # DLT probability
    for(j in 1:J) {
      for(k in 1:K) {
        # DLT probabilities
        p.hat[j,k] <- 1/(1 + exp(-beta[1] - beta[2]*a[j] - beta[3]*b[k] - beta[4]*a[j]*b[k]))
        # de-escalation prob's
        p.d[j,k] <- step(p.hat[j,k]-TTL)
        # escalation prob's
        p.e[j,k] <- step(TTL-p.hat[j,k])
      }
    }
  }"
  
  # estimation
  model.fit <- jags.model(textConnection(jags.script), quiet = TRUE,
                          data=list(y=y.data, n=n.data, a=a, b=b, TTL=TTL, J=J, K=K),
                          n.chains=1, n.adapt=n.adapt)
  update(model.fit, n.adapt, progress.bar="none")
  tt<-jags.samples(model.fit, c('beta', 'p.hat', 'p.e', 'p.d'), n.adapt, progress.bar="none")
  p.hat <- p.e <- p.d <- matrix(nrow=J, ncol=K)
  for (j in 1:J) {
    for (k in 1:K) {
      p.hat[j,k] <- mean(tt$p.hat[j,k,,])
      p.e[j,k] <- mean(tt$p.e[j,k,,])
      p.d[j,k] <- mean(tt$p.d[j,k,,])
    }
  }
  
  
  # escalation rules
  acceptable <- NULL
  # escalate
  if(p.e[comb.curr[1], comb.curr[2]]>c_e) {
    if((comb.curr[1]==J) & (comb.curr[2]==K)) {
      comb.next <- comb.curr
    } else {
      if(comb.curr[1]<J) acceptable <- rbind(acceptable, c(comb.curr[1]+1, comb.curr[2]))
      if(comb.curr[2]<K) acceptable <- rbind(acceptable, c(comb.curr[1], comb.curr[2]+1))
      if((comb.curr[1]<J) & (comb.curr[2]>1)) acceptable <- rbind(acceptable, c(comb.curr[1]+1, comb.curr[2]-1))
      if((comb.curr[1]>1) & (comb.curr[2]<K)) acceptable <- rbind(acceptable, c(comb.curr[1]-1, comb.curr[2]+1))
      
      p.acceptable <- apply(acceptable, 1, function(u) p.hat[u[1], u[2]])
      dist.e <- Vectorize(function(u) {
        out <- u-TTL
        if(out<0) out <- 10
        out
      }, "u")
      comb.next <- acceptable[which.min(dist.e(p.acceptable)),]
    }
  }
  # de-escalate
  if(p.d[comb.curr[1], comb.curr[2]]>c_d) {
    if((comb.curr[1]==1) & (comb.curr[2]==1)) {
      comb.next <- comb.curr
    } else {
      if(comb.curr[1]>1) acceptable <- rbind(acceptable, c(comb.curr[1]-1, comb.curr[2]))
      if(comb.curr[2]>1) acceptable <- rbind(acceptable, c(comb.curr[1], comb.curr[2]-1))
      if((comb.curr[1]<J) & (comb.curr[2]>1)) acceptable <- rbind(acceptable, c(comb.curr[1]+1, comb.curr[2]-1))
      if((comb.curr[1]>1) & (comb.curr[2]<K)) acceptable <- rbind(acceptable, c(comb.curr[1]-1, comb.curr[2]+1))
      
      p.acceptable <- apply(acceptable, 1, function(u) p.hat[u[1], u[2]])
      dist.d <- Vectorize(function(u) {
        out <- TTL-u
        if(out<0) out <- 10
        out
      }, "u")
      comb.next <- acceptable[which.min(dist.d(p.acceptable)),]
    }
  }
  # stay
  if((p.e[comb.curr[1], comb.curr[2]]<=c_e) & (p.d[comb.curr[1], comb.curr[2]]<=c_d)) comb.next <- comb.curr
  
  return(list(comb.next=comb.next, tox.data=y.data, n.tox=sum(y.data), pt.allocation=n.data, 
              ptox.hat=p.hat, p.e=p.e, p.d=p.d))
}


#
# Simulate 1 trial 
#
blrm<-function(p0, TTL, c_e, c_d, a, b, cohortsize, ncohort, start.comb, n.adapt=1e+3){
  which.is.max <- function (x) 
  {
    y <- seq_along(x)[x == max(x)]
    if (length(y) > 1L) 
      sample(y, 1L)
    else y
  }
  J <- length(a)
  K <- length(b)
  
  y <- n <- comb.select <- matrix(0, nrow=J, ncol=K) 
  y.tox <- numeric(ncohort)
  x <- NULL
  
  # start the trial here
  i <- 1
  comb.curr <- start.comb
  
  x <- rbind(x, comb.curr)
  while (i <= ncohort) {
    y.curr <- rbinom(1, cohortsize, p0[comb.curr[1], comb.curr[2]])
    y.tox[i] <- y.curr
    y[comb.curr[1], comb.curr[2]] <- y[comb.curr[1], comb.curr[2]] + y.curr
    n[comb.curr[1], comb.curr[2]] <- n[comb.curr[1], comb.curr[2]] + cohortsize
    
    # implement BLRM
    imp <- blrm.imp(y.data=y, n.data=n, comb.curr=comb.curr, TTL=TTL, a=a, b=b,
                    cohortsize=cohortsize, ncohort=ncohort, c_e=c_e, c_d=c_d, n.adapt=n.adapt)
    
    comb.curr <- imp$comb.next
    x <- rbind(x, comb.curr)
    i <- i + 1
  }
  
  comb.select[comb.curr[1], comb.curr[2]] <- comb.select[comb.curr[1], comb.curr[2]] + 1
  return(list(comb.select=comb.select, tox.data=y, n.tox=sum(y), pt.allocation=n, 
              ptox.hat=imp$ptox.hat, x=x, y.tox=y.tox, p.e=imp$p.e, p.d=imp$p.d))
}

#
#  Simulate multiple trials
#
BLRM.sim <- function(n.sim, p0, TTL, c_e, c_d, a, b, cohortsize, ncohort, start.comb, n.adapt=1e+3) {
  sim <- replicate(n.sim, blrm(p0, TTL, c_e, c_d, a, b, cohortsize, ncohort, start.comb, n.adapt))
  J <- length(a); K <- length(b)
  comb.select_mat <- array(dim=c(n.sim, J, K))
  n.tox_vec <- rep()
  for (r in 1:n.sim) {
    comb.select_mat[r,,] <- sim[,r]$comb.select
    n.tox_vec[r] <- sim[,r]$n.tox
  }
  return(list(comb.select = apply(comb.select_mat, c(2,3), mean),
              comb.select.full = comb.select_mat,
              n.tox = mean(n.tox_vec)))
}

# design parameters
TTL <- 0.3
skeleton <- c(0.2, 0.3, 0.4)
a <- b <- log(skeleton/(1-skeleton))
cohortsize <- 3
ncohort <- 20
c_e <- 0.85
c_d <- 0.45
start.comb <- c(1,1)

# run 1000 simulations
p0 <- matrix(c(0.25, 0.30, 0.50, 0.35, 0.45, 0.55, 0.40, 0.60, 0.65), nrow=length(a))
sim <- BLRM.sim(n.sim=1000, p0, TTL, c_e, c_d, a, b, cohortsize, ncohort, start.comb)
