library(pocrm); library(dfcrm)  # for CRM-related functions
library(nnet)
library(dplyr)
library(combinat)


# two-stage POCRM
pocrm.sim <- function (r, alpha, prior.o, x0, stop, n, theta, nsim, tox.range) {
  sim <- sim1 <- apred <- lik <- pord <- ord <- ahat <- rpred <- next.lev <- n1 <- N <- NULL
  d <- ncol(alpha)                     # number of doses
  s <- nrow(alpha)                     # number of orderings
  nu <- which.min(abs(r - theta))      # true MTC
  orders <- t(apply(alpha, 1, order))  # orderings included
  mod.true <- which(apply(orders, 1, function(v) sum(v==order(r))==length(r)))
  if (nsim > 1) {
    lpocrm <- function(r, alpha, prior.o, x0, stop, n, theta) {
      if (is.vector(alpha)) alpha = t(as.matrix(alpha))  # if only 1 model
      nord.tox <- nrow(alpha)  # =s number of orderings
      mprior.tox <- prior.o    # prior probabilities on the orderings
      bcrml <- function(a, p1, y, n) {
        lik = 0  # log-likelihood
        for (j in 1:length(p1)) {
          lik = lik + y[j] * a * log(p1[j]) + (n[j] - y[j]) * log((1 - p1[j]^a))
        }
        return(lik)
      }
      ncomb <- ncol(alpha)     # =d number of doses
      y <- npts <- ptox.hat <- comb.select <- numeric(ncomb)
      comb.curr <- x0[1]       
      stoprule <- 0
      i <- 1
      # stage1: follow the path defined by x0
      stage1 <- c(x0, rep(ncol(alpha), n - length(x0)))
      while (i <= n) {
        y[comb.curr] <- y[comb.curr] + rbinom(1, 1, r[comb.curr])
        npts[comb.curr] <- npts[comb.curr] + 1
        if (sum(y) == sum(npts)) {
          # if all DLT=1, de-escalate if not already at d1
          comb.curr <- ifelse(comb.curr == 1, comb.curr, comb.curr - 1)
        }
        else if (sum(y) == 0) {
          # if all DLT=0, escalate follow the pre-defined path if not already at d_k
          comb.curr <- ifelse(comb.curr == ncomb, comb.curr, stage1[i + 1])
        }
        else {
          break
        }
        if (any(npts > stop)) {
          stoprule <- 0  # early stop due to the number of patients at any dose exceed *stop*
          break
        }
        i <- i + 1
      }
      
      # stage 2: upon first heterogeneity
      while (sum(npts) <= n) {
        if (sum(y) == 0) {
          stop <- 0 # reach n but still no toxicity
          break
        }
        else {
          like.tox <- est.tox <- rep(0, nord.tox)  # vectors of length S
          mod <- numeric(nord.tox)
          for (k in 1:nord.tox) {
            try <- optim(par=1, fn=bcrml, p1 = alpha[k, ], y = y, n = npts, 
                         method="Brent", lower=0, upper=3, 
                         control=list(maxit=50, fnscale=-1))
            est.tox[k] <- try$par
            if(try$convergence!=0) print("not converge")
            like.tox[k] <- bcrml(a = est.tox[k], p1 = alpha[k, ], y = y, n = npts)
          }
          postprob.tox <- numeric(nord.tox)
          for (i in 1:nord.tox) {
            postprob.tox[i] <- mprior.tox[i] / (mprior.tox[i] + sum(exp(like.tox[-i]-like.tox[i])*mprior.tox[-i]))
          }
          if (nord.tox > 1) {
            mtox.sel <- which.is.max(postprob.tox)  # selected ordering 
          }
          else {
            mtox.sel <- 1
          }
          mod[mtox.sel] <- mod[mtox.sel] + 1
          ptox.hat <- alpha[mtox.sel, ]^est.tox[mtox.sel]  # estimated toxicity probability
          loss <- abs(ptox.hat - theta)
          comb.curr <- which.is.max(-loss)                 # selected current dose
          if (sum(npts) == n) {
            stoprule <- 0
            break
          }
          else {
            y[comb.curr] <- y[comb.curr] + rbinom(1, 1, r[comb.curr])
            npts[comb.curr] <- npts[comb.curr] + 1
            i <- i+1
          }
        }
      }
      if (stoprule == 0) {
        comb.select[comb.curr] = comb.select[comb.curr] + 1    # output a vector with 1 at the selected dose
      }
      
      return(list(MTD.selection = comb.select, 
                  tox.data = y, 
                  patient.allocation = npts))
    }
  }
  lpocrm.sim <- function(nsim) {
    ncomb <- length(r)
    comb.select <- y <- npts <- matrix(nrow = nsim, ncol = ncomb)
    nstop <- 0
    for (i in 1:nsim) {
      if(i%%100==0) cat(paste("i = ", i, "\n"))
      result <- lpocrm(r, alpha, prior.o, x0, stop, n, theta)
      comb.select[i, ] = result$MTD.selection
      y[i, ] = result$tox.data
      npts[i, ] = result$patient.allocation
    }
    return(list(true.prob = r, 
                MTD.selection = colMeans(comb.select), 
                patient.allocation = round(colMeans(npts), 3), 
                tox.data = round(colMeans(y), 3)
    ))
  }
  lpocrm.sim(nsim)
}