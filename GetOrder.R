################################################
################################################
##                                            ##
##  Scenario-specific ordering specification  ##
##  Weishi Chen                               ##
##  Last update: 19 Jun 2026                  ##
##                                            ##
################################################
################################################
library(dplyr)


# Specify orderings consistent under a given set of scenarios
get.order <- function(max.sim=1e+4, all.scen, all.orderings, TTL=0.30, S, order.fix=NULL) {
  # Inputs:
  # ------------------------------
  # max.sim      : integer, maximum number of simulations.
  # all.scen     : matrix, each row contains a toxicity scenario.
  # all.orderings: matrix, each row contains a possible ordering.
  # TTL          : numeric in [0,1], the target toxicity level.
  # S            : integer, number of orderings to include.
  # order.fix    : vector of length <=S, orderings must include.
  
  # Given a set of scenarios, find the orderings in the correct group under each scenario
  FindGroups <- function(all.scen, all.orderings, TTL) {
    # Allow multiple MTCs
    # Allow the MTC with toxicity not exactly the TTL
    # Inputs:
    # --------------------------------------
    # all.scen     : Rxk matrix, each row contains a scenario. 
    # all.orderings: Mxk matrix, each row contains an ordering.
    # TTL          : [0,1], target toxicity level.
    # Outputs:
    # -------------------------------------
    # Correct.group    : R-list, each list contains a vector of the correct group under all MTCs.
    # Correct.group.MTC: R-list, each list contains a list of length equal to the number of MTCs. 
    #                    the correct gourp under each MTC.
    relabel <- function(ordering, label) {
      # Relabel combinations s.t. toxicity increases 1 -> ... -> k.
      new.ordering <- numeric()
      for (i in 1:length(ordering)) {
        new.ordering[which(ordering==i)] <- label[i]
      }
      new.ordering
    }
    group.models <- function(u, MTC.label, scen, new.label, MTC.tox) {
      # Inputs:
      # -------------------------------------
      # u        : k-vector, one ordering (under new label).
      # MTC.label: integer, the new label of the MTC.
      # scen     : toxicity scenario
      # new.label: k-vector, new labels of combinations.
      # MTC.tox  : [0,1], toxicity probability of the MTC, not necessarily the TTL.
      nu <- which(u==MTC.label)
      label <- tibble(label = new.label, tox = scen) %>% as.matrix()
      # toxicity under the specified ordering u
      # note here we are using *scen* not *scen_temp*, i.e. did not add 0.001 to the MTCs.
      scen_u <- numeric()  
      for (i in 1:length(scen)) {
        scen_u[i] <- label[new.label==u[i], "tox"]
      }
      if(nu==length(u) | nu==1) {
        w <- n <- 0
      } else {
        w <- sum(scen_u[1:(nu-1)]>MTC.tox) # number of more toxic dose ordered before MTC
        n <- sum(scen_u[(nu+1):k]<MTC.tox) # number of less toxic dose ordered after MTC
      }
      c(w,n)
    }
    
    k <- ncol(all.orderings)     # number of combo's
    M <- nrow(all.orderings)     # number of orderings
    if(is.vector(all.scen)) all.scen <- matrix(all.scen, nrow=1)
    correct.group <- Group <- list()
    for (c in 1:nrow(all.scen)) {
      scen <- all.scen[c,]
      MTC_c <- which(abs(scen-TTL)==min(abs(scen-TTL)))    # all MTCs
      # look for the correct group under each MTC, 
      # while adding 0.001 to all other MTCs
      correct.group[[c]] <- list()
      for (j in 1:length(MTC_c)) {
        scen_temp <- scen
        scen_temp[MTC_c] <- scen_temp[MTC_c]+0.001
        MTC <- MTC_c[j]
        scen_temp[MTC] <- scen[MTC]
        # the new labels of the combination
        new.label <- order(order(scen_temp))
        # orderings under the new label
        ordering.new <- apply(all.orderings, 1, relabel, label=new.label) %>% t()
        # compute the groups
        MTC.label <- which.min(abs(sort(scen_temp)-scen[MTC]))
        MTC_j <- apply(ordering.new, 1, group.models, MTC.label=MTC.label, scen=scen, new.label=new.label, MTC.tox=scen[MTC])
        rownames(MTC_j) <- c("w", "n")
        # correct group
        correct.group[[c]][[j]] <- MTC_j %>% t() %>% as_tibble()%>% mutate(model=1:M, .before=w) %>% 
          filter(w==0, n==0) %>% select(model) %>% unlist() %>% as.vector()
      }
      # the correct group under any of the MTC
      Group[[c]] <- sort(unique(unlist(correct.group[[c]])))
    }
    return(list(Correct.group=Group, Correct.group.MTC=correct.group))
  }
  
  groups <- FindGroups(all.scen=all.scen, all.orderings=all.orderings, TTL=TTL)$Correct.group
  
  M <- nrow(all.orderings)
  
  G <- length(groups)
  out <- numeric(S)
  # maximum number of unique selection
  if(is.null(order.fix)) {
    sample.max <- choose(M, S)
  } else {
    sample.max <- choose(M-length(order.fix), S-length(order.fix))
  }
  max.sim <- min(max.sim, sample.max)
  
  # select orderings
  Continue <- TRUE
  sim <- 1
  seen_signatures <- character(0)
  while (Continue) {
    sample.done <- FALSE
    while(!sample.done){
      if(is.null(order.fix)) {
        sel <- sample(1:M, S, replace = FALSE)
      } else {
        S1 <- length(order.fix)
        sel <- c(order.fix, sample((1:M)[-order.fix], S-S1, replace=FALSE))
      }
      # avoid testing over the same subset of orderings.
      signature <- paste(sort(sel), collapse = "_")
      if (!(signature %in% seen_signatures)) {
        seen_signatures <- c(seen_signatures, signature)
        sample.done <- TRUE
      }
    }
    
    y <- rep(TRUE, G)
    for (j in 1:G) {
      if(sum(sel %in% groups[[j]])==0) {
        y[j] <- FALSE
      }
    }
    if(all(y)==TRUE) {
      out <- all.orderings[sort(sel),]
      Continue <- FALSE
    } 
    if(sim>=max.sim) {
      out <- NULL
      cat(paste0("No ", S, " orderings can be found, output ", length(out), " orderings.\n"))
      Continue <- FALSE
    } 
    sim <- sim + 1
  }
  return(out)
}



# For 3x3 grid on 19 scenarios
load("scen3x3.RData")
source("ListOrder.R")

orderings <- ListOrder(dim=c(3,3))

scen.spec <- get.order(max.sim=500, all.scen=scen, all.orderings=orderings, TTL=0.3, S=3)
