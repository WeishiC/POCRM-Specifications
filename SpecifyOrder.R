#####################################
#####################################
##                                 ##
##  The Adding-Refining Algorithm  ##
##  Weishi Chen                    ##
##  Last update: 30th Mar 2025     ##
##                                 ##
#####################################
#####################################
library(pocrm); library(dfcrm)  # for CRM-related functions
library(nnet)
library(combinat)               # for listing orderings
library(dplyr); library(tidyr)
library(ggplot2)

########################################################
#
#  List orderings using the simulation method
#  Update: 14th Feb 2025
dim1 <- 3
dim2 <- 4
dim3 <- 2
nsims<-10000

orders.store <- t(sapply(1:nsims, function(sim){
  y <- array(dim=c(dim1, dim2, dim3))
  y[1,1,1] <- 0
  y[1,1,2:dim3] <- sort(runif(dim3-1, min=y[1,1,1]))
  y[1,2:dim2,1] <- sort(runif(dim2-1, min=y[1,1,1]))
  y[2:dim1,1,1] <- sort(runif(dim1-1, min=y[1,1,1]))
  for (i in 2:dim1) {
    for (j in 2:dim2) {
      y[i,j,1] <- runif(1, min=max(y[i, j-1, 1], y[i-1, j, 1]))
    }
  }
  for (index3 in 2:dim3) {
    for (j in 2:dim2) {
      y[1,j,index3] <- runif(1, min=max(y[1,j-1,index3], y[1,j,index3-1]))
    }
    for (i in 2:dim1) {
      y[i,1,index3] <- runif(1, min=max(y[i-1,1,index3], y[i,1,index3-1]))
    }
    for (i in 2:dim1) {
      for (j in 2:dim2) {
        y[i,j,index3] <- runif(1, min=max(y[i, j-1, index3], y[i-1, j, index3], y[i, j, index3-1]))
      }
    }
  }
  y <- as.vector(y)
  y <- y[c(1, 4, 5, 8, 9, 11, 12, 17, 20, 21, 23, 24)]
  order(as.vector(t(y)))
}))
all.orderings <- unique(orders.store)
M <- nrow(all.orderings)
L <- ncol(all.orderings)
colnames(all.orderings) <- paste0("d", 1:L)
all.orderings <- all.orderings %>% as_tibble() %>%
  arrange(d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12)
all.orderings <- all.orderings %>% as.matrix()

#
#  Simulate toxicity scenarios given the MTC
#
sim.scen <- function(sim, MTC, dim, TTL=0.25){
  # Inputs:
  #--------------------------------------------------------
  # sim: integer, index of simulation
  # MTC: vector, the MTC.
  # dim: vector, the dimension of the grid.
  # TTL: target toxicity level, value in [0,1].
  # -------------------------------------------------------
  dim1 <- dim[1]; dim2 <- dim[2]; dim3 <- dim[3]
  y <- array(dim=c(dim[1], dim[2], dim[3]))
  fix <- MTC
  y[fix[1], fix[2], fix[3]] <- TTL
  if(fix[1]>1) y[1:(fix[1]-1),fix[2],fix[3]] <-sort(runif(fix[1]-1, max=TTL))
  if(fix[1]<dim1) y[(fix[1]+1):dim1,fix[2],fix[3]] <-sort(runif(dim1-fix[1], min=TTL))
  if(fix[2]>1) y[fix[1],1:(fix[2]-1),fix[3]] <-sort(runif(fix[2]-1, max=TTL))
  if(fix[2]<dim2) y[fix[1],(fix[2]+1):dim2,fix[3]] <-sort(runif(dim2-fix[2], min=TTL))
  if(fix[3]>1) y[fix[1],1:fix[2],1:(fix[3]-1)] <-sort(runif(fix[3]-1, max=TTL))
  if(fix[3]<dim3) y[fix[1],fix[2],(fix[3]+1):dim3] <-sort(runif(dim3-fix[3], min=TTL))
  if(fix[3]==1) {
    if(fix[1]>1) {
      for (i in rev(1:(fix[1]-1))) {
        if(fix[2]>1){
          for (j in rev(1:(fix[2]-1))) {
            y[i,j,1] <- runif(1, max=min(y[i, j+1, 1], y[i+1, j, 1]))
          }
        }
        if(fix[2]<dim2) {
          for (j in (fix[2]+1):dim2) {
            y[i,j,1] <- runif(1, min=y[i,j-1,1], max=y[i+1,j,1])
          }
        }
      }
    }
    if(fix[1]<dim1) {
      for (i in (fix[1]+1):dim1) {
        if(fix[2]>1){
          for (j in rev(1:(fix[2]-1))) {
            y[i,j,1] <- runif(1, min=y[i-1,j,1], max=y[i, j+1, 1])
          }
        }
        if(fix[2]<dim2) {
          for (j in (fix[2]+1):dim2) {
            y[i,j,1] <- runif(1, min=max(y[i, j-1, 1], y[i-1, j, 1]))
          }
        }
      }
    }
    for (index3 in (fix[3]+1):dim3) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          y[i,fix[2],index3] <- runif(1, min=y[i,fix[2],index3-1], max=y[i+1,fix[2],index3])
        }
      }
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          y[i,fix[2],index3] <- runif(1, min=max(y[i,fix[2],index3-1], y[i-1,fix[2],index3]))
        }
      }
      if(fix[2]>1) {
        for (j in rev(1:(fix[2]-1))) {
          y[fix[1],j,index3] <- runif(1, min=y[fix[1],j,index3-1], max=y[fix[1],j+1,index3])
        }
      }
      if(fix[2]<dim2) {
        for (j in (fix[2]+1):dim2) {
          y[fix[1],j,index3] <- runif(1, min=max(y[fix[1],j,index3-1], y[fix[1],j-1,index3]))
        }
      }
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, min=y[i,j,index3-1], max=min(y[i, j+1, index3], y[i+1, j, index3]))
            }
          }
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=max(y[i,j-1,index3], y[i,j,index3-1]), max=y[i+1,j,index3])
            }
          }
        }
      }
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, min=max(y[i-1,j,index3], y[i,j,index3-1]), max=y[i, j+1, index3])
            }
          }
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=max(y[i, j-1, index3], y[i-1, j, index3], y[i,j,index3-1]))
            }
          }
        }
      }
    }
  }
  if(fix[3]>1) {
    if(fix[1]>1) {
      for (i in rev(1:(fix[1]-1))) {
        if(fix[2]>1){
          for (j in rev(1:(fix[2]-1))) {
            y[i,j,fix[3]] <- runif(1, max=min(y[i, j+1, fix[3]], y[i+1, j, fix[3]]))
          }
        }
        if(fix[2]<dim2) {
          for (j in (fix[2]+1):dim2) {
            y[i,j,fix[3]] <- runif(1, min=y[i,j-1,fix[3]], max=y[i+1,j,fix[3]])
          }
        }
      }
    }
    if(fix[1]<dim1) {
      for (i in (fix[1]+1):dim1) {
        if(fix[2]>1){
          for (j in rev(1:(fix[2]-1))) {
            y[i,j,fix[3]] <- runif(1, min=y[i-1,j,fix[3]], max=y[i, j+1, fix[3]])
          }
        }
        if(fix[2]<dim2) {
          for (j in (fix[2]+1):dim2) {
            y[i,j,fix[3]] <- runif(1, min=max(y[i, j-1, fix[3]], y[i-1, j, fix[3]]))
          }
        }
      }
    }
    for (index3 in rev(1:(fix[3]-1))) {
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          y[i,fix[2],index3] <- runif(1, max=min(y[i+1,fix[2],index3], y[i,fix[2],index3+1]))
        }
      }
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          y[i,fix[2],index3] <- runif(1, min=y[i-1,fix[2],index3], max=y[i,fix[2],index3+1])
        }
      }
      if(fix[2]>1) {
        for (j in rev(1:(fix[2]-1))) {
          y[fix[1],j,index3] <- runif(1, max=min(y[fix[1],j+1,index3], y[fix[1],j,index3+1]))
        }
      }
      if(fix[2]<dim2) {
        for (j in (fix[2]+1):dim2) {
          y[fix[1],j,index3] <- runif(1, min=y[fix[1],j-1,index3], max=y[fix[1],j,index3+1])
        }
      }
      if(fix[1]>1) {
        for (i in rev(1:(fix[1]-1))) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, max=min(y[i, j+1, index3], y[i+1, j, index3], y[i,j,index3+1]))
            }
          }
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=y[i,j-1,index3], max=min(y[i+1,j,index3], y[i,j,index3+1]))
            }
          }
        }
      }
      if(fix[1]<dim1) {
        for (i in (fix[1]+1):dim1) {
          if(fix[2]>1){
            for (j in rev(1:(fix[2]-1))) {
              y[i,j,index3] <- runif(1, min=y[i-1,j,index3], max=min(y[i, j+1, index3], y[i,j,index3+1]))
            }
          }
          if(fix[2]<dim2) {
            for (j in (fix[2]+1):dim2) {
              y[i,j,index3] <- runif(1, min=max(y[i, j-1, index3], y[i-1, j, index3]), max=y[i,j,index3+1])
            }
          }
        }
      }
    }
  }
  
  
  y <- as.vector(y)
  y <- y[c(1, 4, 5, 8, 9, 11, 12, 17, 20, 21, 23, 24)]
  y
}


# Setting
MTC_list <- matrix(c(1, 1, 1,
                     1, 2, 1, 
                     2, 2, 1, 
                     2, 3, 1, 
                     3, 3, 1,
                     2, 4, 1,
                     3, 4, 1,
                     2, 2, 2, 
                     2, 3, 2,
                     3, 3, 2, 
                     2, 4, 2, 
                     3, 4, 2), byrow=TRUE, ncol=3)
nsims <- 10000
dim <- c(3, 4, 2)
TTL <- 0.25
L <- 12

#
# Adding step
#
Adding <- function(nsims, MTC_list, dim, TTL, L=12) {
  all.scen <- list()
  order_selected <- list()
  order_selected[[1]] <- 1:L
  u <- 1
  for(c in 1:L){
    print(c)
    all.scen[[c]] <- list()
    scen.store <- t(sapply(1:nsims, sim.scen, MTC=MTC_list[c,], dim, TTL=0.25))
    colnames(scen.store) <- paste0("d", 1:L)
    label.store <- apply(scen.store, 1, function(u) which(sort(u)==TTL))
    label.unique <- unique(label.store)
    scen_curr <- scen.store %>% as_tibble() %>% mutate(label=label.store)
    scen.index <- 1
    for (i in label.unique) {
      scen_temp <- scen_curr %>% mutate(n.scen=1:nsims) %>% filter(label==i) %>% as.matrix()
      scen.needed <- t(apply(matrix(scen_temp[,1:L], ncol=L), 1, function(u) order(u)))[,1:i] %>%
        unique() %>% as.matrix()
      
      
      if(ncol(scen.needed)==1) {
        scen_include <- numeric()
        scen.orders.needed <- apply(scen_temp, 1, function(u) all.equal(sort(order(as.vector(unlist(u)))[1:i]),sort(scen.needed)))
        scen.include_temp <- scen_temp %>% as_tibble() %>% mutate(include=scen.orders.needed) %>%
          filter(include==TRUE) %>% select(n.scen) %>% unlist() %>% as.vector()
        scen_include <- scen.include_temp[1]
        all.scen[[c]][[scen.index]] <- scen_curr[scen_include,1:L]
        scen.index <- scen.index + 1
        
        f1 <- function(x,y) all(x[1:i]==y[1:i])
        consis <- mapply(f1, order_selected, list(scen.needed))
        if(sum(consis)==0) {
          u <- u+1
          orders.needed <- apply(scen_temp, 1, function(u) all.equal(sort(order(as.vector(unlist(u)))[1:i]),sort(scen.needed)))
          which(orders.needed==TRUE)
          include_temp <- scen_temp %>% as_tibble() %>% mutate(include=orders.needed) %>%
            filter(include==TRUE) %>% select(-label, -n.scen, -include) %>% as.matrix()
          order(include_temp[1,1:L])
          order_selected[[u]] <- order(include_temp[1,1:L])
        }
      } else {
        if(c>1) {
          scen.needed <- cbind(t(apply(matrix(scen.needed[,1:(i-1)], ncol=(i-1)), 1, sort)) %>% unique(), c)
          colnames(scen.needed) <- NULL
        }
        scen_include <- numeric()
        for (j in 1:nrow(scen.needed)) {
          scen.orders.needed <- apply(scen_temp, 1, function(u) all.equal(sort(order(as.vector(unlist(u)))[1:i]),sort(scen.needed[j,])))
          scen.include_temp <- scen_temp %>% as_tibble() %>% mutate(include=scen.orders.needed) %>%
            filter(include==TRUE) %>% select(n.scen) %>% unlist() %>% as.vector()
          scen_include[j] <- scen.include_temp[1]
        }
        all.scen[[c]][[scen.index]] <- scen_curr[scen_include,1:L]
        scen.index <- scen.index + 1
        
        for (j in 1:nrow(scen.needed)) {
          f1 <- function(x,y) all(x[1:i]==y[1:i])
          consis <- mapply(f1, order_selected, list(scen.needed[j,]))
          if(sum(consis)==0) {
            u <- u+1
            orders.needed <- apply(scen_temp, 1, function(u) all.equal(sort(order(as.vector(unlist(u)))[1:i]),sort(scen.needed[j,])))
            include_temp <- scen_temp %>% as_tibble() %>% mutate(include=orders.needed) %>%
              filter(include==TRUE) %>% select(-label, -n.scen, -include) %>% as.matrix()
            order_selected[[u]] <- order(include_temp[1,1:L])
          }
        }
      }
    }
  }
  
  # collect results
  order.selected.mat <- order_selected %>% unlist() %>% matrix(byrow=TRUE, ncol=L)
  all.scen.temp <- all.scen %>% unlist(recursive = FALSE)
  all.scen.mat <- matrix(NA, nrow=1, ncol=L)
  colnames(all.scen.mat) <- paste0("d", 1:L)
  for (i in 1:length(all.scen.temp)) {
    all.scen.mat <- rbind(all.scen.mat, all.scen.temp[[i]])
  }
  all.scen.mat <- all.scen.mat[-1,]
  all.scen.mat <- all.scen.mat %>% as.matrix()
  colnames(all.scen.mat) <- NULL
  return(list(OrderScen=all.scen.mat, order.adding=order.selected.mat))
}

# the order-scenarios stored in *OrderScen*
# the orderings included by the adding step stored in *OrderAdding*
adding.out <- Adding(nsims, MTC_list, dim, TTL, L=12)
OrderScen <- adding.out$OrderScen
OrderAdding <- adding.out$order.adding



#
#  Refining step
#
Refining <- function(max.sim, all.scen, all.orderings, S, order.fix=NULL, TTL=0.25) {
  relabel <- function(ordering, label) {
    new.ordering <- numeric()
    for (i in 1:length(ordering)) {
      new.ordering[which(ordering==i)] <- label[i]
    }
    new.ordering
  }
  group.models <- function(u, MTC) {
    # u: k-vector, one complete ordering
    nu <- which(u==MTC)
    if(nu>1) {
      w <- sum(u[1:(nu-1)]>MTC) # number of more toxic dose ordered before MTC
    } else {
      w <- 0
    }
    if(nu<k) {
      n <- sum(u[(nu+1):k]<MTC) # number of less toxic dose ordered after MTC
    } else {
      n <- 0
    }
    
    c(w,n)
  }
  k <- ncol(all.orderings)
  M <- nrow(all.orderings)
  G <- nrow(all.scen)
  out <- numeric(S)
  # select orderings
  Continue <- TRUE
  sim <- 1
  while (Continue) {
    if(is.null(order.fix)) {
      sel <- sample(1:M, S, replace = FALSE)
    } else {
      S1 <- length(order.fix)
      sel <- c(order.fix, sample((1:M)[-order.fix], S-S1, replace=FALSE))
    }
    ordering.sel <- all.orderings[sel,]
    y <- rep(TRUE, G)
    for (c in 1:G) {
      new.label <- order(order(all.scen[c,]))
      ordering.new <- apply(ordering.sel, 1, relabel, label=new.label) %>% t()
      MTC_j <- apply(ordering.new, 1, group.models, MTC=which(sort(all.scen[c,])==TTL))
      rownames(MTC_j) <- c("w", "n")
      MTC_j <- MTC_j %>% t() %>% as_tibble() %>% mutate(model=sel, .before=w)
      consis_c <- MTC_j %>% filter(w==0, n==0) %>% nrow()
      if(consis_c==0) y[c] <- FALSE
    }
    if(all(y)==TRUE) {
      out <- sel
      Continue <- FALSE
    } 
    if(sim>=max.sim) {
      Continue <- FALSE
    } 
    sim <- sim + 1
    if((sim %% 50)==0) print(sim)
  }
  out
}

OrderRefine <- Refining(max.sim=10000, OrderScen, OrderAdding, S=9, order.fix=NULL, TTL=0.25)
# If you get a vector of 0's, that means it's not possible to find consistent orderings
# with only *S* orderings, try to increase *max.sim* to 10^6, if still not possible, 
# increase *S*.


# change back to the index among all 148 orderings
order.refined <- OrderAdding[OrderRefine,]
ordering_filter <- function(ordering1, ordering2) {
  M1 <- nrow(ordering1); M2 <- nrow(ordering2)
  overlap <- numeric(M1)
  for (i in 1:M1) {
    O1 <- ordering1[i,]
    colnames(ordering2) <- paste0("d", 1:ncol(ordering2))
    overlap[i] <- cbind(order=1:M2, ordering2) %>% as_tibble() %>% filter(
      d1==O1[1], d2==O1[2], d3==O1[3], d4==O1[4],
      d5==O1[5], d6==O1[6], d7==O1[7], d8==O1[8],
      d9==O1[9], d10==O1[10], d11==O1[11], d12==O1[12]
    ) %>% select(order) %>% as.integer()
  }
  overlap
}
ordering_filter(order.refined%>%as.matrix(), all.orderings%>%as.matrix())
