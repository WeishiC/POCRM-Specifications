################################################
################################################
##                                            ##
##  Scenario-agnostic ordering specification  ##
##  Weishi Chen                               ##
##  Last update: 19 Jun 2026                  ##
##                                            ##
################################################
################################################

library(dplyr)

get.order.scen <- function(nsims=1e+4, dim, TTL) {
  # MTC_list: matrix, each row gives the location of the MTC
  # labels: vector of length=nrow(MTC_list), labels of the MTCs
  
  sim.scen <- function(MTC, dim, TTL) {
    dim1 <- dim[1]; dim2 <- dim[2]
    y <- matrix(nrow=dim1, ncol=dim2)
    y[MTC[1], MTC[2]] <- TTL
    # MTC row
    if(MTC[1]>1) y[1:(MTC[1]-1),MTC[2]] <- sort(runif(MTC[1]-1, max=TTL))
    if(MTC[1]<dim1) y[(MTC[1]+1):dim1,MTC[2]] <- sort(runif(dim1-MTC[1], min=TTL))
    # MTC col
    if(MTC[2]>1) y[MTC[1],1:(MTC[2]-1)] <- sort(runif(MTC[2]-1, max=TTL))
    if(MTC[2]<dim2) y[MTC[1],(MTC[2]+1):dim2] <- sort(runif(dim2-MTC[2], min=TTL))
    
    # bottom-left
    if((MTC[1]>1) & (MTC[2]>1)) {
      for (i in rev(1:(MTC[1]-1))) {
        for (j in rev(1:(MTC[2]-1))) {
          y[i,j] <- runif(1, max=min(y[i+1,j], y[i,j+1]))
        }
      }
    }
    # bottom-right
    if((MTC[1]<dim1) & (MTC[2]>1)) {
      for (i in (MTC[1]+1):dim1) {
        for (j in rev(1:(MTC[2]-1))) {
          y[i,j] <- runif(1, min=y[i-1,j], max=y[i,j+1])
        }
      }
    }
    # top-left
    if((MTC[1]>1) & (MTC[2]<dim2)) {
      for (i in rev(1:(MTC[1]-1))) {
        for (j in (MTC[2]+1):dim2) {
          y[i,j] <- runif(1, min=y[i,j-1], max=y[i+1,j])
        }
      }
    }
    # top-right
    if((MTC[1]<dim1) & (MTC[2]<dim2)) {
      for (i in (MTC[1]+1):dim1) {
        for (j in (MTC[2]+1):dim2) {
          y[i,j] <- runif(1, min=max(y[i-1,j], y[i,j-1]))
        }
      }
    }
    return(as.vector(y))
  }
  
  MTC_list <- cbind(rep(1:dim[1],dim[2]), rep(1:dim[2], each=dim[1]))
  L <- prod(dim)
  labels <- 1:L
  
  all.scen <- list()
  order_selected <- list()
  order_selected[[1]] <- 1:L
  if(is.vector(MTC_list)) MTC_list <- matrix(MTC_list, nrow=1)
  u <- 1
  for(c in 1:nrow(MTC_list)){
    print(c)
    # simulate nsims scenarios with MTC at specified place  
    all.scen[[c]] <- list()
    scen.store <- t(sapply(1:nsims, function(i) {
      sim.scen(MTC=MTC_list[c,], dim=dim, TTL=TTL)
    }))
    
    colnames(scen.store) <- paste0("d", 1:L)
    # labels of the MTCs
    label.store <- apply(scen.store, 1, function(u) which(sort(u)==TTL))
    label.unique <- unique(label.store)
    scen_curr <- scen.store %>% as_tibble() %>% mutate(label=label.store)
    scen.index <- 1
    for (i in label.unique) {
      # simulated scenarios with the specified label
      scen_temp <- scen_curr %>% mutate(n.scen=1:nsims) %>% filter(label==i) %>% as.matrix()
      # needed scenario, i.e. one for each below MTC set
      scen.needed <- t(apply(matrix(scen_temp[,1:L], ncol=L), 1, function(u) order(u)))[,1:i] %>%
        unique() 
      if(is.vector(scen.needed)) scen.needed <- matrix(scen.needed, nrow=1)
      # the order among the below MTC set doesn't matter
      if(i>1) {
        scen.needed <- cbind(t(apply(matrix(scen.needed[,1:(i-1)], ncol=(i-1)), 1, sort)) %>% unique(), labels[c])
        colnames(scen.needed) <- NULL
      }
      # add one order-scenario for each below MTC set
      scen_include <- numeric()
      for (j in 1:nrow(scen.needed)) {
        scen.orders.needed <- apply(scen_temp, 1, function(u) all.equal(sort(order(as.vector(unlist(u)))[1:i]),sort(scen.needed[j,])))
        scen.include_temp <- scen_temp %>% as_tibble() %>% mutate(include=scen.orders.needed) %>%
          filter(include==TRUE) %>% select(n.scen) %>% unlist() %>% as.vector()
        scen_include[j] <- scen.include_temp[1]
      }
      all.scen[[c]][[scen.index]] <- scen_curr[scen_include,1:L]
      scen.index <- scen.index + 1
      # under each below MTC set
      for (j in 1:nrow(scen.needed)) {
        # test if the selected orderings has already included one from correct group
        consis <- mapply(function(x, y){
          isTRUE(all.equal(x[1:i], y[1:i]))
        }, order_selected, list(scen.needed))
        # if not, add one from the correct group
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
  # collect results
  # orderings added
  order.selected.mat <- order_selected %>% unlist() %>% matrix(byrow=TRUE, ncol=L)
  # order scenarios
  all.scen.temp <- all.scen %>% unlist(recursive = FALSE)
  all.scen.mat <- matrix(NA, nrow=1, ncol=L)
  colnames(all.scen.mat) <- paste0("d", 1:L)
  for (i in 1:length(all.scen.temp)) {
    all.scen.mat <- rbind(all.scen.mat, all.scen.temp[[i]])
  }
  all.scen.mat <- all.scen.mat[-1,]
  all.scen.mat <- all.scen.mat %>% as.matrix()
  colnames(all.scen.mat) <- NULL
  # number of order-scenario by MTC location
  n.scen <- sapply(all.scen, function(list){
    u <- 0
    for (i in 1:length(list)) {
      u <- u + nrow(list[[i]])
    }
    u
  })
  # total number of order-scenario
  n <- sum(n.scen)
  return(list(OrderScen=all.scen.mat, order.adding=order.selected.mat, n.scen=n.scen, n=n))
}


# on 3x3 grid

# obtain all order-scenarios
sim <- get.order.scen(nsims=1000, dim=c(3,3), TTL=0.3)
order.scen <- sim$OrderScen     

# obtain all orderings
source("ListOrder.R")
orderings <- ListOrder(dim=c(3,3)) 

# specify orderings
source("GetOrder.R")
scen.agn <- get.order(all.scen=order.scen, all.orderings=orderings, TTL=0.30, S=6)
