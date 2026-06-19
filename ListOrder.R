############################################################
############################################################
##                                                        ##
##  List all possible orderings on two-dimensional grids  ##
##  Weishi Chen                                           ##
##  Last update: 19 Jun 2026                              ##
##                                                        ##
############################################################
############################################################
library(dplyr)


ListOrder <- function(dim, nsims=5000) {
  # Inputs:
  # -------------------------------------------
  # dim  : 2-vector, the dimension of the grid
  # nsims: integer, the number of simulated toxicity scenarios
  
  if(length(dim)!=2) stop("This function only works for two-dimensional grids")
  if((dim[1]<2) | (dim[2]<2)) stop("Both agents need to have at least 2 levels")
  
  # Simulate a random toxicity scenario R
  # and extract the ordering of R
  sim.scen <- function(dim) {
    R <- matrix(nrow=dim[1], ncol=dim[2])
    R[1,1] <- 0
    R[2:dim[1], 1] <- sort(runif(dim[1]-1))
    R[1, 2:dim[2]] <- sort(runif(dim[2]-1))
    for (i in 2:dim[1]) {
      for (j in 2:dim[2]) {
        R[i,j] <- runif(1, min=max(R[i, j-1], R[i-1, j]))
      }
    }
    O <- order(as.vector(R))
    return(O)
  }
  
  # Simulate *nsims* orderings and take the unique ones
  order.sim <- t(replicate(nsims, sim.scen(dim=dim)))
  orderings <- unique(order.sim)
  
  # arrange the orderings
  colnames(orderings) <- paste0("d", 1:prod(dim))
  orderings <- orderings %>% as_tibble() %>% 
    arrange(pick(starts_with("d"))) %>% 
    as.matrix()
  
  return(orderings = orderings)
}


# to obtain all 42 orderings on the 3x3 grid
# Each row gives a possible ordering
# Within each row, the numbers 1, 2, ..., k corresponds to d_{1,1}, d_{2,1}, ..., d_{c,1}, d_{1,2}, ..., d_{c,2}, ..., d_{c,r}.
ListOrder(dim=c(3, 3))
