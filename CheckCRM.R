library(pocrm)

CheckCRM <- function(R, alpha, ordering, TTL=0.3) {
  J <- length(R)   # number of dose levels
  out <- list()
  out$R <- R
  MTD <- which.min(abs(sort(R)-TTL))  # the true MTD
  
  # calculate a_j
  a_j <- log(R[ordering])/log(alpha)
  a_j_grid <- getwm(matrix(ordering, nrow=1), a_j)
  out$a_j <- as.vector(a_j)
  
  b <- numeric(J+1)
  b[1] <- 0
  b[J+1] <- 100
  
  # calculate the interval H_j's
  CRM <- function(d, a) d^a
  for (j in 2:J) {
    f.b <- function(b) (CRM(alpha[j-1], b)+CRM(alpha[j], b)-2*TTL)^2
    b[j] <- optimise(f.b, interval = c(0, 10))$minimum
  }
  out$b <- b
  consis1 <- (a_j[MTD]>b[MTD]) & (a_j[MTD]<b[MTD+1])
  consis_detail <- consis1
  consis2 <- ifelse(MTD>1, all(sapply(1:(MTD-1), function(j) a_j[j]>b[j+1])), TRUE)
  if(MTD>1) consis_detail <- c(sapply(1:(MTD-1), function(j) a_j[j]>b[j+1]), consis_detail)
  consis3 <- ifelse(MTD<J, all(sapply((MTD+1):J, function(j) a_j[j]<b[j])), TRUE)
  if(MTD<J) consis_detail <- c(consis_detail, sapply((MTD+1):J, function(j) a_j[j]<b[j]))
  consis <- all(consis1, consis2, consis3)
  out$consis <- consis
  out$detail <- consis_detail
  return(out)
}
