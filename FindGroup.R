FindGroups <- function(all.scen, all.orderings, TTL, ids=NULL) {
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
  if(is.null(ids)) ids <- 1:nrow(all.orderings)
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
    Group[[c]] <- ids[sort(unique(unlist(correct.group[[c]])))]
  }
  return(Group)
}

# groups <- FindGroups(all.scen=all.scen, all.orderings=all.orderings, TTL=TTL)$Correct.group