## =====================================================================
## Check of Condition 2 of Theorem 2 ONLY (no MTC-limit check).
##
## Choice of restriction on the directions zeta (argument
## `restriction`):
##   "tox"      (default) a combination in W2 receives no patients whenever
##              the combination in T2 ordered immediately before it (under
##              the incorrect ordering) receives no patients; similarly, a
##              combination in T1 receives no patients whenever its adjacent
##              combination in W1 (the W1 combination ordered before it,
##              possibly via a run of consecutive T1 combinations) receives
##              no patients; and, within each run of two or more consecutive
##              T1 combinations, the first combination of the run receives at
##              least as many patients as each of the others;
##   "adjacent"  an over-toxic combination receives no more patients than the
##              nearest earlier combination below the TTL; and, within each run
##              of two or more consecutive T1 combinations, the first combination
##              of the run receives at least as many patients as each of the others;
##   "none"      no restriction.
##
## For an incorrect ordering O_m, the condition requires an included
## correct ordering O_t such that, for every allocation eta supported on
## S(m;R) = T(m;R) ∪ W(m;R) ∪ {MTC},
##
##   LHS(eta) = sum_{(i,j) in S°} eta_ij [ f(Rhat^(t)_ij(eta), R_ij)
##                                       - f(Rhat^(m)_ij(eta), R_ij) ] > 0,
##
## where S° = S \ {MTC}. Allocations are generated along rays
##   eta = (1 - u) * zeta + u * e_MTC,
## with zeta a direction on S° and u the share of patients at the MTC.
##
## Excluded by construction:
##  * eta = e_MTC (u = 1): all patients at the MTC;
##  * all patients at a single combination of S° (u = 0, zeta a vertex):
##    both orderings fit it exactly, so LHS = 0 by construction.
##
## Conventions
##  * Combinations are indexed 1..k; a scenario R is a length-k vector.
##  * An ordering is a length-k vector of combination indices, least to
##    most toxic; `orderings` is an M x k matrix.
##  * The skeleton `skel` (increasing, length k) assigns skel[s] to the
##    combination in position s of an ordering.
##  * Single-MTC scenarios only (the MTC is the combination equal to TTL).
## =====================================================================

f_part <- function(p, R) R * log(p) + (1 - R) * log(1 - p)

## ---------------------------------------------------------------------
## Sets in Definitions 2 and 3
## ---------------------------------------------------------------------
## get.S(): the set S(m;R) for ONE ordering (Section 7.2 of the paper).
##  R        : true toxicity probabilities (length k, combination-index order)
##  ordering : one ordering (length-k vector of combination indices, least to
##             most toxic), typically an incorrect ordering
##  TTL      : target toxicity level
##  cnames   : optional labels of the combinations (e.g. "(4,1)"); if given, the
##             index vectors below are named by these labels
## Returns a list with
##  S        : S(m;R) = T(m;R) ∪ W(m;R) ∪ {MTC}, as combination indices
##  S0       : S(m;R) without the MTC
##  MTC      : index of the MTC
##  correct  : TRUE if the ordering belongs to the correct group (Definition 2),
##             in which case T, W and S0 are empty and S is the MTC only
##  T1, T2, W1, W2, T, W : the sets of Definition 3 (T = T1 ∪ T2, W = W1 ∪ W2)
get.S <- function(R, ordering, TTL = 0.3, cnames = NULL) {
  nu <- which(abs(R - TTL) < 1e-10)
  if (length(nu) != 1) stop("Scenario must have exactly one MTC (toxicity == TTL).")
  below <- which(R < R[nu]); nb <- length(below)
  p <- order(ordering)                            # p[c] = position of combination c
  correct <- ordering[nb + 1] == nu && setequal(ordering[seq_len(nb)], below)
  T1 <- which(R > R[nu] & p < p[nu])              # more toxic, ordered before the MTC
  T2 <- which(R < R[nu] & p > p[nu])              # less toxic, ordered after the MTC
  W1 <- setdiff(which((p + 1) %in% p[T1]), T1)    # immediately before a T1 combination
  W2 <- setdiff(which((p - 1) %in% p[T2]), T2)    # immediately after a T2 combination
  T <- sort(c(T1, T2)); W <- sort(c(W1, W2)); S0 <- sort(union(T, W))
  out <- list(S = sort(c(S0, nu)), S0 = S0, MTC = nu, correct = correct,
              T1 = T1, T2 = T2, W1 = W1, W2 = W2, T = T, W = W)
  if (!is.null(cnames)) {
    if (length(cnames) != length(R)) stop("cnames must have one label per combination.")
    for (nm in setdiff(names(out), "correct")) out[[nm]] <- setNames(out[[nm]], cnames[out[[nm]]])
  }
  out
}

## get.S_all(): get.S() for several orderings.
##  orderings : matrix with one ordering per row
##  ids       : labels of the orderings (e.g. c(1, 11, 36, ...)); default 1..M
##  cnames    : optional labels of the combinations
## Returns a data frame with one row per ordering: whether it belongs to the
## correct group, and the sets of Definition 3 and S0, written as labels (if
## cnames is given) or combination indices ("-" for an empty set). The full
## get.S() output of every ordering is attached as attr(, "S").
get.S_all <- function(R, orderings, ids = seq_len(nrow(orderings)), TTL = 0.3, cnames = NULL) {
  res <- lapply(seq_len(nrow(orderings)), function(m) get.S(R, orderings[m, ], TTL, cnames))
  lab <- function(v) if (!length(v)) "-" else if (is.null(cnames)) paste(v, collapse = ",") else paste(cnames[v], collapse = ",")
  out <- data.frame(ordering = ids,
                    correct  = sapply(res, `[[`, "correct"),
                    T1 = sapply(res, function(r) lab(r$T1)), W1 = sapply(res, function(r) lab(r$W1)),
                    T2 = sapply(res, function(r) lab(r$T2)), W2 = sapply(res, function(r) lab(r$W2)),
                    S0 = sapply(res, function(r) lab(r$S0)), stringsAsFactors = FALSE)
  attr(out, "MTC") <- if (is.null(cnames)) res[[1]]$MTC else cnames[res[[1]]$MTC]
  attr(out, "S") <- setNames(res, ids)
  out
}
## Internal: MTC, included correct orderings, positions, and the sets of every
## incorrect ordering (via get.S()).
get_sets <- function(R, orderings, TTL = 0.3) {
  pos <- t(apply(orderings, 1, order))            # pos[m, c] = position of combination c
  all <- lapply(seq_len(nrow(orderings)), function(m) get.S(R, orderings[m, ], TTL))
  correct <- which(sapply(all, `[[`, "correct"))
  sets <- list()
  for (m in setdiff(seq_len(nrow(orderings)), correct))
    sets[[as.character(m)]] <- all[[m]][c("T", "W", "T1", "W1", "T2", "W2")]
  list(nu = all[[1]]$MTC, correct = correct, sets = sets, pos = pos)
}

## ---------------------------------------------------------------------
## Grid of allocations along rays
## ---------------------------------------------------------------------
## all points of the simplex in d dimensions with coordinates in {0, 1/h, ..., 1}
simplex_lattice <- function(d, h) {
  if (d == 1) return(matrix(1, 1, 1))
  comp <- function(n, d) {
    if (d == 1) return(matrix(n, 1, 1))
    do.call(rbind, lapply(0:n, function(i) cbind(i, comp(n - i, d - 1))))
  }
  comp(h, d) / h
}
## default lattice resolution by number of combinations in S°
default_h <- function(d) if (d <= 2) 20 else if (d == 3) 10 else if (d == 4) 6 else if (d == 5) 4 else 3

## Pairs (toxic, lower) for the zeta restriction: for each combination of S0
## with toxicity above the TTL, `lower` is the nearest combination ordered
## before it under the incorrect ordering whose toxicity is below the TTL.
## Directions must satisfy zeta[toxic] <= zeta[lower].
toxic_pairs <- function(Sc, R, pos_m, TTL) {
  ord <- order(pos_m)                                  # combinations in the order of O_m
  out <- NULL
  for (c in Sc[R[Sc] > TTL + 1e-10]) {
    prev <- rev(ord[seq_len(pos_m[c] - 1)])            # combinations before c, nearest first
    low  <- prev[R[prev] < TTL - 1e-10][1]
    if (!is.na(low)) {
      if (!(low %in% Sc)) warning("Lower-toxicity neighbour of combination ", c, " is not in S0; restriction ignored.")
      else out <- rbind(out, c(toxic = c, lower = low))
    }
  }
  out
}

## Rules (dep, cond) for restriction "tox": the combination `dep` receives no
## patients whenever the combination `cond` receives none.
##  * W2 rule: dep in W2, cond = the T2 combination ordered immediately before it;
##  * T1 rule: dep in T1, cond = the W1 combination ordered before it, reached by
##    walking back through any consecutive T1 combinations (by definition the
##    first non-T1 combination before a T1 combination is in W1).
w2_rules <- function(s, pos_m) {
  ord <- order(pos_m)
  out <- NULL
  for (c in s$W2) {
    pred <- ord[pos_m[c] - 1]
    if (pred %in% s$T2) out <- rbind(out, c(dep = c, cond = pred))
  }
  out
}
t1_rules <- function(s, pos_m) {
  ord <- order(pos_m)
  out <- NULL
  for (c in s$T1) {
    p <- pos_m[c] - 1
    while (p >= 1 && ord[p] %in% s$T1) p <- p - 1      # skip consecutive T1 combinations
    if (p >= 1 && ord[p] %in% s$W1) out <- rbind(out, c(dep = c, cond = ord[p]))
  }
  out
}
## Rules (toxic, lower) for runs of consecutive T1 combinations under the
## incorrect ordering: within a run of length >= 2, the first combination of the
## run (`lower`) receives at least as many patients as each other one (`toxic`),
## i.e. zeta[toxic] <= zeta[lower].
t1_run_rules <- function(s, pos_m) {
  ord <- order(pos_m)
  inT1 <- ord %in% s$T1                                # by position
  out <- NULL; p <- 1
  while (p <= length(ord)) {
    if (inT1[p]) {
      q <- p; while (q < length(ord) && inT1[q + 1]) q <- q + 1
      if (q > p) for (r in (p + 1):q) out <- rbind(out, c(toxic = ord[r], lower = ord[p]))
      p <- q + 1
    } else p <- p + 1
  }
  out
}

## Build the restriction for one incorrect ordering
make_restriction <- function(restriction, s, Sc, R, pos_m, TTL) {
  switch(restriction,
         tox      = list(zero = rbind(w2_rules(s, pos_m), t1_rules(s, pos_m)),
                         le   = t1_run_rules(s, pos_m)),
         adjacent = list(le = rbind(toxic_pairs(Sc, R, pos_m, TTL), t1_run_rules(s, pos_m))),
         none     = list(),
         stop("restriction must be one of \"tox\", \"adjacent\", \"none\"."))
}
restriction_text <- function(rs, cnames) {
  txt <- c(if (!is.null(rs$le))   paste0("zeta", cnames[rs$le[, "toxic"]], " <= zeta", cnames[rs$le[, "lower"]]),
           if (!is.null(rs$zero)) paste0("zeta", cnames[rs$zero[, "cond"]], " = 0 => zeta", cnames[rs$zero[, "dep"]], " = 0"))
  if (length(txt)) paste(txt, collapse = "; ") else "none"
}

make_grid <- function(k, nu, Sc, u_grid, h = NULL, n_dir = 20, rs = list()) {
  d  <- length(Sc)
  hh <- if (is.null(h)) default_h(d) else h
  ok <- function(Z) {                                  # rows satisfying the zeta restriction
    keep <- rep(TRUE, nrow(Z))
    if (!is.null(rs$le)) {
      jt <- match(rs$le[, "toxic"], Sc); jl <- match(rs$le[, "lower"], Sc)
      keep <- keep & apply(Z, 1, function(z) all(z[jt] <= z[jl] + 1e-12))
    }
    if (!is.null(rs$zero)) {
      jw <- match(rs$zero[, "dep"], Sc); jt <- match(rs$zero[, "cond"], Sc)
      keep <- keep & apply(Z, 1, function(z) all(!(z[jt] < 1e-12 & z[jw] > 1e-12)))
    }
    keep
  }
  L <- simplex_lattice(d, hh); L <- L[ok(L), , drop = FALSE]
  Rn <- NULL                                           # random directions, by rejection
  for (tries in 1:50) {
    if (!is.null(Rn) && nrow(Rn) >= n_dir) break
    x <- matrix(rgamma(20 * max(n_dir, 1) * d, 1), 20 * max(n_dir, 1), d); x <- x / rowSums(x)
    Rn <- rbind(Rn, x[ok(x), , drop = FALSE])
  }
  Z <- rbind(L, Rn[seq_len(min(n_dir, nrow(Rn))), , drop = FALSE])
  E <- NULL
  for (u in u_grid) {
    Zu <- Z
    if (u == 0) Zu <- Z[apply(Z, 1, max) < 1 - 1e-12, , drop = FALSE]   # drop exact ties
    Eu <- matrix(0, nrow(Zu), k)
    Eu[, Sc] <- (1 - u) * Zu
    Eu[, nu] <- Eu[, nu] + u
    E <- rbind(E, Eu)
  }
  E
}

## ---------------------------------------------------------------------
## Limiting MLEs for one ordering at many allocations (rows of E).
## The score is decreasing in a, so bisection is used for all rows at once.
## Returns the fitted probabilities, an nrow(E) x k matrix.
## ---------------------------------------------------------------------
fit_at_allocations <- function(E, R, sk, lo = 1e-3, hi = 40, iter = 60) {
  ls <- log(sk)
  lo <- rep(lo, nrow(E)); hi <- rep(hi, nrow(E))
  for (it in seq_len(iter)) {
    a  <- (lo + hi) / 2
    P  <- exp(outer(a, ls))                                   # P[r, c] = sk[c]^a[r]
    sc <- rowSums(E * sweep(-P, 2, R, "+") / (1 - P) * matrix(ls, nrow(E), length(sk), byrow = TRUE))
    up <- sc > 0
    lo[up] <- a[up]; hi[!up] <- a[!up]
  }
  exp(outer((lo + hi) / 2, ls))
}

## ---------------------------------------------------------------------
## Condition 2 for one scenario.
##  skeleton  : toxicity skeleton (increasing, length k)
##  R         : true toxicities (length k)
##  orderings : M x k matrix of the INCLUDED orderings
##  ids       : labels of the orderings (e.g. c(1, 11, 36, ...)); default 1..M
##  u_grid    : shares of patients at the MTC (must be in (0, 1)); default 0.1, 0.2, ..., 0.9,
##              0.95, 0.98, 0.995
##  u_min     : optional lower bound on the MTC share (default 0 = full condition)
##  h, n_dir  : lattice resolution and number of extra random directions
##  tol       : numerical tolerance; LHS >= -tol counts as satisfied
##  restriction : "tox" (default), "adjacent" or "none"; see the header.
##  (previous argument restrict_toxic, kept for reference:) directions zeta must give each combination of S0
##              with toxicity above the TTL no more patients than the nearest combination
##              ordered before it under the incorrect ordering with toxicity below the TTL
## Returns a data frame with one row per incorrect ordering:
##  T, W, the correct orderings, the best correct ordering t, the smallest
##  LHS over the grid for that t (the global minimum, not the first failure),
##  whether the condition holds, and the allocation at the minimum.
## ---------------------------------------------------------------------
CheckPOCRM <- function(skeleton, R, orderings, ids = seq_len(nrow(orderings)), TTL = 0.3,
                       u_grid = c(seq(0.1, 0.9, by = 0.1), 0.95, 0.98, 0.995),
                       u_min = 0, h = NULL, n_dir = 20, tol = 1e-3, seed = 1,
                       cnames = NULL, restriction = "tox") {
  set.seed(seed)
  k  <- length(R)
  if (is.null(cnames)) cnames <- paste0("c", seq_len(k))
  st <- get_sets(R, orderings, TTL); nu <- st$nu
  u_grid <- u_grid[u_grid >= u_min & u_grid < 1]
  lab <- function(v) if (length(v)) paste(cnames[v], collapse = ",") else "-"
  ct <- st$correct
  out <- NULL
  for (m in setdiff(seq_len(nrow(orderings)), ct)) {
    s  <- st$sets[[as.character(m)]]
    Sc <- sort(union(s$T, s$W))
    if (!length(Sc)) next
    if (!length(ct)) {
      out <- rbind(out, data.frame(m = ids[m], T = lab(s$T), W = lab(s$W), correct = "-",
                                   best_t = NA, min_LHS = NA, holds = FALSE, worst_eta = "no correct ordering included"))
      next
    }
    rs <- make_restriction(restriction, s, Sc, R, st$pos[m, ], TTL)
    E  <- make_grid(k, nu, Sc, u_grid, h, n_dir, rs)
    Fm <- fit_at_allocations(E, R, skeleton[st$pos[m, ]])
    Dm <- f_part(Fm[, Sc, drop = FALSE], matrix(R[Sc], nrow(E), length(Sc), byrow = TRUE))
    mins <- sapply(ct, function(t) {
      Ft  <- fit_at_allocations(E, R, skeleton[st$pos[t, ]])
      Dt  <- f_part(Ft[, Sc, drop = FALSE], matrix(R[Sc], nrow(E), length(Sc), byrow = TRUE))
      lhs <- rowSums(E[, Sc, drop = FALSE] * (Dt - Dm))
      c(min(lhs), which.min(lhs))
    })
    b    <- which.max(mins[1, ])
    eta  <- E[mins[2, b], ]
    keep <- which(eta > 1e-12)
    out  <- rbind(out, data.frame(
      m = ids[m], T = lab(s$T), W = lab(s$W), correct = paste(ids[ct], collapse = ","),
      best_t = ids[ct[b]], min_LHS = signif(mins[1, b], 3), holds = mins[1, b] >= -tol,
      worst_eta = paste0(cnames[keep], "=", round(eta[keep], 3), collapse = " ")))
  }
  out
}

## Wrapper over several scenarios (rows of `scen`)
CheckPOCRM_all <- function(skeleton, scen, orderings, ids = seq_len(nrow(orderings)), ...) {
  res <- lapply(seq_len(nrow(scen)), function(s) {
    r <- CheckPOCRM(skeleton, scen[s, ], orderings, ids, ...)
    if (is.null(r)) return(NULL)
    cbind(scenario = s, r)
  })
  do.call(rbind, res)
}

## =====================================================================
## Condition 2 for ONE incorrect ordering under one scenario.
##  m          : label of the incorrect ordering to check (an element of `ids`)
##  return_grid: if TRUE, also return every allocation checked and the LHS of
##               (2) for each correct ordering (useful for plotting/diagnosis)
## Other arguments as in CheckPOCRM().
## Returns a list with
##  ordering, holds, best_t, min_LHS   : overall result
##  T, W, S0                           : the sets of Definition 3 (S0 = T ∪ W)
##  by_t                               : for EACH included correct ordering t,
##                                       the smallest LHS, whether it holds,
##                                       and the allocation where the minimum occurs
##  grid (optional)                    : allocations (rows) and LHS for each t
## =====================================================================
CheckPOCRM_one <- function(skeleton, R, orderings, m, ids = seq_len(nrow(orderings)), TTL = 0.3,
                           u_grid = c(seq(0.1, 0.9, by = 0.1), 0.95, 0.98, 0.995),
                           u_min = 0, h = NULL, n_dir = 20, tol = 1e-6, seed = 1,
                           cnames = NULL, return_grid = FALSE, restriction = "tox") {
  set.seed(seed)
  k <- length(R)
  if (is.null(cnames)) cnames <- paste0("c", seq_len(k))
  lab <- function(v) if (length(v)) paste(cnames[v], collapse = ",") else "-"
  row <- match(m, ids)
  if (is.na(row)) stop("Ordering ", m, " is not among the included orderings (ids).")
  st <- get_sets(R, orderings, TTL); nu <- st$nu; ct <- st$correct
  if (row %in% ct) {
    message("Ordering ", m, " belongs to the correct group under this scenario; Condition 2 does not apply.")
    return(invisible(NULL))
  }
  s  <- st$sets[[as.character(row)]]
  Sc <- sort(union(s$T, s$W))
  res <- list(ordering = m, MTC = cnames[nu], T = lab(s$T), W = lab(s$W), S0 = lab(Sc),
              correct = ids[ct])
  if (!length(Sc)) {
    message("S0 is empty for ordering ", m, ": nothing to check (condition holds trivially).")
    return(c(res, list(holds = TRUE)))
  }
  if (!length(ct)) {
    message("No ordering from the correct group is included: Condition 1 fails.")
    return(c(res, list(holds = FALSE)))
  }
  u_grid <- u_grid[u_grid >= u_min & u_grid < 1]
  rs <- make_restriction(restriction, s, Sc, R, st$pos[row, ], TTL)
  E  <- make_grid(k, nu, Sc, u_grid, h, n_dir, rs)
  Rm <- matrix(R[Sc], nrow(E), length(Sc), byrow = TRUE)
  Dm <- f_part(fit_at_allocations(E, R, skeleton[st$pos[row, ]])[, Sc, drop = FALSE], Rm)
  LHS <- sapply(ct, function(t) {
    Dt <- f_part(fit_at_allocations(E, R, skeleton[st$pos[t, ]])[, Sc, drop = FALSE], Rm)
    rowSums(E[, Sc, drop = FALSE] * (Dt - Dm))
  })
  LHS <- matrix(LHS, nrow(E), length(ct)); colnames(LHS) <- paste0("t", ids[ct])
  fmt_eta <- function(e) { kk <- which(e > 1e-12); paste0(cnames[kk], "=", round(e[kk], 3), collapse = " ") }
  by_t <- data.frame(t = ids[ct], min_LHS = signif(apply(LHS, 2, min), 3),
                     holds = apply(LHS, 2, min) >= -tol,
                     worst_eta = apply(LHS, 2, function(x) fmt_eta(E[which.min(x), ])))
  b <- which.max(by_t$min_LHS)
  restr <- restriction_text(rs, cnames)
  res <- c(res, list(restriction = restr, holds = any(by_t$holds), best_t = by_t$t[b], min_LHS = by_t$min_LHS[b],
                     n_allocations = nrow(E), by_t = by_t))
  if (return_grid) {
    colnames(E) <- cnames
    res$grid <- data.frame(u_MTC = E[, nu], E[, Sc, drop = FALSE], LHS, check.names = FALSE)
  }
  res
}