## =====================================================================
## Check of Condition 3 of Theorem 2: the CRM consistency condition
## (Lemma 1, Cheung and Chappell 2002) under every INCLUDED ordering that
## belongs to the correct group.
##
## For an ordering t with the MTC in position l, and skeleton values
## alpha(1) < ... < alpha(k) assigned by position:
##   a_i = log R_t(i) / log alpha(i), where R_t(i) is the true toxicity of
##         the combination in position i of ordering t;
##   b_i solves alpha(i-1)^b + alpha(i)^b = 2 * TTL, i = 2, ..., k,
##         with b_1 = -Inf and b_{k+1} = +Inf;
##   H_i = (b_i, b_{i+1}).
## Lemma 1 requires
##   (MTC)        a_l in H_l, i.e. b_l < a_l < b_{l+1};
##   (below MTC)  a_i > b_{i+1}  for i = 1, ..., l-1;
##   (above MTC)  a_i < b_i      for i = l+1, ..., k.
##
## Conventions (as in check_condition2_only.R)
##  * Combinations are indexed 1..k; a scenario R is a length-k vector.
##  * An ordering is a length-k vector of combination indices, from least
##    to most toxic; `orderings` is an M x k matrix of the INCLUDED orderings.
##  * The skeleton `skel` is increasing, of length k.
##  * Single-MTC scenarios (exactly one combination with toxicity == TTL).
## =====================================================================

## boundaries b_2, ..., b_k of Lemma 1
crm_boundaries <- function(skel, TTL = 0.3) {
  sapply(2:length(skel), function(i)
    uniroot(function(b) skel[i - 1]^b + skel[i]^b - 2 * TTL,
            lower = 1e-8, upper = 1e3, tol = 1e-12)$root)
}

## correct group (Definition 2) among the given orderings
correct_orderings <- function(R, orderings, TTL = 0.3) {
  nu <- which(abs(R - TTL) < 1e-10)
  if (length(nu) != 1) stop("Scenario must have exactly one MTC (toxicity == TTL).")
  below <- which(R < R[nu]); nb <- length(below)
  which(apply(orderings, 1, function(o) o[nb + 1] == nu && setequal(o[seq_len(nb)], below)))
}

## ---------------------------------------------------------------------
## Lemma 1 for ONE ordering.
## Returns a list with: holds, MTC position, the a_i's, the intervals, and
## a data frame of violated requirements (position, combination, a_i,
## required bound).
## ---------------------------------------------------------------------
CheckCRM <- function(R, alpha, ordering, TTL = 0.3, cnames = NULL, digits = 3) {
  k <- length(alpha)
  if (is.null(cnames)) cnames <- paste0("c", seq_len(k))
  Rt <- R[ordering]                                  # true toxicity by position
  l  <- which(abs(Rt - TTL) < 1e-10)
  if (length(l) != 1) stop("Scenario must have exactly one MTC (toxicity == TTL).")
  a  <- log(Rt) / log(alpha)
  b  <- c(-Inf, crm_boundaries(alpha, TTL), Inf)      # b[i] = b_i, i = 1, ..., k+1
  viol <- NULL
  add <- function(i, type, bound) viol <<- rbind(viol, data.frame(
    position = i, combination = cnames[ordering[i]], R = Rt[i], a = round(a[i], digits),
    requirement = type, bound = round(bound, digits)))
  if (!(a[l] > b[l]))     add(l, "MTC: a_l > b_l", b[l])
  if (!(a[l] < b[l + 1])) add(l, "MTC: a_l < b_{l+1}", b[l + 1])
  for (i in seq_len(l - 1))       if (!(a[i] > b[i + 1])) add(i, "below MTC: a_i > b_{i+1}", b[i + 1])
  if (l < k) for (i in (l + 1):k) if (!(a[i] < b[i]))     add(i, "above MTC: a_i < b_i", b[i])
  return(
    list(holds = is.null(viol), MTC_position = l,
       a = setNames(round(a, digits), cnames[ordering]),
       H = data.frame(i = seq_len(k), lower = round(b[1:k], digits), upper = round(b[2:(k + 1)], digits)),
       violations = viol)
    )
}

## ---------------------------------------------------------------------
## Condition 3 for one scenario: Lemma 1 under every included correct ordering.
##  ids    : labels of the included orderings (e.g. c(1, 11, 36, ...))
##  cnames : labels of the combinations (e.g. "(4,1)")
## Returns a data frame with one row per included correct ordering, and the
## violations (if any) as an attribute "violations".
## ---------------------------------------------------------------------
CheckCRM_all <- function(R, alpha, orderings, ids = seq_len(nrow(orderings)), TTL = 0.3,
                             cnames = NULL) {
  ct <- correct_orderings(R, orderings, TTL)
  if (!length(ct)) {
    warning("No included ordering belongs to the correct group: Condition 1 fails.")
    return(NULL)
  }
  res <- lapply(ct, function(t) CheckCRM(R, alpha, orderings[t, ], TTL, cnames))
  out <- data.frame(t = ids[ct], MTC_position = sapply(res, `[[`, "MTC_position"),
                    holds = sapply(res, `[[`, "holds"),
                    n_violations = sapply(res, function(r) if (is.null(r$violations)) 0L else nrow(r$violations)))
  v <- do.call(rbind, lapply(seq_along(ct), function(j)
    if (!is.null(res[[j]]$violations)) cbind(t = ids[ct[j]], res[[j]]$violations)))
  attr(out, "violations") <- v
  out
}