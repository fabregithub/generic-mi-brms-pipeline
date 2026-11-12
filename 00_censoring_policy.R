# ---------------------------------------------------------------------------
# 00_censoring_policy.R  --  turn raw per-participant limits into leftcens bounds
#                            under an explicitly NAMED censoring policy.
#
# Why this file exists
# --------------------
# A lab delivers, per analyte and per participant: an interval (lo, hi) and two
# limits, MDL and LCMRL. Which limit a value is "censored at" is NOT a property
# of the data -- it is an analysis choice, and different arms of the same study
# make different ones. Itoh et al. substitute limit/2; our pipeline keeps the
# interval. Freezing that choice at extraction time would make the two arms
# indistinguishable, which is the whole contrast.
#
# So: the extraction carries the FINEST bounds the lab can give, and this file
# applies a policy on top. The engine downstream
# (00_censored_exposure.R -> leftcens::impute_censored_conditional) already
# consumes two-sided intervals natively; nothing there needs to change.
#
# The one-way street
# ------------------
#   three_tier  (0,MDL) | (MDL,LCMRL) | point      <- finest
#        |  collapse
#        v
#   lcmrl       (0,LCMRL) | point                  <- coarser, DNQ band lost
#        |  collapse
#        v
#   substitute  point at limit/divisor everywhere  <- coarsest (Itoh's arm A)
#
# substitute(limit = "mdl") is a HYBRID, not a fourth level: sub-MDL rows become
# a point at MDL/2, while DNQ rows stay intervals, because a DNQ value is known
# to lie above the MDL and MDL/2 would fall below its own lower bound.
#
# Coarsening is always possible; refining never is. ce_apply_policy() therefore
# REFUSES to return bounds finer than what it was given, rather than inventing
# the missing distinction.
#
# Written 2026-11-11. Self-test at the bottom: Rscript 00_censoring_policy.R
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

CE_POLICIES <- c("three_tier", "lcmrl", "substitute")

#' Describe what censoring structure a dataset's bounds actually carry.
#'
#' Call this on the real extraction BEFORE choosing a policy. It answers the
#' question that decides whether arm B can exploit the DNQ band at all.
ce_describe_bounds <- function(df, analytes,
                               lo = "_lo", hi = "_hi",
                               mdl = "_mdl", lcmrl = "_lcmrl",
                               tol = 1e-8) {
  near <- function(a, b) abs(a - b) <= tol + 1e-6 * pmax(abs(a), abs(b))
  out <- lapply(analytes, function(a) {
    L <- df[[paste0(a, lo)]]; H <- df[[paste0(a, hi)]]
    M <- df[[paste0(a, mdl)]]; C <- df[[paste0(a, lcmrl)]]
    point <- near(L, H)
    nd    <- !point & L <= tol & near(H, M)
    dnq   <- !point & near(L, M) & near(H, C)
    coll  <- !point & L <= tol & near(H, C)      # (0, LCMRL): already collapsed
    other <- !(point | nd | dnq | coll)
    data.frame(analyte = a, n = length(L),
               quantified = sum(point), nondetect_mdl = sum(nd),
               dnq = sum(dnq), collapsed_lcmrl = sum(coll), unclassified = sum(other),
               finest = if (sum(dnq) > 0) "three_tier" else if (sum(coll) > 0 || sum(nd) > 0)
                          (if (sum(coll) > 0) "lcmrl" else "three_tier") else "all_quantified",
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  attr(out, "overall_finest") <-
    if (any(out$dnq > 0)) "three_tier" else if (any(out$collapsed_lcmrl > 0)) "lcmrl" else "all_quantified"
  out
}

#' Apply a named censoring policy, rewriting <A>_lo / <A>_hi in place and
#' adding a derived <A>_censtype.
#'
#' @param policy "three_tier" keeps the finest structure (censor at MDL, DNQ as
#'   a bounded interval); "lcmrl" treats everything below LCMRL as one
#'   non-detect; "substitute" replaces every censored row with a point value at
#'   limit/divisor (Itoh's arm A).
#' @param limit which limit "lcmrl"/"substitute" use. Itoh say "MRL" -- confirm
#'   against their section 2.2 before running arm A.
ce_apply_policy <- function(df, analytes,
                            policy = c("three_tier", "lcmrl", "substitute"),
                            limit = c("lcmrl", "mdl"), divisor = 2,
                            lo = "_lo", hi = "_hi", mdl = "_mdl", lcmrl = "_lcmrl",
                            tol = 1e-8, strict = TRUE) {
  policy <- match.arg(policy); limit <- match.arg(limit)
  near <- function(a, b) abs(a - b) <= tol + 1e-6 * pmax(abs(a), abs(b))

  have <- ce_describe_bounds(df, analytes, lo, hi, mdl, lcmrl, tol)
  if (strict && policy == "three_tier" && attr(have, "overall_finest") == "lcmrl")
    stop("ce_apply_policy(): asked for 'three_tier' but the bounds supplied are already ",
         "collapsed to (0, LCMRL) for at least one analyte. The MDL/LCMRL distinction was ",
         "discarded upstream and cannot be recovered here. Either obtain finer bounds or ",
         "use policy = 'lcmrl'.", call. = FALSE)
  if (any(have$unclassified > 0))
    warning("ce_apply_policy(): ", sum(have$unclassified), " row(s) have bounds matching no ",
            "recognised pattern. Inspect with ce_describe_bounds() before trusting the result.",
            call. = FALSE)

  for (a in analytes) {
    cl <- paste0(a, c(lo, hi, mdl, lcmrl))
    if (!all(cl %in% names(df))) stop("ce_apply_policy(): missing column(s) for '", a, "'", call. = FALSE)
    L <- df[[cl[1]]]; H <- df[[cl[2]]]; M <- df[[cl[3]]]; C <- df[[cl[4]]]
    if (anyNA(M) || anyNA(C))
      stop("ce_apply_policy(): '", a, "' has missing MDL/LCMRL. Both must be populated for ",
           "every participant, censored or not.", call. = FALSE)
    if (any(C < M - tol)) stop("ce_apply_policy(): '", a, "' has LCMRL < MDL", call. = FALSE)
    if (any(H < L - tol)) stop("ce_apply_policy(): '", a, "' has hi < lo", call. = FALSE)

    point <- near(L, H)
    thr   <- if (limit == "lcmrl") C else M

    if (policy == "three_tier") {
      newL <- L; newH <- H                                  # pass through unchanged
    } else if (policy == "lcmrl") {
      keep <- point & (L >= C - tol)                        # quantified at or above LCMRL
      newL <- ifelse(keep, L, 0)
      newH <- ifelse(keep, H, C)
    } else {                                                 # substitute
      sub <- thr / divisor
      if (limit == "lcmrl") {
        keep <- point & (L >= C - tol)                      # Itoh's arm A: everything
        newL <- ifelse(keep, L, sub)                        # below LCMRL -> LCMRL/2,
        newH <- ifelse(keep, H, sub)                        # DNQ band included
      } else {
        # Substituting at the MDL leaves the DNQ band undefined: those rows are
        # known to lie ABOVE the MDL, so MDL/2 would sit below their own lower
        # bound. The only coherent reading is to substitute the sub-MDL rows and
        # leave DNQ as the interval it is. The result is therefore a hybrid, not
        # a fully substituted dataset -- which is itself worth reporting.
        isnd <- !point & L <= tol
        newL <- ifelse(isnd, sub, L)
        newH <- ifelse(isnd, sub, H)
      }
    }

    df[[cl[1]]] <- newL; df[[cl[2]]] <- newH
    p2 <- near(newL, newH)
    df[[paste0(a, "_censtype")]] <-
      ifelse(p2 & !point, "substituted",                    # was an interval, now a point
      ifelse(p2, "quantified",
      ifelse(newL <= tol, "nondetect", "dnq")))
    df[[a]] <- ifelse(p2, newL, NA_real_)                   # engine reads bounds; this is for QC
  }
  attr(df, "censoring_policy") <- list(policy = policy, limit = limit,
                                       divisor = if (policy == "substitute") divisor else NA_real_,
                                       applied = Sys.time())
  df
}

# ---------------------------------------------------------------------------
# Self-test. Verifies the three policies against values constructed so the
# right answer is known, and verifies that refinement is refused.
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  set.seed(11)
  n <- 2000; A <- "pfoa"
  M <- runif(n, 0.05, 0.15); C <- M * runif(n, 1.5, 3)
  x <- exp(rnorm(n, log(0.12), 1.1))
  nd <- x < M; dq <- x >= M & x < C
  d <- data.frame(pfoa_lo = ifelse(nd, 0, ifelse(dq, M, x)),
                  pfoa_hi = ifelse(nd, M, ifelse(dq, C, x)),
                  pfoa_mdl = M, pfoa_lcmrl = C)
  fails <- 0L
  chk <- function(ok, msg) { if (!ok) { fails <<- fails + 1L; cat("  FAIL  ", msg, "\n") } else cat("  ok    ", msg, "\n") }

  s <- ce_describe_bounds(d, A)
  chk(attr(s, "overall_finest") == "three_tier", "describe: recognises three-tier input")
  chk(s$unclassified == 0, "describe: every row classified")
  chk(s$quantified + s$nondetect_mdl + s$dnq == n, "describe: states partition the rows")

  t3 <- ce_apply_policy(d, A, "three_tier")
  chk(identical(t3$pfoa_lo, d$pfoa_lo) && identical(t3$pfoa_hi, d$pfoa_hi),
      "three_tier: bounds pass through unchanged")
  chk(sum(t3$pfoa_censtype == "dnq") == sum(dq), "three_tier: DNQ band retained")
  chk(all(t3$pfoa_lo[t3$pfoa_censtype == "nondetect"] == 0), "three_tier: ND lower bound is 0")

  tl <- ce_apply_policy(d, A, "lcmrl")
  chk(sum(tl$pfoa_censtype == "dnq") == 0, "lcmrl: DNQ band collapsed away")
  chk(sum(tl$pfoa_censtype == "nondetect") == sum(nd | dq), "lcmrl: ND count = ND + DNQ")
  chk(all(tl$pfoa_hi[tl$pfoa_censtype == "nondetect"] == C[nd | dq]), "lcmrl: upper bound is LCMRL")
  chk(all(tl$pfoa_lo[tl$pfoa_censtype == "quantified"] == x[!nd & !dq]), "lcmrl: quantified untouched")

  ts <- ce_apply_policy(d, A, "substitute", limit = "lcmrl", divisor = 2)
  chk(all(ts$pfoa[nd | dq] == (C / 2)[nd | dq]), "substitute: censored rows set to LCMRL/2")
  chk(all(ts$pfoa_lo == ts$pfoa_hi), "substitute: every row is a point value")
  chk(sum(ts$pfoa_censtype == "substituted") == sum(nd | dq), "substitute: censtype records it")

  tm <- ce_apply_policy(d, A, "substitute", limit = "mdl", divisor = 2)
  chk(sum(tm$pfoa_censtype == "substituted") == sum(nd),
      "substitute(mdl): only sub-MDL rows substituted, DNQ kept as values")
  chk(!isTRUE(all.equal(sum(ts$pfoa_censtype == "substituted"),
                        sum(tm$pfoa_censtype == "substituted"))),
      "substitute: the choice of limit changes the result (so it must be declared)")

  coarse <- ce_apply_policy(d, A, "lcmrl")
  r <- tryCatch({ ce_apply_policy(coarse, A, "three_tier"); FALSE },
                error = function(e) grepl("cannot be recovered", conditionMessage(e)))
  chk(isTRUE(r), "refinement from collapsed bounds is refused, not faked")

  cat(sprintf("\n%d failure(s)\n", fails))
  quit(status = if (fails > 0L) 1L else 0L)
}
