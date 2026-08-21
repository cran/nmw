# ctl2nmw: auto-generate the nmw PRED() + initial-estimate block from a
# closed-form NONMEM $PRED control stream.  Reads NONMEM *input* (control
# stream); cf. nm_parse.R which reads NONMEM *output*.
#
# Scope (v1): a single $PROBLEM with a closed-form $PRED block (F = <algebra>,
# Y = <error model in F and EPS>).  ODE models ($PK/$DES/ADVAN/$SUBROUTINE),
# IF/THEN branching, MU-referencing, mixtures, and FIXed parameters are out of
# scope and raise an informative error rather than silently mis-generating.

## ------------------------------------------------------------------ helpers

## Functions R's stats::deriv() can differentiate (its derivative table).
.nm_deriv_funcs <- c("exp", "log", "log10", "sqrt", "sin", "cos", "tan",
                     "sinh", "cosh", "asin", "acos", "atan",
                     "pnorm", "dnorm", "gamma", "lgamma", "digamma",
                     "trigamma", "psigamma")
.nm_ok_calls <- c("+", "-", "*", "/", "^", "(", .nm_deriv_funcs)

.nm_read_lines <- function(ctl) {
  if (is.character(ctl) && length(ctl) == 1 && file.exists(ctl))
    ctl <- readLines(ctl, warn = FALSE)
  unlist(strsplit(as.character(ctl), "\r?\n"))   # also splits embedded newlines
}

## Canonicalize a record tag (uppercase, no '$') to a short label.
.nm_canon_tag <- function(tag) {
  t <- toupper(tag)
  exact <- c(PROB = "PROB", PROBLEM = "PROB", PRO = "PROB",
             INPUT = "INPUT", INP = "INPUT", INPT = "INPUT",
             DATA = "DATA", DAT = "DATA",
             PRED = "PRED", PRE = "PRED",
             THETA = "THETA", THE = "THETA", THETAS = "THETA",
             OMEGA = "OMEGA", OME = "OMEGA", OMEGAS = "OMEGA",
             SIGMA = "SIGMA", SIG = "SIGMA", SIGMAS = "SIGMA",
             EST = "EST", ESTM = "EST", ESTIMATE = "EST", ESTIMATION = "EST",
             ESTIM = "EST", ESTIMAT = "EST",
             COV = "COV", COVARIANCE = "COV", COVR = "COV",
             TABLE = "TABLE", TAB = "TABLE",
             SUBROUTINE = "SUB", SUBROUTINES = "SUB", SUBS = "SUB", SUB = "SUB",
             MODEL = "MODEL", MOD = "MODEL",
             PK = "PK", DES = "DES", ERROR = "ERROR", ERR = "ERROR",
             MIX = "MIX", MIXTURE = "MIX", PRIOR = "PRIOR",
             SIM = "SIM", SIMULATION = "SIM", SIMULATE = "SIM",
             MSFI = "MSFI", ABBREVIATED = "ABBR", ABBR = "ABBR",
             AES = "AES", AESINITIAL = "AESI", INFN = "INFN", CONTR = "CONTR")
  if (!is.na(exact[t])) return(unname(exact[t]))
  t
}

## Split control-stream lines into an ordered list of {tag, body, line}.
.nm_split_records <- function(lines) {
  clean <- gsub("\t", " ", sub(";.*$", "", lines))       # strip comments/tabs
  is_start <- grepl("^\\s*\\$", clean)
  recs <- list()
  cur <- NULL; body <- character(0); startline <- NA_integer_
  flush <- function() {
    if (!is.null(cur))
      recs[[length(recs) + 1L]] <<- list(tag = cur, body = body, line = startline)
  }
  for (i in seq_along(clean)) {
    ln <- clean[i]
    if (is_start[i]) {
      flush()
      raw_tag <- sub("^\\s*\\$([A-Za-z0-9]+).*$", "\\1", ln)
      cur <- .nm_canon_tag(raw_tag)
      body <- sub("^\\s*\\$[A-Za-z0-9]+[ \t]*", "", ln)   # remainder on tag line
      startline <- i
    } else if (!is.null(cur)) {
      body <- c(body, ln)
    }
  }
  flush()
  recs
}

.nm_recs_of <- function(recs, tag) {
  lapply(Filter(function(r) r$tag == tag, recs), function(r) r$body)
}

## NONMEM abbreviated code -> R source string (one RHS at a time).
.nm_translate <- function(s) {
  s <- sub(";.*$", "", s)
  s <- gsub("THETA\\s*\\(\\s*([0-9]+)\\s*\\)", "TH\\1", s, ignore.case = TRUE)
  s <- gsub("\\bETA\\s*\\(\\s*([0-9]+)\\s*\\)", "ETA\\1", s, ignore.case = TRUE)
  s <- gsub("\\b(EPS|ERR)\\s*\\(\\s*([0-9]+)\\s*\\)", "EPS\\2", s, ignore.case = TRUE)
  s <- gsub("\\*\\*", "^", s)
  for (f in c("LOG10", "LOG", "EXP", "SQRT", "SINH", "COSH", "SIN", "COS",
              "TAN", "ASIN", "ACOS", "ATAN", "GAMMA", "LGAMMA")) {
    s <- gsub(paste0("\\b", f, "\\b"), tolower(f), s, ignore.case = TRUE)
  }
  s
}

.nm_num <- function(txt) {
  x <- as.numeric(regmatches(txt,
        gregexpr("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", txt))[[1]])
  x[!is.na(x)]
}

.nm_has_fix <- function(bodies) {
  any(grepl("\\bFIX(ED)?\\b", toupper(unlist(bodies))))
}

## ------------------------------------------------------------- $INPUT

.nm_parse_input <- function(bodies) {
  toks <- unlist(strsplit(trimws(paste(unlist(bodies), collapse = " ")), "[ ,]+"))
  toks <- toks[nzchar(toks)]
  reserved <- c("ID", "TIME", "DV", "AMT", "MDV", "EVID", "CMT", "RATE",
                "SS", "II", "ADDL", "DATE", "DAT1", "DAT2", "DAT3")
  all <- character(0); drop <- logical(0); syn <- list()
  for (tk in toks) {
    if (grepl("=", tk)) {
      lr <- strsplit(tk, "=")[[1]]
      L <- toupper(lr[1]); R <- toupper(lr[2])
      if (R %in% c("DROP", "SKIP")) {
        all <- c(all, L); drop <- c(drop, TRUE)
      } else {
        eff <- if (R %in% reserved) R else L
        all <- c(all, eff); drop <- c(drop, FALSE)
        syn[[length(syn) + 1L]] <- c(L, R)
      }
    } else {
      all <- c(all, toupper(tk)); drop <- c(drop, FALSE)
    }
  }
  list(all = all, used = all[!drop], drop = drop, syn = syn)
}

## ------------------------------------------------------------- $THETA

.nm_parse_theta <- function(bodies) {
  txt <- paste(unlist(bodies), collapse = " ")
  toks <- regmatches(txt,
    gregexpr("\\([^)]*\\)|[Xx][0-9]+|[A-Za-z]+|[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?",
             txt))[[1]]
  init <- lb <- ub <- numeric(0)
  widened <- integer(0)                        # unbounded thetas whose LB became -1e6
  clampInf <- function(v, sign) if (is.na(v)) sign * 1e6 else v
  i <- 1L
  while (i <= length(toks)) {
    tk <- toks[i]
    spec <- NULL; unb <- FALSE
    if (grepl("^\\(", tk)) {
      parts <- trimws(strsplit(gsub("[()]", "", tk), ",")[[1]])
      nums <- suppressWarnings(as.numeric(parts))
      # map INF-like tokens
      for (k in seq_along(parts)) if (is.na(nums[k]) && grepl("INF", toupper(parts[k])))
        nums[k] <- if (grepl("-", parts[k])) -1e6 else 1e6
      nums <- nums[!is.na(nums)]
      if (length(nums) == 1) {                 # "(init)": unbounded in NONMEM
        spec <- c(if (nums[1] > 0) 0 else -1e6, nums[1], 1e6); unb <- nums[1] <= 0
      }
      else if (length(nums) == 2) spec <- c(nums[1], nums[2], 1e6)
      else if (length(nums) >= 3) spec <- c(nums[1], nums[2], nums[3])
    } else if (grepl("^[-+]?[0-9.]", tk)) {    # bare init: unbounded in NONMEM
      iv <- as.numeric(tk)
      spec <- c(if (iv > 0) 0 else -1e6, iv, 1e6); unb <- iv <= 0
    } else { i <- i + 1L; next }               # stray keyword
    reps <- 1L
    if (i + 1L <= length(toks) && grepl("^[Xx][0-9]+$", toks[i + 1L])) {
      reps <- as.integer(sub("^[Xx]", "", toks[i + 1L])); i <- i + 1L
    }
    for (r in seq_len(reps)) {
      lb <- c(lb, spec[1]); init <- c(init, spec[2]); ub <- c(ub, spec[3])
      if (unb) widened <- c(widened, length(init))
    }
    i <- i + 1L
  }
  list(init = init, lb = lb, ub = ub, widened = widened)
}

## ------------------------------------------------------- $OMEGA / $SIGMA

.nm_parse_block <- function(body) {
  up <- toupper(paste(body, collapse = " "))
  if (grepl("\\bSAME\\b", up)) return(list(kind = "SAME", n = NA, vals = numeric(0)))
  clean <- gsub("BLOCK\\s*\\([0-9]+\\)|DIAGONAL\\s*\\([0-9]+\\)|\\bFIX(ED)?\\b|\\bSAME\\b",
                " ", up)
  nums <- .nm_num(clean)
  mB <- regmatches(up, regexpr("BLOCK\\s*\\(\\s*[0-9]+\\s*\\)", up))
  mD <- regmatches(up, regexpr("DIAGONAL\\s*\\(\\s*[0-9]+\\s*\\)", up))
  if (length(mB) > 0)
    list(kind = "BLOCK", n = as.integer(gsub("[^0-9]", "", mB)), vals = nums)
  else if (length(mD) > 0)
    list(kind = "DIAG", n = as.integer(gsub("[^0-9]", "", mD)), vals = nums)
  else
    list(kind = "DIAG", n = length(nums), vals = nums)
}

.nm_parse_omega <- function(bodies) {
  blocks <- lapply(bodies, .nm_parse_block)
  ns <- vapply(blocks, function(b) if (identical(b$kind, "SAME")) NA_integer_ else b$n, 0L)
  # resolve SAME dims from previous
  prev_n <- NA_integer_
  for (k in seq_along(blocks)) {
    if (identical(blocks[[k]]$kind, "SAME")) { blocks[[k]]$n <- prev_n }
    prev_n <- blocks[[k]]$n
  }
  n_eta <- sum(vapply(blocks, function(b) b$n, 0L))
  M <- matrix(0, n_eta, n_eta)
  off <- 0L; prev <- NULL
  for (b in blocks) {
    if (identical(b$kind, "SAME")) b <- prev
    idx <- off + seq_len(b$n)
    if (b$kind == "DIAG") {
      for (r in seq_len(b$n)) M[idx[r], idx[r]] <- b$vals[r]
    } else {                                    # BLOCK: lower triangle, row-major
      k <- 1L
      for (r in seq_len(b$n)) for (c in seq_len(r)) {
        M[idx[r], idx[c]] <- M[idx[c], idx[r]] <- b$vals[k]; k <- k + 1L
      }
    }
    off <- off + b$n; prev <- b
  }
  single_full <- length(blocks) == 1 && blocks[[1]]$kind == "BLOCK"
  list(mat = M, n_eta = n_eta, densified = !single_full)
}

.nm_parse_sigma <- function(bodies) {
  blocks <- lapply(bodies, .nm_parse_block)
  diag_vals <- numeric(0); had_offdiag <- FALSE
  for (b in blocks) {
    if (b$kind == "BLOCK") {
      k <- 1L
      for (r in seq_len(b$n)) for (c in seq_len(r)) {
        if (r == c) diag_vals <- c(diag_vals, b$vals[k])
        else if (!isTRUE(all.equal(b$vals[k], 0))) had_offdiag <- TRUE
        k <- k + 1L
      }
    } else {
      diag_vals <- c(diag_vals, b$vals)
    }
  }
  list(mat = diag(diag_vals, nrow = length(diag_vals)),
       n_eps = length(diag_vals), had_offdiag = had_offdiag)
}

## ------------------------------------------------------------- $EST

.nm_parse_est <- function(bodies) {
  if (length(bodies) == 0) return(list(method = "ZERO", inter = FALSE, posthoc = FALSE))
  up <- toupper(paste(unlist(bodies[[length(bodies)]]), collapse = " "))
  m <- regmatches(up, regexpr("METH(OD)?\\s*=\\s*[A-Z0-9]+", up))
  mval <- if (length(m) > 0) sub(".*=\\s*", "", m) else ""
  laplace <- grepl("\\bLAPL", up) || grepl("\\bLAPLAC", up)
  inter   <- grepl("\\bINTER", up)
  posthoc <- grepl("\\bPOSTHOC\\b", up)
  is_cond <- mval %in% c("COND", "CONDITIONAL", "1")
  is_zero <- mval %in% c("ZERO", "FO", "0", "")
  known   <- mval %in% c("COND", "CONDITIONAL", "1", "ZERO", "FO", "0", "")
  if (!known)
    stop(sprintf("ctl2nmw: $EST METHOD=%s is not supported (v1 handles FO/FOCE(I)/LAPL only).",
                 mval), call. = FALSE)
  method <- if (laplace) "LAPL" else if (is_cond) "COND" else "ZERO"
  list(method = method, inter = inter, posthoc = posthoc)
}

## ------------------------------------------------------------- ParseCtl

ParseCtl <- function(ctl) {
  lines <- .nm_read_lines(ctl)
  recs  <- .nm_split_records(lines)
  tags  <- vapply(recs, function(r) r$tag, "")

  ## scope guards ---------------------------------------------------------
  forbidden <- intersect(tags, c("SUB", "MODEL", "PK", "DES", "MIX", "PRIOR"))
  if (length(forbidden) > 0)
    stop(sprintf(paste0("ctl2nmw: $%s found. ODE/advanced models have no closed-form ",
                        "prediction and are out of scope (v1 supports closed-form $PRED)."),
                 paste(forbidden, collapse = ", $")), call. = FALSE)
  if (!"PRED" %in% tags)
    stop("ctl2nmw: no $PRED record found (v1 requires a closed-form $PRED block).",
         call. = FALSE)
  if (.nm_has_fix(.nm_recs_of(recs, "THETA")) ||
      .nm_has_fix(.nm_recs_of(recs, "OMEGA")) ||
      .nm_has_fix(.nm_recs_of(recs, "SIGMA")))
    stop("ctl2nmw: FIXed parameters are not supported in v1 (nmw estimates every parameter).",
         call. = FALSE)

  pred_lines <- unlist(.nm_recs_of(recs, "PRED"))
  bad <- grepl("\\bIF\\b|\\bTHEN\\b|\\bELSE\\b|\\bENDIF\\b|\\bWHILE\\b|\\.EQ\\.|\\.NE\\.|\\.GT\\.|\\.LT\\.|\\.GE\\.|\\.LE\\.|\\bMU_",
               toupper(pred_lines))
  if (any(bad))
    stop(sprintf("ctl2nmw: conditional/MU-referencing logic in $PRED is unsupported (line: %s).",
                 trimws(pred_lines[which(bad)[1]])), call. = FALSE)

  input <- .nm_parse_input(.nm_recs_of(recs, "INPUT"))
  th    <- .nm_parse_theta(.nm_recs_of(recs, "THETA"))
  om    <- .nm_parse_omega(.nm_recs_of(recs, "OMEGA"))
  sg    <- .nm_parse_sigma(.nm_recs_of(recs, "SIGMA"))
  est   <- .nm_parse_est(.nm_recs_of(recs, "EST"))

  ## nmw's InitStep() scales theta with a logit-type transform, which needs
  ## LB < init < UB strictly; a violating spec would yield NaN alpha silently.
  bad <- which(!(th$lb < th$init & th$init < th$ub))
  if (length(bad) > 0)
    stop(sprintf(paste0("ctl2nmw: $THETA %d initial estimate (%g) is not strictly inside ",
                        "its bounds (%g, %g); nmw's scaling requires LB < init < UB. ",
                        "Give $THETA bounds that strictly bracket the initial estimate."),
                 bad[1], th$init[bad[1]], th$lb[bad[1]], th$ub[bad[1]]), call. = FALSE)

  prob  <- .nm_recs_of(recs, "PROB")
  notes <- character(0)
  if (om$densified)
    notes <- c(notes, paste("$OMEGA is diagonal/block-diagonal; nmw estimates a FULL block,",
                            "so degrees of freedom and OFV differ from the NONMEM original."))
  if (isTRUE(sg$had_offdiag))
    notes <- c(notes, "$SIGMA has off-diagonal terms; nmw uses the diagonal only (correlated residual error is dropped).")
  if (length(th$widened) > 0)
    notes <- c(notes, sprintf(paste("$THETA %s: unbounded in the control stream with a non-positive",
                                    "initial estimate; LB was set to -1e6 (nmw's usual default is 0)",
                                    "so that LB < init < UB holds."),
                              paste(th$widened, collapse = ", ")))
  if (est$method != "ZERO" && !est$inter)
    notes <- c(notes, "$EST is COND/LAPL without INTERACTION; nmw always applies interaction, so the objective differs from a non-INTER FOCE run.")
  if (est$posthoc)
    notes <- c(notes, "$EST POSTHOC: call PostHocEta() after EstStep() to obtain post-hoc ETAs (EstStep does not run it automatically).")

  structure(list(
    prob_title = if (length(prob)) trimws(paste(unlist(prob[[1]]), collapse = " ")) else "",
    input = input,
    n_theta = length(th$init), theta_init = th$init, lb = th$lb, ub = th$ub,
    n_eta = om$n_eta, om_init = om$mat, om_densified = om$densified,
    n_eps = sg$n_eps, sg_init = sg$mat,
    method = est$method, inter = est$inter, posthoc = est$posthoc,
    pred_lines = pred_lines, notes = notes
  ), class = "nmctl")
}

## ------------------------------------------------ $PRED math + codegen

.nm_inline <- function(expr, defs) {
  if (length(defs) == 0) return(expr)
  repeat {
    new <- do.call("substitute", list(expr, defs))
    if (identical(new, expr)) break
    expr <- new
  }
  expr
}

.nm_idx <- function(syms, prefix)
  sort(as.integer(sub(paste0("^", prefix), "",
       syms[grepl(paste0("^", prefix, "[0-9]+$"), syms)])))

.nm_check_funcs <- function(expr, where) {
  calls <- setdiff(all.names(expr), all.vars(expr))
  bad <- setdiff(calls, .nm_ok_calls)
  if (length(bad) > 0)
    stop(sprintf("ctl2nmw: function(s) %s in %s are not in R's deriv() table (cannot auto-differentiate).",
                 paste(bad, collapse = ", "), where), call. = FALSE)
  invisible(TRUE)
}

## Build the FGD/H/PRED source (and the pieces the caller needs) from a parsed CTL.
.nm_build_pred <- function(cs, dose = NULL) {
  nTheta <- cs$n_theta; nEta <- cs$n_eta; nEps <- cs$n_eps; method <- cs$method

  ## parse $PRED assignments into an ordered name -> language map
  defs <- list(); order <- character(0)
  for (ln in cs$pred_lines) {
    if (!grepl("=", ln)) next
    lhs <- trimws(sub("=.*$", "", ln))
    if (!grepl("^[A-Za-z][A-Za-z0-9_]*$", lhs)) next
    rhs <- .nm_translate(sub("^[^=]*=", "", ln))
    e <- tryCatch(parse(text = rhs)[[1]], error = function(z) NULL)
    if (is.null(e)) stop(sprintf("ctl2nmw: cannot parse $PRED line: %s", trimws(ln)), call. = FALSE)
    if (lhs %in% names(defs))
      stop(sprintf(paste0("ctl2nmw: '%s' is assigned more than once in $PRED (line: %s); ",
                          "sequential reassignment is out of scope in v1 -- ",
                          "write each variable as a single expression."),
                   lhs, trimws(ln)), call. = FALSE)
    defs[[lhs]] <- e; order <- unique(c(order, lhs))
  }
  if (length(defs) == 0) stop("ctl2nmw: no assignments found in $PRED.", call. = FALSE)

  ## definitions are inlined substitutionally, so every symbol must be fully
  ## defined before it is used; a self or forward reference would make the
  ## repeated substitute() in .nm_inline() non-terminating
  nm_defs <- names(defs)
  for (k in seq_along(defs)) {
    fwd <- intersect(all.vars(defs[[k]]), nm_defs[seq(k, length(nm_defs))])
    if (length(fwd) > 0)
      stop(sprintf(paste0("ctl2nmw: '%s' references '%s' before it is defined in $PRED; ",
                          "reorder the assignments (each symbol must be defined before use)."),
                   nm_defs[k], fwd[1]), call. = FALSE)
  }

  has_eps <- function(x) any(grepl("^EPS[0-9]+$", all.vars(x)))
  eps_names <- names(defs)[vapply(defs, has_eps, logical(1))]
  if (length(eps_names) == 0)
    stop("ctl2nmw: no EPS/ERR term found in $PRED (need an error model line, e.g. Y = F + F*EPS(1) + EPS(2)).",
         call. = FALSE)
  Yname <- eps_names[length(eps_names)]

  ## literal-constant defs (e.g. DOSE = 320)
  is_lit <- function(x) length(all.vars(x)) == 0 &&
    is.numeric(tryCatch(eval(x), error = function(z) NA))
  lit_names <- names(defs)[vapply(defs, is_lit, logical(1))]
  lit_vals  <- vapply(defs[lit_names], function(x) as.numeric(eval(x)), 0)
  names(lit_vals) <- lit_names
  if (!is.null(dose) && "DOSE" %in% lit_names) lit_vals["DOSE"] <- dose

  ## structural prediction symbol (F, else last non-eps assignment)
  nonEps <- setdiff(order, eps_names)
  Pname <- if ("F" %in% names(defs)) "F"
           else if (length(nonEps)) nonEps[length(nonEps)]
           else stop("ctl2nmw: could not identify the structural prediction (define F in $PRED).",
                     call. = FALSE)

  intermediates <- setdiff(names(defs), c(eps_names, Pname, lit_names))
  inl <- defs[intermediates]

  F_flat <- .nm_inline(defs[[Pname]], inl)
  Y_H    <- .nm_inline(defs[[Yname]], inl)
  if (Pname != "F") Y_H <- do.call("substitute", list(Y_H, setNames(list(quote(F)), Pname)))
  if (!("F" %in% all.vars(Y_H)))
    stop("ctl2nmw: the error model does not reference the structural prediction; cannot build H.",
         call. = FALSE)

  .nm_check_funcs(F_flat, "$PRED F")
  .nm_check_funcs(Y_H, "$PRED error model")

  allsym <- union(all.vars(F_flat), all.vars(Y_H))
  ## cross-checks
  if (length(.nm_idx(allsym, "ETA")) && max(.nm_idx(allsym, "ETA")) > nEta)
    stop(sprintf("ctl2nmw: $PRED uses ETA(%d) but $OMEGA has only %d ETA(s).",
                 max(.nm_idx(allsym, "ETA")), nEta), call. = FALSE)
  if (length(.nm_idx(allsym, "EPS")) && max(.nm_idx(allsym, "EPS")) > nEps)
    stop(sprintf("ctl2nmw: error model uses EPS(%d) but $SIGMA has only %d EPS(s).",
                 max(.nm_idx(allsym, "EPS")), nEps), call. = FALSE)
  if (length(.nm_idx(allsym, "TH")) && max(.nm_idx(allsym, "TH")) > nTheta)
    stop(sprintf("ctl2nmw: $PRED uses THETA(%d) but $THETA declares only %d.",
                 max(.nm_idx(allsym, "TH")), nTheta), call. = FALSE)

  data_items <- setdiff(allsym,
    c(grep("^(TH|ETA|EPS)[0-9]+$", allsym, value = TRUE), lit_names, "F"))

  ## paste helpers that stay length-0 when the input is empty
  ## (base paste0() drops zero-length args, e.g. paste0("TH", integer(0)) == "TH")
  P  <- function(pre, x, post = "") if (length(x)) paste0(pre, x, post) else character(0)
  Plit <- function(nm, v) if (length(nm)) paste0(nm, " = ", v[nm]) else character(0)

  ## ---- FGD (structural derivatives) ----
  Fsym <- all.vars(F_flat)
  f_th   <- .nm_idx(Fsym, "TH")
  f_lit  <- intersect(lit_names, Fsym)
  f_data <- intersect(data_items, Fsym)
  eta_nm <- paste0("ETA", seq_len(nEta))
  fgd_args <- c(P("TH", f_th), eta_nm, f_lit, f_data)
  Ftxt <- paste(deparse(F_flat, width.cutoff = 500L), collapse = " ")
  fgd_src <- sprintf(
    "FGD = deriv(~%s,\n            c(%s),\n            function.arg = c(%s),\n            func = TRUE, hessian = %s)",
    Ftxt, paste0('"', eta_nm, '"', collapse = ", "),
    paste0('"', fgd_args, '"', collapse = ", "),
    if (method == "LAPL") "TRUE" else "FALSE")

  ## ---- H (residual error derivatives), F kept symbolic ----
  Ysym <- all.vars(Y_H)
  h_th   <- .nm_idx(Ysym, "TH")
  h_data <- intersect(data_items, Ysym)
  eps_nm <- paste0("EPS", seq_len(nEps))
  h_args <- c("F", P("TH", h_th), h_data, eps_nm)
  Ytxt <- paste(deparse(Y_H, width.cutoff = 500L), collapse = " ")
  h_src <- sprintf(
    "H = deriv(~%s,\n          c(%s),\n          function.arg = c(%s), func = TRUE)",
    Ytxt, paste0('"', eps_nm, '"', collapse = ", "),
    paste0('"', h_args, '"', collapse = ", "))

  ## ---- PRED wiring ----
  fgd_call <- c(P("THETA[", f_th, "]"),
                paste0("ETA[", seq_len(nEta), "]"),
                Plit(f_lit, lit_vals),
                P('DATAi[, "', f_data, '"]'))
  h_call <- c("FGDres", P("THETA[", h_th, "]"),
              P('DATAi[, "', h_data, '"]'), rep("0", nEps))
  Gnm <- paste0("G", seq_len(nEta)); Hnm <- paste0("H", seq_len(nEps))
  Dnm <- character(0); Dextract <- character(0)
  for (i in seq_len(nEta)) for (j in seq_len(i)) {
    Dnm <- c(Dnm, paste0("D", i, j)); Dextract <- c(Dextract, sprintf("Dres[, %d, %d]", i, j))
  }
  col_base <- c("F", Gnm, Hnm)
  pred_src <- paste0(
    "PRED = function(THETA, ETA, DATAi)\n{\n",
    "  FGDres = FGD(", paste(fgd_call, collapse = ", "), ")\n",
    "  Gres = attr(FGDres, \"gradient\")\n",
    "  Hres = attr(H(", paste(h_call, collapse = ", "), "), \"gradient\")\n",
    "  if (e$METHOD == \"LAPL\") {\n",
    "    Dres = attr(FGDres, \"hessian\")\n",
    "    Res = cbind(FGDres, Gres, Hres, ", paste(Dextract, collapse = ", "), ")\n",
    "    colnames(Res) = c(", paste0('"', c(col_base, Dnm), '"', collapse = ", "), ")\n",
    "  } else {\n",
    "    Res = cbind(FGDres, Gres, Hres)\n",
    "    colnames(Res) = c(", paste0('"', col_base, '"', collapse = ", "), ")\n",
    "  }\n",
    "  return(Res)\n}")

  model_src <- paste(fgd_src, h_src, pred_src, sep = "\n\n")
  list(model_src = model_src, data_items = data_items,
       lit_vals = lit_vals, f_data = f_data, h_data = h_data)
}

## ------------------------------------------------------------- ctl2nmw

ctl2nmw <- function(ctl, eval = FALSE, dose = NULL, file = NULL) {
  cs <- ParseCtl(ctl)
  bp <- .nm_build_pred(cs, dose = dose)

  ## init block (runnable + documented)
  vecstr <- function(v) paste0("c(", paste(format(v, trim = TRUE), collapse = ", "), ")")
  matstr <- function(M) sprintf("matrix(c(%s), nrow = %d)",
                                paste(format(as.vector(M), trim = TRUE), collapse = ", "), nrow(M))
  used <- cs$input$used
  needcols <- unique(c("ID", "TIME", "DV", bp$f_data, bp$h_data))
  init_src <- paste0(
    "## ---- data (you must supply DataAll) ----\n",
    "# DataAll: a data.frame whose columns are named per $INPUT.\n",
    "# $INPUT columns (DROP removed): ", paste(used, collapse = ", "), "\n",
    "# nmw addresses columns BY NAME (not $INPUT order); required here: ",
    paste(needcols, collapse = ", "), "\n\n",
    "## ---- initial estimates (from CTL) ----\n",
    "THETAinit = ", vecstr(cs$theta_init), "\n",
    "LB = ", vecstr(cs$lb), "\n",
    "UB = ", vecstr(cs$ub), "\n",
    "OMinit = ", matstr(cs$om_init), "\n",
    "SGinit = ", matstr(cs$sg_init), "\n",
    "METHOD = \"", cs$method, "\"\n")
  run_src <- paste0(
    "## ---- run ----\n",
    "# InitStep(DataAll, THETAinit, OMinit, SGinit, LB, UB, Pred = PRED, METHOD = METHOD)\n",
    "# EstStep(); CovStep()",
    if (cs$posthoc) "; PostHocEta()" else "", "; TabStep()")
  note_src <- if (length(cs$notes))
    paste0("## ---- notes ----\n", paste0("# - ", cs$notes, collapse = "\n"), "\n\n") else ""

  header <- paste0("## Auto-generated by nmw::ctl2nmw()",
                   if (nzchar(cs$prob_title)) paste0("  [", cs$prob_title, "]") else "", "\n\n")
  code <- paste0(header, note_src, init_src, "\n",
                 "## ---- model (auto-generated PRED) ----\n", bp$model_src, "\n\n", run_src, "\n")

  if (!is.null(file)) writeLines(code, file)

  PRED <- NULL
  if (isTRUE(eval)) {
    env <- new.env(parent = topenv(environment()))
    eval(parse(text = bp$model_src), envir = env)
    PRED <- get("PRED", envir = env)
  }

  structure(list(
    code = code, PRED = PRED,
    THETAinit = cs$theta_init, LB = cs$lb, UB = cs$ub,
    OMinit = cs$om_init, SGinit = cs$sg_init, METHOD = cs$method,
    nTheta = cs$n_theta, nEta = cs$n_eta, nEps = cs$n_eps,
    input = cs$input, notes = cs$notes, parse = cs
  ), class = "nmwCode")
}

print.nmwCode <- function(x, ...) {
  cat(x$code)
  invisible(x)
}
