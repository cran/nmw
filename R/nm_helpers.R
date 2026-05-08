# =============================================================================
# nm_helpers.R  --  General-purpose helpers for NONMEM dataset construction
#
# Designed as composable building blocks.  Base R only.
#
# Pipeline:
#   build_nm_dataset(DM, EX, PC, ...)   --> NM data.frame (mass units)
#     |--> merge_cov_locf(...)          --  add time-varying covariates
#     |--> add_crcl_cg(...)             --  compute CRCL
#   nm_to_molar(df, mw_dose, mw_dv)     --> NM data.frame (molar units)
#
# Saving to CSV is the caller's responsibility.
# =============================================================================


# ---- 1. coding helpers ------------------------------------------------------

#' Code SEX as Integer
#'
#' Code SEX as integer: M -> 0, F -> 1, otherwise -> 2.
#' Numeric input is returned unchanged (cast to integer).
#'
#' @param x character or numeric vector of sex values.
#' @return integer vector.
#' @export
code_sex <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  s   <- toupper(trimws(as.character(x)))
  out <- rep(2L, length(s))
  out[s == "M" | s == "MALE"]   <- 0L
  out[s == "F" | s == "FEMALE"] <- 1L
  out
}

#' Code RACE as Integer
#'
#' Code RACE as integer: ASIAN -> 1, WHITE -> 2, BLACK -> 3, otherwise -> 4.
#' Numeric input is returned unchanged (cast to integer).
#'
#' @param x character or numeric vector of race values.
#' @return integer vector.
#' @export
code_race <- function(x) {
  if (is.numeric(x)) return(as.integer(x))
  s   <- toupper(trimws(as.character(x)))
  out <- rep(4L, length(s))
  out[substr(s, 1, 5) == "ASIAN"] <- 1L
  out[substr(s, 1, 5) == "WHITE"] <- 2L
  out[substr(s, 1, 5) == "BLACK"] <- 3L
  out
}

#' Cockcroft-Gault Creatinine Clearance
#'
#' Cockcroft-Gault creatinine clearance (mL/min).
#'
#' @param age numeric, age in years.
#' @param bwt numeric, body weight in kg.
#' @param sex integer/numeric, 0 = male, 1 = female.
#' @param crea numeric, serum creatinine in mg/dL.
#' @return numeric vector of CRCL in mL/min.
#' @export
crcl_cg <- function(age, bwt, sex, crea) {
  (140 - age) * bwt * (1 - 0.15 * sex) / (72 * crea)
}

#' Add CRCL Column via Cockcroft-Gault
#'
#' Add a CRCL column to a NM-style data.frame using Cockcroft-Gault.
#' Column names are configurable.
#'
#' @param df data.frame containing AGE, BWT, SEX, CREA columns.
#' @param age_col name of the age column (years). Default "AGE".
#' @param bwt_col name of the body-weight column (kg). Default "BWT".
#' @param sex_col name of the sex column (0=M, 1=F). Default "SEX".
#' @param crea_col name of the creatinine column (mg/dL). Default "CREA".
#' @param out_col name of the output column. Default "CRCL".
#' @return data.frame with the CRCL column added.
#' @export
add_crcl_cg <- function(df,
                        age_col  = "AGE",
                        bwt_col  = "BWT",
                        sex_col  = "SEX",
                        crea_col = "CREA",
                        out_col  = "CRCL") {
  df[[out_col]] <- crcl_cg(df[[age_col]], df[[bwt_col]],
                           df[[sex_col]], df[[crea_col]])
  df
}


# ---- 2. datetime helpers ----------------------------------------------------

#' Combine Date and Time Strings into POSIXct
#'
#' Combine "YYYY-MM-DD" and "HH:MM" character vectors into POSIXct.
#'
#' @param dat2 character vector of dates in "YYYY-MM-DD" format.
#' @param time character vector of times in "HH:MM" format.
#' @param tz time zone. Default "UTC".
#' @return POSIXct vector.
#' @export
dat2_time_to_posix <- function(dat2, time, tz = "UTC") {
  s <- paste(as.character(dat2), as.character(time))
  as.POSIXct(s, format = "%Y-%m-%d %H:%M", tz = tz)
}

#' Parse SDTM-Style DTC Strings into POSIXct
#'
#' Parse SDTM-style datetime strings into POSIXct.
#' Accepts "YYYY-MM-DDTHH:MM:SS", "YYYY-MM-DDTHH:MM",
#' "YYYY-MM-DD HH:MM:SS", "YYYY-MM-DD HH:MM", or "YYYY-MM-DD".
#' POSIXct / Date inputs are returned coerced to POSIXct.
#'
#' @param x character vector (or POSIXct/Date) of datetime values.
#' @param tz time zone. Default "UTC".
#' @return POSIXct vector.
#' @export
parse_dtc <- function(x, tz = "UTC") {
  if (inherits(x, "POSIXct")) return(x)
  if (inherits(x, "Date"))    return(as.POSIXct(x, tz = tz))
  s   <- as.character(x)
  s   <- gsub("T", " ", s)
  out <- rep(as.POSIXct(NA, tz = tz), length(s))
  ix  <- !is.na(s) & nchar(s) > 0
  if (!any(ix)) return(out)

  has_hms <- grepl("\\d{2}:\\d{2}:\\d{2}", s)
  has_hm  <- grepl("\\d{2}:\\d{2}",        s) & !has_hms

  out[has_hms] <- as.POSIXct(s[has_hms], format = "%Y-%m-%d %H:%M:%S", tz = tz)
  out[has_hm]  <- as.POSIXct(s[has_hm],  format = "%Y-%m-%d %H:%M",    tz = tz)
  rest <- ix & is.na(out)
  if (any(rest))
    out[rest] <- as.POSIXct(s[rest], format = "%Y-%m-%d", tz = tz)
  out
}


# ---- 3. LOCF / LOCB ---------------------------------------------------------

#' LOCF / LOCB Lookup
#'
#' For each target time, return the most-recent source value at
#' source_t <= target_t (LOCF). If no preceding source exists, return the
#' earliest source_v at source_t > target_t (LOCB). Empty source returns
#' \code{fallback} for all targets.
#'
#' @param target_t vector of target times (numeric or POSIXct).
#' @param source_t vector of source times (same type as target_t).
#' @param source_v numeric vector of source values, parallel to source_t.
#' @param fallback value used when no source is available. Default NA.
#' @return numeric vector, length(target_t).
#' @export
locf_value <- function(target_t, source_t, source_v, fallback = NA_real_) {
  n <- length(target_t)
  if (length(source_t) == 0) return(rep(fallback, n))
  ord <- order(source_t)
  st  <- source_t[ord]
  sv  <- source_v[ord]
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    tt <- target_t[i]
    if (is.na(tt)) next
    le <- which(st <= tt)
    if (length(le) > 0) {
      out[i] <- sv[max(le)]
    } else {
      ge <- which(st > tt)
      if (length(ge) > 0) out[i] <- sv[min(ge)]
    }
  }
  out[is.na(out)] <- fallback
  out
}

#' Merge Time-Varying Covariates by LOCF
#'
#' Merge time-varying covariate columns into a NM-style data.frame using LOCF
#' (with LOCB for records before the first source). Subjects/records still NA
#' after LOCF/LOCB can be filled with the column-wise median of \code{cov_df}.
#'
#' @param nm data.frame with key columns \code{subj_col}, \code{dat_col},
#'   \code{tim_col}.
#' @param cov_df data.frame with \code{subj_col}, \code{time_col}, and one or
#'   more value columns.
#' @param value_cols character vector of value columns in \code{cov_df} to
#'   bring into \code{nm}.
#' @param time_col name of the time column in \code{cov_df}, parsed by
#'   \code{\link{parse_dtc}}.
#' @param median_fb logical. If TRUE (default), fill remaining NA values with
#'   the column-wise median of \code{cov_df}.
#' @param verbose logical. Print median fill counts. Default FALSE.
#' @param subj_col subject ID column name. Default "SUBJID".
#' @param dat_col date column name in \code{nm}. Default "DAT2".
#' @param tim_col time column name in \code{nm}. Default "TIME".
#' @return \code{nm} with \code{value_cols} added.
#' @export
merge_cov_locf <- function(nm, cov_df, value_cols, time_col,
                           median_fb = TRUE, verbose = FALSE,
                           subj_col = "SUBJID",
                           dat_col  = "DAT2",
                           tim_col  = "TIME") {
  if (is.null(cov_df) || nrow(cov_df) == 0) {
    for (vc in value_cols) nm[[vc]] <- NA_real_
    return(nm)
  }
  nm_pt  <- dat2_time_to_posix(nm[[dat_col]], nm[[tim_col]])
  cov_pt <- parse_dtc(cov_df[[time_col]])
  for (vc in value_cols) nm[[vc]] <- NA_real_

  subjs <- unique(nm[[subj_col]])
  for (sid in subjs) {
    nm_idx  <- which(nm[[subj_col]] == sid)
    cov_idx <- which(cov_df[[subj_col]] == sid & !is.na(cov_pt))
    for (vc in value_cols) {
      vv   <- suppressWarnings(as.numeric(cov_df[[vc]][cov_idx]))
      keep <- !is.na(vv)
      nm[[vc]][nm_idx] <- locf_value(
        target_t = nm_pt[nm_idx],
        source_t = cov_pt[cov_idx][keep],
        source_v = vv[keep],
        fallback = NA_real_
      )
    }
  }

  if (median_fb) {
    for (vc in value_cols) {
      if (any(is.na(nm[[vc]]))) {
        med <- median(suppressWarnings(as.numeric(cov_df[[vc]])), na.rm = TRUE)
        if (verbose) {
          n_na <- sum(is.na(nm[[vc]]))
          cat(sprintf("  [merge_cov_locf] %s: %d records filled with median = %g\n",
                      vc, n_na, med))
        }
        nm[[vc]][is.na(nm[[vc]])] <- med
      }
    }
  }
  nm
}


# ---- 4. record builders -----------------------------------------------------

#' Build NONMEM Dose Records
#'
#' Build NONMEM dose records from an EX-like data.frame.
#' EX must contain: SUBJID, DAT2, TIME, AMT, RATE.
#' Optional: CMT (defaults to \code{dose_cmt}).
#' Output rows have MDV = 1, DV = 0.
#'
#' @param EX data.frame with dose information.
#' @param dose_cmt default compartment for dose records lacking a CMT column.
#' @return data.frame of dose records with columns
#'   SUBJID, DAT2, TIME, AMT, RATE, CMT, DV, MDV.
#' @export
build_dose_records <- function(EX, dose_cmt = 1L) {
  cmt <- if ("CMT" %in% names(EX)) as.integer(EX$CMT) else rep(as.integer(dose_cmt), nrow(EX))
  data.frame(
    SUBJID = as.character(EX$SUBJID),
    DAT2   = as.character(EX$DAT2),
    TIME   = as.character(EX$TIME),
    AMT    = as.numeric(EX$AMT),
    RATE   = as.numeric(EX$RATE),
    CMT    = cmt,
    DV     = 0,
    MDV    = 1L,
    stringsAsFactors = FALSE
  )
}

#' Build NONMEM Observation Records
#'
#' Build NONMEM observation records from a PC-like data.frame.
#' PC must contain: SUBJID, DAT2, TIME, DV, CMT.
#' Output rows have AMT = 0, RATE = 0; DV that is NA or 0 -> MDV = 1.
#'
#' @param PC data.frame with observation information.
#' @return data.frame of observation records with columns
#'   SUBJID, DAT2, TIME, AMT, RATE, CMT, DV, MDV.
#' @export
build_obs_records <- function(PC) {
  dv  <- suppressWarnings(as.numeric(PC$DV))
  mdv <- ifelse(is.na(dv) | dv == 0, 1L, 0L)
  dv[is.na(dv)] <- 0
  data.frame(
    SUBJID = as.character(PC$SUBJID),
    DAT2   = as.character(PC$DAT2),
    TIME   = as.character(PC$TIME),
    AMT    = 0,
    RATE   = 0,
    CMT    = as.integer(PC$CMT),
    DV     = dv,
    MDV    = mdv,
    stringsAsFactors = FALSE
  )
}


# ---- 5. main NM dataset builder --------------------------------------------

#' Build a NONMEM-Format Dataset from DM/EX/PC
#'
#' Build a NONMEM-format dataset from cleaned DM/EX/PC tables.
#' Time-varying covariates (VS, LB, ...) are NOT merged here.  Use
#' \code{\link{merge_cov_locf}} afterwards to add them.  CRCL is added
#' separately via \code{\link{add_crcl_cg}}.  Molar conversion is done with
#' \code{\link{nm_to_molar}}.
#'
#' Records are sorted by SUBJID, DAT2, TIME, CMT, MDV, AMT.  At tied
#' (SUBJID, DAT2, TIME, CMT, MDV), AMT = 0 (observation) sorts before AMT > 0
#' (dose) so a pre-dose observation precedes the dose given at the same minute.
#'
#' @param DM data.frame with one row per subject.  Must contain SUBJID; any
#'   other columns are merged in as subject-level constants.
#' @param EX dose records (see \code{\link{build_dose_records}}).
#' @param PC observation records (see \code{\link{build_obs_records}}).
#' @param IDs character vector of SUBJIDs to keep.  NULL (default) is the
#'   intersect of EX and PC where DV > 0.
#' @param id_prefix character.  When \code{id_func} is NULL,
#'   ID = paste0(id_prefix, SUBJID).
#' @param id_func function(SUBJID) -> character ID.  Overrides \code{id_prefix}.
#' @param dose_cmt default compartment for dose records lacking a CMT column.
#' @param verbose print progress. Default FALSE.
#' @return data.frame with columns
#'   ID, SUBJID, DAT2, TIME, AMT, RATE, CMT, DV, MDV, <DM-columns ...>.
#' @export
build_nm_dataset <- function(DM, EX, PC,
                             IDs       = NULL,
                             id_prefix = "",
                             id_func   = NULL,
                             dose_cmt  = 1L,
                             verbose   = FALSE) {

  # ---- 5a. validate ----
  req <- list(
    DM = "SUBJID",
    EX = c("SUBJID", "DAT2", "TIME", "AMT", "RATE"),
    PC = c("SUBJID", "DAT2", "TIME", "DV", "CMT")
  )
  dfs <- list(DM = DM, EX = EX, PC = PC)
  for (nm in names(req)) {
    miss <- setdiff(req[[nm]], names(dfs[[nm]]))
    if (length(miss) > 0)
      stop(sprintf("%s: missing column(s): %s", nm, paste(miss, collapse = ", ")))
  }

  # ---- 5b. resolve subject IDs ----
  if (is.null(IDs)) {
    ids_ex <- unique(as.character(EX$SUBJID))
    pc_dv  <- suppressWarnings(as.numeric(PC$DV))
    ids_pc <- unique(as.character(PC$SUBJID[!is.na(pc_dv) & pc_dv > 0]))
    IDs    <- intersect(ids_ex, ids_pc)
  }
  if (length(IDs) == 0) stop("No subjects to include (empty IDs).")
  if (verbose) cat(sprintf("[build_nm_dataset] %d subjects\n", length(IDs)))

  DM <- DM[as.character(DM$SUBJID) %in% IDs, , drop = FALSE]
  DM <- DM[!duplicated(DM$SUBJID), , drop = FALSE]
  EX <- EX[as.character(EX$SUBJID) %in% IDs, , drop = FALSE]
  PC <- PC[as.character(PC$SUBJID) %in% IDs, , drop = FALSE]

  # ---- 5c. dose + observation records ----
  EXr <- build_dose_records(EX, dose_cmt = dose_cmt)
  PCr <- build_obs_records(PC)
  nm  <- rbind(EXr, PCr)
  nm  <- nm[order(nm$SUBJID, nm$DAT2, nm$TIME, nm$CMT, nm$MDV, nm$AMT), ,
            drop = FALSE]
  if (verbose) cat(sprintf("[build_nm_dataset] %d records (%d dose + %d obs)\n",
                            nrow(nm), nrow(EXr), nrow(PCr)))

  # ---- 5d. merge subject-level DM (any columns the caller has) ----
  nm <- merge(nm, DM, by = "SUBJID", all.x = TRUE, sort = FALSE)
  nm <- nm[order(nm$SUBJID, nm$DAT2, nm$TIME, nm$CMT, nm$MDV, nm$AMT), ,
           drop = FALSE]

  # ---- 5e. ID column ----
  if (!is.null(id_func)) {
    nm$ID <- id_func(nm$SUBJID)
  } else {
    nm$ID <- paste0(id_prefix, nm$SUBJID)
  }

  front <- c("ID", "SUBJID", "DAT2", "TIME", "AMT", "RATE", "CMT", "DV", "MDV")
  rest  <- setdiff(names(nm), front)
  nm    <- nm[, c(front, rest), drop = FALSE]
  rownames(nm) <- NULL
  nm
}


# ---- 6. molar conversion ----------------------------------------------------

#' Convert NM Dataset from Mass to Molar Units
#'
#' Convert AMT/RATE/DV columns of a NM data.frame from mass units to molar units.
#'
#' Typical use:
#' \itemize{
#'   \item AMT, RATE in mg, DV in ng/mL, MW in g/mol
#'   \item mg / MW * 1000  ->  umol
#'   \item ng/mL / MW      ->  umol/L
#' }
#'
#' For multiple analytes (e.g., parent at CMT=1 and metabolite at CMT=3),
#' supply \code{mw_dv} as a named numeric vector with names equal to the CMT
#' values (as character):  \code{c(`1` = 394.47, `3` = 380.47)}.
#'
#' @param df data.frame.
#' @param mw_dose numeric scalar.  MW for AMT and RATE columns.
#' @param mw_dv numeric scalar OR named numeric.
#'   Scalar: applied to all DV regardless of CMT.
#'   Named: names are CMT values; rows whose CMT is not in names() are left
#'   unchanged.  Default = mw_dose.
#' @param amt_cols columns to convert with mw_dose.  Default c("AMT", "RATE").
#' @param dv_col DV column.  Default "DV".
#' @param cmt_col CMT column.  Default "CMT".
#' @param dose_factor multiplier applied after AMT/RATE divided by mw_dose.
#'   1000 = mg -> umol (default).
#' @param dv_factor multiplier applied after DV divided by mw_dv.
#'   1 = ng/mL -> umol/L (default).
#' @return data.frame with converted columns; other columns unchanged.
#' @export
nm_to_molar <- function(df, mw_dose, mw_dv = mw_dose,
                        amt_cols    = c("AMT", "RATE"),
                        dv_col      = "DV",
                        cmt_col     = "CMT",
                        dose_factor = 1000,
                        dv_factor   = 1) {
  out <- df

  for (a in amt_cols) {
    if (a %in% names(out))
      out[[a]] <- out[[a]] / mw_dose * dose_factor
  }

  if (length(mw_dv) == 1 && is.null(names(mw_dv))) {
    out[[dv_col]] <- out[[dv_col]] / mw_dv * dv_factor
  } else {
    if (is.null(names(mw_dv)))
      stop("mw_dv must be a named numeric when length > 1.")
    cmt_vals <- as.character(out[[cmt_col]])
    for (k in names(mw_dv)) {
      mask <- cmt_vals == k
      if (any(mask))
        out[[dv_col]][mask] <- out[[dv_col]][mask] / mw_dv[[k]] * dv_factor
    }
  }
  out
}
