\name{NEWS}
\title{News for Package \pkg{nmw}}
\section{Version 0.3.1 (2026-06-09)}{
  \itemize{
    \item{\code{SumOut()} bug fix: the line \code{RunFolders <- MatchEnd(Folders, RunExt)} previously assigned the logical vector returned by \code{MatchEnd()} to \code{RunFolders} instead of using it as a mask to subset \code{Folders}.  This caused \code{strsplit()} inside \code{GetModelNames()} to fail with "non-character argument" whenever any run folders were present.  Now \code{RunFolders <- Folders[MatchEnd(Folders, RunExt)]}.}
    \item{\code{SumOut()} default \code{RunExt} changed from \code{".R75"} to \code{".R76"} to match NONMEM 7.6 (the current production release).  Users on NONMEM 7.5 should pass \code{RunExt = ".R75"} explicitly.}
    \item{\code{SumOut()} now warns and returns \code{invisible(NULL)} when no matching run folders are found, instead of erroring on a zero-length loop.}
    \item{\code{GetProbVal()} bug fix: \code{grep(Tag, PROB)} returns \code{integer(0)} when the requested tag (e.g. \code{"P:"} or \code{"F:"}) is absent — for example when the \code{$PROB} record is longer than 80 characters and NM-TRAN truncates it inside the \code{FCON} file, dropping the tag.  The previous \code{if (Loc > 0)} test raised "argument is of length zero".  The test was changed to \code{if (length(Loc) > 0 && Loc[1] > 0)}, returning an empty string when the tag is absent.}
    \item{\code{man/SumOut.Rd} updated to reflect the new default \code{RunExt} and a richer description of the return value.}
    \item{\code{MatchEnd()} bug fix: replaced \code{for (i in 1:n)} with \code{for (i in seq_len(n))} so that calling on an empty character vector no longer iterates over \code{c(1, 0)} and tries to read \code{NA} elements.}
    \item{Regression tests added in \code{tests/testthat/test-utils_data.R} covering the \code{SumOut()} empty-folder, \code{GetProbVal()} missing-tag, and \code{MatchEnd()} empty-input cases.}
    \item{\code{report_ebe()} (S5 \dQuote{ETA vs Categorical Variables}) bug fix: when the number of ETAs times the number of categorical covariates produced a panel grid too dense for boxplot margins, \code{plot.new()} raised \dQuote{figure margins too large} and the PDF was truncated mid-page.  The categorical-covariate panels are now split across multiple PDF pages (cap of ~35 panels per page, configurable through the internal \code{MAX_CELLS_PER_PAGE} constant), with page numbers appended to the panel title when more than one page is produced.}
  }
}

\section{Version 0.3.0 (2026-05-08)}{
  \itemize{
    \item{General-purpose helpers for NONMEM dataset construction added (\code{build_nm_dataset}, \code{build_dose_records}, \code{build_obs_records}).}
    \item{Time-varying covariate merge by LOCF/LOCB: \code{locf_value}, \code{merge_cov_locf}.}
    \item{Datetime parsers: \code{parse_dtc}, \code{dat2_time_to_posix}.}
    \item{Coding helpers: \code{code_sex}, \code{code_race}.}
    \item{Cockcroft-Gault CRCL: \code{crcl_cg}, \code{add_crcl_cg}.}
    \item{Mass to molar unit conversion: \code{nm_to_molar}.}
  }
}

\section{Version 0.2.0 (2026-03-08)}{
  \itemize{
    \item{Post-run processing of NONMEM output with PDF diagnostic reports merged from nma package.}
    \item{8 report types: OFV, Parameters, Predictions, Residuals, EBE, Input, Output, IndiPK.}
    \item{New utility functions for NONMEM data processing and parsing.}
    \item{Switched from exportPattern to explicit exports in NAMESPACE.}
  }
}

\section{Version 0.1.6 (2026-03-08)}{
  \itemize{
    \item{Typos in documentation corrected.}
  }
}

\section{Version 0.1.5 (2023-05-09)}{
  \itemize{
    \item{AddCox function for adding a covariate to existing NONMEM dataset.}
    \item{CombDmExPc function for combining demographics(DM), dosing(EX), and DV(PC) tables.}
    \item{TrimOut function to beautify the original NONMEM OUTPUT}
    \item{Errata in manual corrected.}
  }
}

\section{Version 0.1.4 (2018-04-10)}{
  \itemize{
    \item{Errata in manual corrected.}
  }
}

\section{Version 0.1.3 (2017-10-23)}{
  \itemize{
    \item{\code{TabStep()} is added.}
    \item{Speed is increased.}
    \item{References are updated.}
    \item{Vignette is updated.}
  }
}

\section{Version 0.1.2 (2017-10-20)}{
  \itemize{
    \item{This package will be used at the workshop on 2017-10-25 at Kyoto Japan 1st APC satellite workshop.}
    \item{Bug at calculating objective function value for "LAPL" method is fixed.}
    \item{Bug at scaling other than zero lower bound transformation fixed.}
  }
}

\section{Version 0.1.1 (2017-03-14)}{
  \itemize{
    \item{This package will be used at the workshop within US FDA on 2017-03-21.}
    \item{For the initialization step, \code{Pred} function name should be given by the user. It does not read directly from the file any more.}
  }
}

\section{Version 0.1.0 (2017-02-21)}{
  \itemize{
    \item{This package was born and used at the Workshop held by Kyun-Seop Bae MD PhD <k@acr.kr>.}
    \item{\code{InitStep} used \code{Pred} function within a file.}
  }
}
