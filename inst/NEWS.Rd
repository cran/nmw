\name{NEWS}
\title{News for Package \pkg{nmw}}
\section{Version 0.5.1 (2026-08-21)}{
  \itemize{
    \item{\code{nmw_report_resid()} now omits unavailable residual diagnostics
      cleanly when \code{sdtab} lacks \code{WRES}, \code{CWRES}, or
      \code{IWRE}, or when such a column contains no finite values. Binomial
      tests, extreme-value tables, distribution plots, and subject plots are
      generated only for usable residuals; missing panels contain an
      explanatory placeholder rather than failing on non-finite plot limits.}
    \item{Diagnostic report functions now document their invisible summary
      return values. The \code{model} argument documentation for
      \code{nmw_report_output()} and \code{nmw_report_indipk()} now states that
      it is accepted only for a uniform report API because those reports read
      the fixed files \file{PRINT.OUT} and \file{patab}.}
    \item{The minimum \pkg{simPDF} version is now 0.1.1, the first CRAN
      release, and the minimum R version is aligned with \pkg{simPDF} at 3.6.0.}
    \item{Regression tests now exercise PDF generation and the missing
      \code{IWRE} case using the bundled NONMEM example.}
  }
}
\section{Version 0.5.0 (2026-07-11)}{
  \itemize{
    \item{The eight diagnostic report functions (\code{nmw_report_ofv},
      \code{nmw_report_param}, \code{nmw_report_pred}, \code{nmw_report_resid},
      \code{nmw_report_ebe}, \code{nmw_report_input}, \code{nmw_report_output},
      \code{nmw_report_indipk}) are re-implemented on the \pkg{simPDF} report
      engine (now in \code{Imports}).  The statistics/parsing is unchanged; the
      former dead-reckoned \code{pdf()} + \code{PrinTxt}/\code{split.screen}
      layout is replaced by a measured flowing layout, which eliminates the
      text/matrix overlap that appeared in \code{S2-Parameters.PDF} for models
      with many random effects (validated overlap-free across hundreds of real
      runs, up to 14 \code{ETA}s).  Each report gained \code{file} and
      \code{model} arguments (backward compatible).}
    \item{\code{nmw_report_pred} no longer errors on subjects/compartments with
      no observation records (the "zero-length labels" crash); \code{nmw_report_resid}
      tolerates an \code{sdtab} that is missing a residual column (e.g. \code{IWRE}).}
    \item{The low-level layout helpers \code{PrepPDF}, \code{PrinTxt},
      \code{AddPage}, \code{PrinMTxt} and \code{ClosePDF} are removed (superseded
      by \pkg{simPDF}); \code{DxPlotPost} is retained.}
  }
}
\section{Version 0.4.0 (2026-07-08)}{
  \itemize{
    \item{New function \code{ctl2nmw()}: reads a closed-form \code{$PRED} NONMEM
      control stream and auto-generates the \code{PRED} function (with its
      \code{deriv()} companions) and the initial estimates
      (\code{THETAinit}, \code{OMinit}, \code{SGinit}, \code{LB}, \code{UB},
      \code{METHOD}) needed by \code{InitStep()}.  This removes the main
      inconvenience of the package: hand-writing the \code{PRED} function.  The
      generated \code{PRED} reproduces the hand-written example \code{PRED}
      numerically (verified to floating-point tolerance for \code{FO},
      \code{FOCE-INTER} and \code{Laplace} on the bundled THEO control streams).}
    \item{The generator classifies free symbols over both the structural model
      \code{F} and the error model \code{Y}, so a \code{THETA} or data item that
      appears only in the residual error (e.g. a combined
      additive-plus-proportional \code{W} term) is wired into \code{H} correctly.
      Random-effect subscripts are ordered numerically (so \code{ETA10} follows
      \code{ETA9}), and \code{$PRED} \code{ETA}/\code{EPS} indices are
      cross-checked against \code{$OMEGA}/\code{$SIGMA}.}
    \item{Out-of-scope control streams fail with an informative message rather
      than generating incorrect code: ODE models
      (\code{$PK}/\code{$DES}/\code{$SUBROUTINE}/\code{ADVAN}),
      \code{IF}/\code{THEN} branching, MU-referencing, mixtures, unsupported
      estimation methods, and \code{FIX}ed parameters.  Assigning the same
      symbol more than once in \code{$PRED} (sequential reassignment) or
      referencing a symbol before it is defined is likewise rejected, as is a
      \code{$THETA} initial estimate on or outside its bounds
      (\code{InitStep()} scaling needs \code{LB < init < UB}).  A diagonal
      \code{$OMEGA} is promoted to a full block, a non-diagonal \code{$SIGMA}
      is reduced to its diagonal, and an unbounded \code{$THETA} with a
      non-positive initial estimate gets \code{LB = -1e6} instead of the usual
      default \code{0}, each with a note.}
    \item{The environment \code{e} is exported again.  It had been exported up
      to 0.1.6 via \code{exportPattern()} and was inadvertently dropped in
      0.2.0's switch to explicit exports.  Without it, a user-written
      \code{PRED} function referencing \code{e$METHOD} -- the documented
      workflow, and also the code that \code{ctl2nmw()} generates -- stopped
      with \dQuote{object 'e' not found} when \code{EstStep()} ran it, and
      reading results with \code{get("EBE", envir=e)} failed likewise.}
    \item{Internal refactoring (prepared as an unreleased 0.3.2) is included
      in this release.  It changes no other public API and no numerical
      result of the estimation steps
      (\code{FO}, \code{FOCE-INTER}, \code{Laplace}) or of any diagnostic
      report; verified by comparing the estimation output and all eight report
      PDFs against the pre-refactor sources on the bundled example run.}
    \item{De-duplication: the covariate-classification block repeated three
      times inside \code{report_ebe} is now the internal helper
      \code{ClassifyCovariates()}; the shared record-count extraction is
      \code{GetNRecFromXML()}; \code{NMVarStat()} reuses the
      \code{AllNA()}/\code{AllSame()} helpers.}
    \item{Bug fix: the \dQuote{Sigmas} section of \code{nmw_report_param}
      (S2-Parameters) no longer errors for models with two or more
      \code{SIGMA} elements (e.g. a combined additive-plus-proportional error
      model); a test that mixed \code{&&} and \code{&} produced a condition of
      length greater than one.}
    \item{Bug fix: \code{PanelCor()} no longer errors when an explicit
      \code{cex.cor} value is supplied (an undefined \code{cex} was referenced).}
    \item{\code{InitStep()}: the documented \code{LB}/\code{UB} defaults
      (\code{0} and \code{1e6}) are now applied before they are used, so
      omitting these arguments works as documented instead of raising an error.
      Passing \code{LB}/\code{UB} explicitly behaves exactly as before.}
    \item{Loop bounds in the data-handling utilities (\code{AddDoNoTaLD},
      \code{ExpandDoseHist}, \code{RemoveNA}, \code{RmvFixed}, \code{RmvZero},
      \code{RenCol}) and in the shared report helpers now use
      \code{seq_len()}/\code{seq_along()} so they no longer misbehave on
      zero-row input.}
    \item{Housekeeping: removed leftover working files and dead commented code
      from \code{R/}, replaced \code{T}/\code{F} literals with
      \code{TRUE}/\code{FALSE}, corrected user-visible typos in report titles
      and labels (\dQuote{Individiual}, \dQuote{Estitmate}, \dQuote{Recods}),
      and corrected the \code{TrimOut()} return-value documentation.}
  }
}

\section{Version 0.3.1 (2026-06-09)}{
  \itemize{
    \item{\code{SumOut()} bug fix: the line \code{RunFolders <- MatchEnd(Folders, RunExt)} previously assigned the logical vector returned by \code{MatchEnd()} to \code{RunFolders} instead of using it as a mask to subset \code{Folders}.  This caused \code{strsplit()} inside \code{GetModelNames()} to fail with "non-character argument" whenever any run folders were present.  Now \code{RunFolders <- Folders[MatchEnd(Folders, RunExt)]}.}
    \item{\code{SumOut()} default \code{RunExt} changed from \code{".R75"} to \code{".R76"} to match NONMEM 7.6 (the current production release).  Users on NONMEM 7.5 should pass \code{RunExt = ".R75"} explicitly.}
    \item{\code{SumOut()} now warns and returns \code{invisible(NULL)} when no matching run folders are found, instead of erroring on a zero-length loop.}
    \item{\code{GetProbVal()} bug fix: \code{grep(Tag, PROB)} returns \code{integer(0)} when the requested tag (e.g. \code{"P:"} or \code{"F:"}) is absent -- for example when the \code{$PROB} record is longer than 80 characters and NM-TRAN truncates it inside the \code{FCON} file, dropping the tag.  The previous \code{if (Loc > 0)} test raised "argument is of length zero".  The test was changed to \code{if (length(Loc) > 0 && Loc[1] > 0)}, returning an empty string when the tag is absent.}
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
