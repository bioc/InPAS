% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/24.filter_testOut.R
\name{filter_testOut}
\alias{filter_testOut}
\title{filter 3' UTR usage test results}
\usage{
filter_testOut(
  res,
  gp1,
  gp2,
  outdir = getInPASOutputDirectory(),
  background_coverage_threshold = 2,
  P.Value_cutoff = 0.05,
  adj.P.Val_cutoff = 0.05,
  dPDUI_cutoff = 0.2,
  PDUI_logFC_cutoff = log2(1.5)
)
}
\arguments{
\item{res}{a \linkS4class{UTR3eSet} object, output of \code{\link[=test_dPDUI]{test_dPDUI()}}}

\item{gp1}{tag names involved in group 1. gp1 and gp2 are used for filtering
purpose if both are specified; otherwise only other specified thresholds
are used for filtering.}

\item{gp2}{tag names involved in group 2}

\item{outdir}{A character(1) vector, a path with write permission for storing
InPAS analysis results. If it doesn't exist, it will be created.}

\item{background_coverage_threshold}{background coverage cut off value. for
each group, more than half of the long form should greater than
background_coverage_threshold. for both group, at least in one group, more
than half of the short form should greater than
background_coverage_threshold.}

\item{P.Value_cutoff}{cutoff of P value}

\item{adj.P.Val_cutoff}{cutoff of adjust P value}

\item{dPDUI_cutoff}{cutoff of dPDUI}

\item{PDUI_logFC_cutoff}{cutoff of PDUI log2 transformed fold change}
}
\value{
A data frame converted from an object of
\link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges}.
}
\description{
filter results of \code{\link[=test_dPDUI]{test_dPDUI()}}
}
\examples{
library(limma)
path <- system.file("extdata", package = "InPAS")
load(file.path(path, "eset.MAQC.rda"))
tags <- colnames(eset@PDUI)
g <- factor(gsub("\\\\..*$", "", tags))
design <- model.matrix(~ -1 + g)
colnames(design) <- c("Brain", "UHR")
contrast.matrix <- makeContrasts(
  contrasts = "Brain-UHR",
  levels = design
)
res <- test_dPDUI(
  eset = eset,
  method = "limma",
  normalize = "none",
  design = design,
  contrast.matrix = contrast.matrix
)
filter_testOut(res,
  gp1 = c("Brain.auto", "Brain.phiX"),
  gp2 = c("UHR.auto", "UHR.phiX"),
  background_coverage_threshold = 2,
  P.Value_cutoff = 0.05,
  adj.P.Val_cutoff = 0.05,
  dPDUI_cutoff = 0.3,
  PDUI_logFC_cutoff = .59
)
}
\seealso{
\code{\link[=test_dPDUI]{test_dPDUI()}}
}
\author{
Jianhong Ou, Haibo Liu
}
