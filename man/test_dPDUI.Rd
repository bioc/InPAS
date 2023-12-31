% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/22.test_dPDUI.R
\name{test_dPDUI}
\alias{test_dPDUI}
\title{do test for dPDUI}
\usage{
test_dPDUI(
  eset,
  sqlite_db,
  outdir = getInPASOutputDirectory(),
  method = c("limma", "fisher.exact", "singleSample", "singleGroup"),
  normalize = c("none", "quantiles", "quantiles.robust", "mean", "median"),
  design,
  contrast.matrix,
  coef = 1,
  robust = FALSE,
  ...
)
}
\arguments{
\item{eset}{An object of \linkS4class{UTR3eSet}. It is an output of
\code{\link[=get_UTR3eSet]{get_UTR3eSet()}}}

\item{sqlite_db}{A path to the SQLite database for InPAS, i.e. the output of
setup_sqlitedb().}

\item{outdir}{A character(1) vector, a path with write permission for storing
InPAS analysis results. If it doesn't exist, it will be created.}

\item{method}{A character(1), indicating the method for testing dPDUI. It can
be "limma", "fisher.exact", "singleSample", or "singleGroup"}

\item{normalize}{A character(1), indicating the normalization method. It can
be "none", "quantiles", "quantiles.robust", "mean", or "median"}

\item{design}{a design matrix of the experiment, with rows corresponding to
samples and columns to coefficients to be estimated. Defaults to the unit
vector meaning that the samples are treated as replicates. see
\code{\link[stats:model.matrix]{stats::model.matrix()}}. Required for limma-based analysis.}

\item{contrast.matrix}{a numeric matrix with rows corresponding to
coefficients in fit and columns containing contrasts. May be a vector if
there is only one contrast. see \code{\link[limma:makeContrasts]{limma::makeContrasts()}}. Required for
limma-based analysis.}

\item{coef}{column number or column name specifying which coefficient or
contrast of the linear model is of interest. see more \code{\link[limma:toptable]{limma::topTable()}}.
default value: 1}

\item{robust}{A logical(1) vector, indicating whether the estimation of the empirical Bayes prior
parameters should be robustified against outlier sample variances.}

\item{...}{other arguments are passed to lmFit}
}
\value{
An object of \linkS4class{UTR3eSet}, with the last element \code{testRes} containing the test results in a matrix.
}
\description{
do test for dPDUI
}
\details{
if method is "limma", design matrix and contrast is required. if
method is "fisher.exact", gp1 and gp2 is required.
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
  sqlite_db,
  method = "limma",
  normalize = "none",
  design = design,
  contrast.matrix = contrast.matrix
)
}
\seealso{
\code{\link[=run_singleSampleAnalysis]{run_singleSampleAnalysis()}}, \code{\link[=run_singleGroupAnalysis]{run_singleGroupAnalysis()}},
\code{\link[=run_fisherExactTest]{run_fisherExactTest()}}, \code{\link[=run_limmaAnalysis]{run_limmaAnalysis()}}
}
\author{
Jianhong Ou, Haibo Liu
}
