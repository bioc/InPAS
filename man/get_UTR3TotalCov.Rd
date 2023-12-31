% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/10.get_UTR3TotalCov.R
\name{get_UTR3TotalCov}
\alias{get_UTR3TotalCov}
\title{extract coverage of 3' UTR for CP sites prediction}
\usage{
get_UTR3TotalCov(
  chr.utr3,
  chr.totalCov,
  gcCompensation = NA,
  mappabilityCompensation = NA,
  FFT = FALSE,
  fft.sm.power = 20
)
}
\arguments{
\item{chr.utr3}{An object of \link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges}. It must be an
element of the output of \code{\link[=extract_UTR3Anno]{extract_UTR3Anno()}} for a given chromosome.}

\item{chr.totalCov}{total coverage for each condition of a given chromosome. It
must be an output of \code{\link[=get_totalCov]{get_totalCov()}}}

\item{mappabilityCompensation}{mappability compensation vector. Not support
yet.}

\item{FFT}{Use FFT smooth or not.}

\item{fft.sm.power}{the cut-off frequency of FFT smooth.}

\item{gcCompensationensation}{GC compensation vector. Not support yet.}
}
\value{
path to a file storing the UTR3 total coverage for a given
chromosome/scaffold
}
\description{
extract 3' UTR coverage from totalCov according to the
\link[GenomicRanges:GRanges-class]{GenomicRanges::GRanges} object utr3.
}
\author{
Jianhong Ou
}
\keyword{internal}
