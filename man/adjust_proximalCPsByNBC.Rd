% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/15.adjust_proximalCPsByNBC.R
\name{adjust_proximalCPsByNBC}
\alias{adjust_proximalCPsByNBC}
\title{adjust the proximal CP sites by using Naive Bayes classifier from cleanUpdTSeq}
\usage{
adjust_proximalCPsByNBC(
  idx.list,
  cov_diff.list,
  seqnames,
  starts,
  strands,
  genome,
  classifier,
  classifier_cutoff,
  shift_range,
  search_point_START,
  step = 1
)
}
\arguments{
\item{idx.list}{the offset of positions of CP sites}

\item{cov_diff.list}{the MSE values}

\item{seqnames}{a character(n) vector, the chromosome/scaffolds' names}

\item{starts}{starts}

\item{strands}{strands}

\item{genome}{a \link[BSgenome:BSgenome-class]{BSgenome::BSgenome} object}

\item{classifier}{cleanUpdTSeq classifier}

\item{classifier_cutoff}{cutoff value of the classifier}

\item{shift_range}{the searching range for the better CP sites}

\item{search_point_START}{just in case there is no better CP sites}

\item{step}{adjusting step, default 1, means adjust by each base by
cleanUpdTSeq.}
}
\value{
the offset of positions of CP sites after filter
}
\description{
adjust the proximal CP sites by using Naive Bayes classifier from cleanUpdTSeq
}
\details{
the step for calculating is 10, can not do every base base it is
really very slow.
}
\seealso{
\code{\link[=adjust_proximalCPsByPWM]{adjust_proximalCPsByPWM()}}, \code{\link[=get_PAscore2]{get_PAscore2()}}
}
\author{
Jianhong Ou
}
\keyword{internal}
