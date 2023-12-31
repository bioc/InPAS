% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/15.find_valleyBySpline.R
\name{find_valleyBySpline}
\alias{find_valleyBySpline}
\title{Find major valleys after spline smoothing}
\usage{
find_valleyBySpline(
  x,
  ss,
  se = length(x),
  nknots = ceiling((se - ss + 1)/1000 * 10),
  n = -1,
  min.dist = 200,
  filter.last = TRUE,
  DIST2END = 1200,
  plot = FALSE
)
}
\arguments{
\item{x}{A vector of numeric(n), containing MSEs for a given range}

\item{ss}{An positive integer, search start site relative to the leftmost base}

\item{se}{An positive integer, search end site relative to the leftmost base}

\item{nknots}{An positive integer, the number of knots for smoothing using
spline\code{\link[stats:smooth.spline]{stats::smooth.spline()}}. By default, set to 10 knots per kb.}

\item{n}{An integer, specifying the number of location where MSE are local
minima (candidate CP sites). If set to -1, return all candidate CP sites.}

\item{min.dist}{An integer, minimal distance allowed between two adjacent
candidate CP sites otherwise collapsed by selecting the one with lower MSE.}

\item{filter.last}{A logical(1), whether to filter out the last valley, which
is likely the 3' end of the longer 3' UTR if no novel distal CP site is
detected and the 3' end excluded by setting cutEnd/search_point_END is small.}

\item{DIST2END}{An integer, specifying a cutoff of the distance between last valley
and the end of the 3' UTR (where MSE of the last base is calculated). If
the last valley is closer to the end than the specified distance, it will be
not be considered because it is very likely due to RNA coverage decay at the
end of mRNA. Default is 1200. User can consider a value between 1000 and
1500, depending on the library preparation procedures: RNA fragmentation and
size selection.}

\item{plot}{A logical(1), whether to plot the MSE profile and the candidate
valleys.}
}
\value{
A vector of integer.
}
\description{
Find major valleys after spline smoothing
}
\keyword{internal}
