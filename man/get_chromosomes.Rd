% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04.get_chromosomes.R
\name{get_chromosomes}
\alias{get_chromosomes}
\title{Identify chromosomes/scaffolds for CP site discovery}
\usage{
get_chromosomes(utr3, sqlite_db)
}
\arguments{
\item{utr3}{An object of \link[GenomicRanges:GRangesList-class]{GenomicRanges::GRangesList}. An output of
\code{\link[=extract_UTR3Anno]{extract_UTR3Anno()}}.}

\item{sqlite_db}{A path to the SQLite database for InPAS, i.e. the output of
\code{\link[=setup_sqlitedb]{setup_sqlitedb()}}.}
}
\value{
A vector of characters, containing names of chromosomes/scaffolds
for CP site discovery
}
\description{
Identify chromosomes/scaffolds which have both coverage and annotated
3' utr3 for CP site discovery
}
\examples{
library(BSgenome.Mmusculus.UCSC.mm10)
genome <- BSgenome.Mmusculus.UCSC.mm10
data(utr3.mm10)
utr3 <- split(utr3.mm10, seqnames(utr3.mm10), drop = TRUE)
bedgraphs <- system.file("extdata", c(
  "Baf3.extract.bedgraph",
  "UM15.extract.bedgraph"
),
package = "InPAS"
)
tags <- c("Baf3", "UM15")
metadata <- data.frame(
  tag = tags,
  condition = c("Baf3", "UM15"),
  bedgraph_file = bedgraphs
)
outdir <- tempdir()
write.table(metadata,
  file = file.path(outdir, "metadata.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

sqlite_db <- setup_sqlitedb(
  metadata = file.path(
    outdir,
    "metadata.txt"
  ),
  outdir
)
addLockName(filename = tempfile())
coverage <- list()
for (i in seq_along(bedgraphs)) {
      coverage[[tags[i]]] <- get_ssRleCov(
      bedgraph = bedgraphs[i],
      tag = tags[i],
      genome = genome,
      sqlite_db = sqlite_db,
      outdir = outdir,
      chr2exclude = "chrM"
    )
}
get_chromosomes(utr3, sqlite_db)
}
