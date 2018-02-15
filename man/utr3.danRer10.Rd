\name{utr3.danRer10}
\alias{utr3.danRer10}
\docType{data}
\title{
    3'UTR annotation for danRer10 obtained from utr3Annotation
}
\description{
    3'UTR annotation obtained from utr3Annotation by 
    TxDb.Drerio.UCSC.danRer10.refGene and org.Dr.eg.db
}
\usage{data(utr3.danRer10)}
\format{
    GRanges with slot start holding the start position of the 3'UTR, slot end 
    holding the end position of the 3'UTR, slot names holding transcripts and 
    gene names of 3'UTR, slot seqnames holding the chromosome location where 
    the 3'UTR is located and slot strand for strand of 3'UTR. In addition, the 
    following variables are included.
    \describe{
        \item{\code{feature}}{should be unknown or proximalCP_XXXXXXXX}
    
        \item{\code{id}}{should be utr3 or next.exon.gap}
    
        \item{\code{exon}}{exon id}
    
        \item{\code{transcript}}{trnscript id}
    
        \item{\code{gene}}{entriz gene id}
    
        \item{\code{symbol}}{gene symbol}
    }
}
\details{
    used in the examples
    Annotation data obtained by:
    library(TxDb.Drerio.UCSC.danRer10.refGene)
    
    library(org.Dr.eg.db)
        
    utr3Annotation(TxDb.Drerio.UCSC.danRer10.refGene, "org.Dr.egSYMBOL")
}
\value{
    an object of GRanges.
}

\examples{
    data(utr3.danRer10)
    head(utr3.danRer10)
}
\keyword{datasets}
