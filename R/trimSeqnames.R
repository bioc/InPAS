#' filter the sequence names from a BSgenome object
#' 
#' only chromosome-level seqnames are kept
#' 
#' @param genome an object of [BSgenome::BSgenome=class]
#' @param removeScaffolds logical(1), whether the scaffolds should be 
#' removed from the genome
#' @import BSgenome
#' @return an character vector of filtered seqnames
#' @keywords internal
#'
#' @examples
#' 1
trimSeqnames <- function(genome, removeScaffolds = FALSE){
    if(!is(genome, "BSgenome")) stop("genome must be an object of BSgenome")
    seqnames <- seqnames(genome)
    if (removeScaffolds)
    {
        seqnames <- seqnames[grepl("^(chr)?(\\d+|[XY])$", seqnames)]
    }
    seqnames
}

