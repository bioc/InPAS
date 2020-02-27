seqLen <- function(genome){
    if(!is(genome, "BSgenome")) stop("genome must be an object of BSgenome")
    seqnames <- trimSeqnames(genome)
    seqlengths(genome)[seqnames]
}