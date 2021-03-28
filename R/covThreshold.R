#' Helper function to get intronic regions
#' 
#' @param txdb An object of [GenomicFeatures::TxDb-class]
#'
#' @return GRanges of intronic regions not overlapping exons
#' @import plyranges GenomicRanges  GenomicFeatures
#' @keywords internal
#'



intronRegion <- function(txdb){
    if(!is(txdb, "TxDb")) stop("txdb must be an object of TxDb")
    exons <- exons(txdb, columns=c("exon_id", "gene_id")) %>%
        plyranges::reduce_ranges_directed()
    genes <- genes(txdb) %>%
        plyranges::reduce_ranges_directed()
    dis <- disjoin(c(genes, exons))
    ol <- findOverlaps(dis, exons)
    intronRegion <- dis[-unique(queryHits(ol))]
}

## automatically determining the long_coverage_threshold by 
## non_zero intergenicRegion coverage
## auto determine the short_coverage_threshold by 
## quantile of non_zero intragenicRegion coverage

#' calculate the cutoff threshold of coverage
#'
#' calculate the cutoff threshold of coverage for long and short isoforms
#' 
#' @param coverage Coverage for each sample, output of [coverageFromBedGraph()]
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param txdb An object of [GenomicFeatures::TxDb-class]
#' @param utr3  Output of [utr3Annotation()]
#' @param chr  Chromosome to be used for calculation, default is "chr1"
#' @param hugeData Is this dataset consume too much memory? if it is 
#' TRUE, the coverage will be saved into tempfiles.
#' @param groupList Group list of tag names 
#'
#' @return  A numeric vector
#' @import GenomicRanges plyranges
#' @keywords internal


covThreshold <- function(coverage, genome, txdb, utr3, 
                         chr="chr1", hugeData, groupList){
    if(!is(txdb, "TxDb")) stop("txdb must be an object of TxDb")
    
    totalCov <- totalCoverage(coverage, genome, hugeData, groupList)
    chr1totCov <- lapply(totalCov, "[[", chr)
    N <- length(chr1totCov)
    if(N>1){
        for(i in 2:N){
            chr1totCov[[1]] <- chr1totCov[[1]] + chr1totCov[[i]]
        }
    }
    chr1totCov <- chr1totCov[[1]]
    intronRegion <- intronRegion(txdb) %>% 
        plyranges::filter(seqnames == chr)
    
    utr3Chr1Region <- utr3 %>% 
        plyranges::filter(seqnames == chr & feature=="utr3")
    
    covBg<-function(.cvg, start, end){
        view <- Views(.cvg, start, end)
        view <- viewApply(view, function(.ele) as.integer(.ele))
        view <- unlist(view)
        view <- view[view>=N]
        floor(quantile(view)[2])
    }
    long_coverage_threshold <- covBg(chr1totCov, 
                                     start(intronRegion), 
                                     end(intronRegion))
    short_coverage_threshold <- covBg(chr1totCov, 
                                      start(utr3Chr1Region), 
                                      end(utr3Chr1Region))
    ceiling(c(long_coverage_threshold/N, short_coverage_threshold/N))
}
