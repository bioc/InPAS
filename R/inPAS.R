#' InPAS: a package for the identification of novel alternative
#' PolyAdenylation Sites (PAS) based on RNA-seq data
#'
#' The InPAS package provides three categories of important functions:
#' parseTxDb, utr3Annotation, coverageCDS, coverageFromBedGraph, coverageRate,
#' getUTR3eSet, utr3UsageEstimation, testUsage, singleSampleAnalyze,
#' singleGroupAnalyze, limmaAnalyze, filterRes, usage4plot, prepare4GSEA,
#' inPAS
#' 
#' @section functions for retrieving 3' UTR annotation:
#' parseTxDb, utr3Annotation
#' @section functions for processing read coverage data:
#' coverageCDS, coverageFromBedGraph, coverageRate
#' @section  functions for alternative polyadenylation site analysis:
#' utr3UsageEstimation, testUsage, singleSampleAnalyze,
#' singleGroupAnalyze, limmaAnalyze, filterRes, usage4plot
#'
#' @docType package
#' @name InPAS
globalVariables(c(".", "Predicted_Proximal_APA", "X2",
                  "X3", "dup.group", "end.utr3.last",
                  "exon_name","exon_rank", "feature",
                  "gene", "start.utr3.last", "transcript",
                  "truncated", "tx_name"))