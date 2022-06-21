#' Class \code{"TFEAresults"}
#' @description An object of class \code{"TFEAresults"}
#'   represents the results of \link{TFEA}.
#' @aliases TFEAresults
#' @rdname TFEAresults-class
#' @slot enrichmentScore \code{"numeric \link[Matrix]{Matrix}"}, specify the
#'  enrichment score for each transcription factor (TF). Every row represents
#'  a TF. The columns represents the accumulated enrichment score for that rank.
#' @slot bindingSites \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}}
#' object. It is keep same length and order as the columns in enrichmentScore.
#' @slot motifID \code{"list"}. The ranks of binding sites for each TF.
#' @slot resultsTable \code{"data.frame"}. The data frame contains the
#' summarized enrichment score, the p-value, and adjuct p-value for each TF.
#' @importFrom methods setClass representation prototype setMethod setAs
#' setReplaceMethod setGeneric as is new slot slot<- show coerce
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges CharacterList
#' @exportClass TFEAresults
#' @examples
#' res <- readRDS(system.file("extdata", "res.rds", package="ATACseqTFEA"))
#' res

setClass("TFEAresults",
         representation = representation(
           enrichmentScore="matrix",
           bindingSites="GRanges",
           motifID="list",
           resultsTable="data.frame"
         ),
         prototype = prototype(
           enrichmentScore=matrix(),
           bindingSites=GRanges(motif=CharacterList(), score=list()),
           motifID=list(),
           resultsTable=data.frame(TF=character(),
                                   enrichmentScore=numeric(),
                                   p_value=numeric(),
                                   adjPval=numeric())
         ),
         validity = function(object){
           if(!all(c("TF", "enrichmentScore", 'p_value', "adjPval") %in%
                   colnames(object@resultsTable))){
             return("resultsTable must contain TF, enrichmentScore, p_value,
                    and adjPval.")
           }
           if(nrow(object@enrichmentScore)!=length(object@motifID)){
             return("The length of motifID should be identical with the number
                    of row of enrichmentScore.")
           }
           if(ncol(object@enrichmentScore)!=length(object@bindingSites)){
             return("The length of bindingSites should be identical with the
                    number of column of enrichmentScore.")
           }
           return(TRUE)
         })

#' @rdname TFEAresults-class
#' @param \dots Each argument in \dots becomes an slot in the new
#' \code{"TFEAresults"}-class.
#' @return A TFEAresults object.
#' @export
TFEAresults <- function(...){
  new("TFEAresults", ...)
}
