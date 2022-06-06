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
#' @import methods
#' @import GenomicRanges
#' @importFrom IRanges CharacterList
#' @exportClass TFEAresults
#' @examples
#' res <- system.file("extdata", "res.rds", package="ATACseqTFEA")
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
#' @param object an object of trackStyle
#' @exportMethod show
#' @importFrom S4Vectors head
#' @aliases show,TFEAresults-method
setMethod("show", "TFEAresults", function(object){
  cat("This is an object of TFEAresults\n")
  cat("slot enrichmentScore (", nrow(object@enrichmentScore), "x",
      ncol(object@enrichmentScore), "):\n")
  show(object@enrichmentScore[seq.int(min(c(3,
                                            nrow(object@enrichmentScore)))),
                              seq.int(min(c(5,
                                            ncol(object@enrichmentScore))))])
  cat("slot bindingSites:\n")
  show(object@bindingSites)
  cat("slot motifID (", length(object@motifID), "):\n")
  if(length(object@motifID)>0){
    cat("[[1]] '", names(object@motifID)[1], "'\n")
    show(head(object@motifID[[1]]))
  }else{
    show(object@motifID)
  }
  cat("slot resultsTable:\n")
  show(head(object@resultsTable))
})

#' as("TFEAresults", "data.frame")
#' @name as
#' @family TFEAresults
#' @rdname TFEAresults-class
#' @aliases coerce,TFEAresults,data.frame-method
#' @exportMethod coerce
setAs(from="TFEAresults", to="data.frame", function(from){
  from@resultsTable
})

#' @rdname TFEAresults-class
#' @param \dots Each argument in \dots becomes an slot in the new
#' \code{"TFEAresults"}-class.
#' @return A TFEAresults object.
#' @export
TFEAresults <- function(...){
  new("TFEAresults", ...)
}

#' @rdname TFEAresults-class
#' @export
#' @param x TFEAresults object.
setMethod("$", "TFEAresults", function(x, name) slot(x, name))
#' @rdname TFEAresults-class
#' @param name A literal character string or a name (possibly backtick quoted).
#' @param value value to replace.
#' @export
setReplaceMethod("$", "TFEAresults",
                 function(x, name, value){
                   slot(x, name, check = TRUE) <- value
                   x
                 })


#' @rdname TFEAresults-class
#' @export
#' @param i,j indices specifying elements to extract or replace.
#' @param exact see \link[base]{Extract}
setMethod("[[", "TFEAresults", function(x, i, j, ..., exact=TRUE) slot(x, i))
#' @rdname TFEAresults-class
#' @export
setReplaceMethod("[[", "TFEAresults",
                 function(x, i, ..., value){
                   slot(x, i, check = TRUE) <- value
                   x
                 })


#' @rdname TFEAresults-class
#' @export
#' @aliases getEnrichmentScore
setGeneric("getEnrichmentScore", function(x)
  standardGeneric("getEnrichmentScore"))
#' @rdname TFEAresults-class
#' @export
#' @aliases getEnrichmentScore,TFEAresults-method
setMethod("getEnrichmentScore", "TFEAresults", function(x)
  slot(x, "enrichmentScore"))
#' @rdname TFEAresults-class
#' @export
#' @aliases getBindingSites
setGeneric("getBindingSites", function(x) standardGeneric("getBindingSites"))
#' @rdname TFEAresults-class
#' @export
#' @aliases getBindingSites,TFEAresults-method
setMethod("getBindingSites", "TFEAresults", function(x)
  slot(x, "bindingSites"))
#' @rdname TFEAresults-class
#' @export
#' @aliases getMotifID
setGeneric("getMotifID", function(x) standardGeneric("getMotifID"))
#' @rdname TFEAresults-class
#' @export
#' @aliases getMotifID,TFEAresults-method
setMethod("getMotifID", "TFEAresults", function(x) slot(x, "motifID"))

