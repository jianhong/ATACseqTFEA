#' The methods for \link{TFEAresults-class}
#' @description The assessment and replacement methods for
#' \link{TFEAresults-class}
#' @name TFEAresults-methods
#' @family TFEAresults
#' @rdname TFEAresults-methods
#' @aliases coerce,TFEAresults,data.frame-method
#' @aliases as,TFEAresults,data.frame-method
#' @exportMethod coerce
#' @examples
#' res <- readRDS(system.file("extdata", "res.rds", package="ATACseqTFEA"))
#' as(res, "data.frame")
setAs(from="TFEAresults", to="data.frame", function(from){
  from@resultsTable
})

#' @rdname TFEAresults-methods
#' @param object an object of TFEAresults
#' @exportMethod show
#' @importFrom S4Vectors head
#' @aliases show,TFEAresults-method
#' @examples
#' res
setMethod("show", "TFEAresults", function(object){
  cat("This is an object of TFEAresults with \n")
  cat("slot enrichmentScore ( matrix: ", nrow(object@enrichmentScore), "x",
      ncol(object@enrichmentScore), "), \n")
  cat("slot bindingSites ( GRanges object with ",
      length(object@bindingSites), " ranges and ",
      ncol(mcols(object@bindingSites)), " metadata columns",
      "), \n")
  cat("slot motifID ( a list of the positions of binding sites for ",
      length(object@motifID), "TFs ), and \n")
  cat("slot resultsTable (", nrow(object@resultsTable), " x ",
      ncol(object@resultsTable), "). Here is the top 2 rows:\n")
  show(head(object@resultsTable, n=2L))
})

#' @rdname TFEAresults-methods
#' @export
#' @param x TFEAresults object.
setMethod("$", "TFEAresults", function(x, name) slot(x, name))
#' @rdname TFEAresults-methods
#' @param name A literal character string or a name (possibly backtick quoted).
#' @param value value to replace.
#' @export
#' @examples
#' head(res$resultsTable)
setReplaceMethod("$", "TFEAresults",
                 function(x, name, value){
                   slot(x, name, check = TRUE) <- value
                   x
                 })


#' @rdname TFEAresults-methods
#' @export
#' @param i,j indices specifying elements to extract or replace.
#' @param \dots Named or unnamed arguments to form a signature.
#' @param exact see \link[base]{Extract}
setMethod("[[", "TFEAresults", function(x, i, j, ..., exact=TRUE) slot(x, i))
#' @rdname TFEAresults-methods
#' @export
#' @examples
#' head(res[["resultsTable"]])
setReplaceMethod("[[", "TFEAresults",
                 function(x, i, ..., value){
                   slot(x, i, check = TRUE) <- value
                   x
                 })

#' @rdname TFEAresults-methods
#' @export
#' @aliases getEnrichmentScore,TFEAresults-method
#' @examples
#' head(getEnrichmentScore(res))
setMethod("getEnrichmentScore", "TFEAresults", function(x)
  slot(x, "enrichmentScore"))

#' @rdname TFEAresults-methods
#' @export
#' @importFrom S4Vectors isSingleString
#' @aliases getBindingSites,TFEAresults-method
setMethod("getBindingSites", "TFEAresults", function(x, TF){
  if(missing(TF)){
    slot(x, "bindingSites")
  }else{
    stopifnot(isSingleString(TF))
    bs <- slot(x, "bindingSites")
    bs <- bs[vapply(bs$motif,
                    FUN=function(.ele) TF %in% .ele,
                    FUN.VALUE = logical(1L))]
    bs$score <- mapply(bs$score, bs$motif, FUN = function(score, motif){
      max(score[motif==TF])
    })
    bs
  }
})

#' @rdname TFEAresults-methods
#' @export
#' @aliases getMotifID,TFEAresults-method
setMethod("getMotifID", "TFEAresults", function(x) slot(x, "motifID"))

