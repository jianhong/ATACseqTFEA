#' @rdname TFEAresults-methods
#' @export
#' @aliases getEnrichmentScore
setGeneric("getEnrichmentScore", function(x)
  standardGeneric("getEnrichmentScore"))

#' @rdname TFEAresults-methods
#' @param TF Transcription factor
#' @export
#' @aliases getBindingSites
setGeneric("getBindingSites", function(x, TF)
  standardGeneric("getBindingSites"))

#' @rdname TFEAresults-methods
#' @export
#' @aliases getMotifID
setGeneric("getMotifID", function(x) standardGeneric("getMotifID"))
