# private utilities
.globalEnv <- new.env(parent = emptyenv())
assign("prefix", "openscoreP_", envir = .globalEnv)

#' @importFrom S4Vectors isSingleInteger
checkInteger <- function(...){
  dots <- list(...)
  if(is.null(names(dots))) names(dots) <- deparse(substitute(...))
  mapply(dots, names(dots), FUN = function(.ele, .para){
    if(!isSingleNumber(.ele)){
      stop(.para, " is not a single numeric value")
    }
    as.integer(.ele)
  })
}

# pseudo log2
pseudoLog2 <- function(x, dividend, divisor, pseudo_count=1){
  log2(x[, dividend]+pseudo_count) - log2(x[, divisor]+pseudo_count)
}

# get score from RangedSummarizedExperiment
getScore <- function(se, dividend, divisor){
  stopifnot(is(se, "RangedSummarizedExperiment"))
  score <- lapply(assays(se), pseudoLog2,
                         dividend=dividend,
                         divisor=divisor)
  score <- do.call(cbind, score)
}
# binding score = proximal/binding
getBindingScore <- function(se){
  getScore(se, dividend='proximalRegion', divisor='bindingSites')
}
# open score = proximal/distal
getOpenScore <- function(se){
  getScore(se, dividend='proximalRegion', divisor='distalRegion')
}
