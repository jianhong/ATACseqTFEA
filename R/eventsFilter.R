#' Filter the RangedSummarizedExperiment objects
#' @description A helper function to subset the counts object outputed by
#' \link{count5ends}.
#' @param se An \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment} object.
#'  Outputs of \link{count5ends}.
#' @param filter An expression which, when evaluated in the context of
#' assays(se), is a logical vector indicating elements or rows to keep.
#' The expression results for each assay will be combined and use `or` operator
#' to filter the counts assays.
#' @importFrom SummarizedExperiment SummarizedExperiment assays
#' @return A RangedSummarizedExperiment object with assays of
#' count matrix with bindingSites, proximalRegion and distalRegion as
#' column names and bindingSites GRanges object as rowRanges.
#' @export
#' @author Jianhong Ou
#' @examples
#' bam <- system.file("extdata",
#'                    "KD.shift.rep1.bam",
#'                    package="ATACseqTFEA")
#' bsl <- system.file("extdata", "bindingSites.rds",
#'                    package="ATACseqTFEA")
#' bindingSites <- readRDS(bsl)
#' ## get the count regions
#' bsEx <- expandBindingSites(bindingSites)
#' ## count reads by 5'ends
#' res <- count5ends(bam, bindingSites=bindingSites,
#'                   bindingSitesWithGap=bsEx$bindingSitesWithGap,
#'                   bindingSitesWithProximal=bsEx$bindingSitesWithProximal,
#'                   bindingSitesWithProximalAndGap=
#'                       bsEx$bindingSitesWithProximalAndGap,
#'                   bindingSitesWithDistal=bsEx$bindingSitesWithDistal)
#'eventsFilter(res, proximalRegion>0)
#'eventsFilter(res, seqnames(res)=="chr1")
#'eventsFilter(res, sample(c(TRUE, FALSE), length(res), replace=TRUE))
#'eventsFilter(res, "proximalRegion>0")
#'filter <- "proximalRegion>0"
#'eventsFilter(res, filter)
#'filter <- sample(c(TRUE, FALSE), length(res), replace=TRUE)
#'eventsFilter(res, filter)
eventsFilter <- function(se, filter){
  stopifnot(is(se, "RangedSummarizedExperiment"))
  filter0 <- substitute(filter)
  if(is.character(filter0)){
    filter <- parse(text = filter0)
  }else{
    if(is.name(filter0)){
      if(is.character(filter)){
        filter <- parse(text = filter)
      }
    }else{
      filter <- filter0
    }
  }
  if(is.logical(filter)){
    filter <- lapply(assays(se), function(.ele)
      .ele[filter, , drop=FALSE])
  }else{
    filter <- lapply(assays(se), function(.ele)
      as.logical(eval(filter, envir=.ele)))
  }
  filter <- do.call(cbind, filter)
  filter <- rowSums(filter) > 0
  return(se[filter])
}

