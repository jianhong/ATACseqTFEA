#' Calculate the weighted binding score
#' @description Use user predefined weight to get the weighted binding score
#' or use open score to weight the binding score.
#' The open score is calculated by the counts of proximal region divided by
#' the counts of distal region.
#' The binding score is calculated by the counts of proximal region divided by
#' the counts of binding region. This value is the measure of avoidance of
#' reads in the binding sites.
#' @param se An \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#'  object. Outputs of \link{countsNormalization}.
#' @param weight If NA, the weight will be calculated by the open score.
#'  See \link{calWeights}.
#'  User can define the weight by a matrix or numeric vector.
#' @param ... The parameters will be passed to \link{calWeights}.
#' @return A RangedSummarizedExperiment object with assays of
#' count matrix with bindingSites, proximalRegion and distalRegion as
#' column names and bindingSites GRanges object as rowRanges.
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges
#' rowData assays<-
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
#' res <- count5ends(bam, positive=0L, negative=0L,
#'                   bindingSites=bindingSites,
#'                   bindingSitesWithGap=bsEx$bindingSitesWithGap,
#'                   bindingSitesWithProximal=bsEx$bindingSitesWithProximal,
#'                   bindingSitesWithProximalAndGap=
#'                       bsEx$bindingSitesWithProximalAndGap,
#'                   bindingSitesWithDistal=bsEx$bindingSitesWithDistal)
#' ## filter 0 counts in proximal
#' se <- eventsFilter(res, proximalRegion>0)
#' ## normalize counts by width of count region
#' se <- countsNormalization(se, proximal=40, distal=40)
#' ## get the weighted binding scores
#' getWeightedBindingScore(se)
getWeightedBindingScore <- function(se, weight=NA, ...){
  stopifnot(is(se, "RangedSummarizedExperiment"))
  if(is.na(weight) | is.null(weight)){
    pre <- get("prefix", envir = .globalEnv)
    cn <- colnames(rowData(se))
    cn <- cn[grepl(pre, cn)]
    if(length(cn)<1){## calculate weights
      se <- calWeights(se, ...)
      cn <- colnames(rowData(se))
      cn <- cn[grepl(pre, cn)]
      if(length(cn)<1){
        stop("Can not get the weights.")
      }
    }
    weight <- as.matrix(rowData(se)[, cn, drop=FALSE])
    colnames(weight) <- sub(pre, "", colnames(weight), fixed = TRUE)
    rowData(se) <- rowData(se)[, !colnames(rowData(se)) %in% cn]
  }

  ## binding score = proximal/binding
  bindingscore <- getBindingScore(se)

  if(is.null(dim(weight))){
    if(length(weight)!=length(se)){
      if(length(weight)!=1){
        stop("The length of 'weight' must equal to the length of inputs.")
      }else{
        weight <- rep(weight, length(se))
      }
    }
  }else{
    if(nrow(weight)!=nrow(se)){
      stop("The row number of 'weight' must equal to the length of inputs.")
    }
    if(ncol(weight)!=length(assays(se))){
      stop("Then column number of 'weight'",
           " must equeal to the length of assays.")
    }
    if(all(colnames(bindingscore) %in% colnames(weight))){
      weight <- weight[, colnames(weight), drop=FALSE]
    }
    if(!identical(rownames(bindingscore), rownames(weight))){
      rownames(weight) <- rownames(bindingscore)
    }
  }

  ## following gaussion distribution, should we use Skewness to value it?
  adj.bindingscore <- bindingscore * weight
  se <- SummarizedExperiment(assays = list(bindingScore=adj.bindingscore),
                             rowRanges = rowRanges(se))

  ## remove zero sample variances
  keep <- apply(adj.bindingscore, 1, unique)
  keep <- lengths(keep)>1
  se[keep]
}
