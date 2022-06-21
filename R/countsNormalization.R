#' Normalize counts by width of count region
#' @description Do normalization by width for counts in binding sites, proximal
#'  and distal regions.
#' @param se An \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#'  object. Outputs of \link{count5ends} or \link{eventsFilter}.
#' @param proximal,distal numeric(1) or integer(1).
#'        bases for open region from binding sites (proximal) and
#'        extended region for background (distal)
#'        of the binding region for aggregate ATAC-seq footprint.
#' @importFrom SummarizedExperiment SummarizedExperiment assays assays<-
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
#' countsNormalization(se, proximal=40, distal=40)
countsNormalization <- function(se, proximal, distal){
  stopifnot(is(se, "RangedSummarizedExperiment"))
  proximal <- checkInteger(proximal)
  distal <- checkInteger(distal)

  bindingSites <- rowRanges(se)
  counts <- assays(se)
  wid <- width(bindingSites)
  names(wid) <- names(bindingSites)
  norm_assay <- data.frame(bindingSites=wid/2,
                           proximalRegion=proximal,
                           distalRegion=distal)
  # rename the norm_assay by counts table
  colnames(norm_assay) <- colnames(counts[[1]])
  MAX_V <- max(c(wid, proximal, distal))
  assays(se) <- lapply(counts, function(.ele){
    round(.ele*MAX_V/norm_assay)
  })
  se
}





