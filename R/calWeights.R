#' Calculate the weights for binding score
#' @description Use open score to calculate the weights for the binding score.
#' The open score is calculated by the counts of the proximal region divided by
#' the counts of the distal region. And the counts RangedSummarizedExperiment
#' will be filtered by the Z-score of the open score.
#' The weight is calculated by converting the Z score to the range of 0-1
#' following the normal distribution.
#' @param se An \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#'  object. Outputs of \link{countsNormalization}.
#' @param openscoreZcutoff Open score Z value cutoff value. Default is 0.
#'   Open score is calculated by the count ratio of
#'   proximal site and distal site.
#' @param ... Not used.
#' @return A RangedSummarizedExperiment object with assays of
#' count matrix with bindingSites, proximalRegion and distalRegion as
#' column names and bindingSites GRanges object with the weights as rowRanges.
#' @importFrom S4Vectors isSingleNumber
#' @importFrom stats sd pnorm
#' @importFrom SummarizedExperiment rowData rowData<- assays<-
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
#' ## calculate the weights
#' calWeights(se)
calWeights <- function(se, openscoreZcutoff=0, ...){
  stopifnot(is(se, "RangedSummarizedExperiment"))
  stopifnot(isSingleNumber(openscoreZcutoff))

  ## open score = proximal/distal
  openscore <- getOpenScore(se)

  ## weight = openscore>0?1-p:0
  openscoreZ <- apply(openscore, 2, function(.ele){
    mu <- mean(.ele, na.rm = TRUE)
    std <- sd(.ele, na.rm = TRUE)
    (.ele - mu)/std
  })
  openscoreP <- apply(openscoreZ, 2, function(.ele){
    ifelse(.ele>0, 1-2*pnorm(-abs(.ele)), 0)
  })
  colnames(openscoreP) <- paste0(get("prefix", envir = .globalEnv),
                                 colnames(openscoreP))
  rowData(se) <- cbind(rowData(se), openscoreP)
  ## remove the values that Z-score <openscoreZcutoff
  keep <- rowSums(openscoreZ>=openscoreZcutoff) > 0
  se[keep]
}
