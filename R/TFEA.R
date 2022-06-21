#' Transcription factor enrichment analysis
#' @description Transcription factor enrichment analysis for
#' ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing).
#' We treat all the binding sites for one TF as a TF set and
#' all the open regions as features for random walking.
#' @param bamExp A vector of characters indicates the file names of
#' experiment bams. The bam file must be the one with shifted reads.
#' @param bamCtl A vector of characters indicates the file names of
#' control bams. The bam file must be the one with shifted reads.
#' @param indexExp,indexCtl The names of the index file of the 'BAM' file
#' being processed; This is given without the '.bai' extension.
#' @param positive,negative integer(1). the size to be shift for
#' positive/negative strand.
#' If the bam file is 5'end shifed files, please set the parameter to 0.
#' @param bindingSites A object of
#' \link[GenomicRanges:GRanges-class]{GenomicRanges} indicates
#' candidate binding sites. The \link{prepareBindingSites} function
#' is a helper function to generate the binding sites.
#' Users can also use other software for example fimo to generate the list.
#' @param proximal,distal numeric(1) or integer(1).
#'        bases for open region from binding sites (proximal) and
#'        extended region for background (distal)
#'        of the binding region for aggregate ATAC-seq footprint.
#' @param gap numeric(1) or integer(1). bases for gaps among binding sites,
#'            proximal, and distal. default is 10L.
#' @param filter An expression which, when evaluated in the context of
#' assays(se), is a logical vector indicating elements or rows to keep.
#' The expression results for each assay will be combined and use `or` operator
#' to filter the counts assays.
#' @param openscoreZcutoff Open score Z value cutoff value. Default is 0.
#'   Open score is calculated by the count ratio of
#'   proximal site and distal site.
#' @param bindingScorePvalCutoff,bindingScoreLog2FCcutoff Binding score cutoff
#'  values. Default is 1 and 0. Binding score is calculated by the count ratio
#'  of proximal site and binding site. The cutoff values are used to decrease
#'  the total number of binding site for ranking. Increasing the `log2FCcutoff`
#'  value and decreasing the P-value cutoff value can greatly decrease the
#'  memory cost and computing time by decreasing the total binding sites.
#' @importFrom stats p.adjust pnorm
#' @importFrom GenomicAlignments readGAlignments summarizeOverlaps
#' @importFrom Rsamtools ScanBamParam BamFile
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats sd
#' @importFrom pracma erf
#' @importFrom BiocGenerics start end width start<- end<-
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowRanges
#' @importFrom Matrix Matrix rowSums
#' @return A \link{TFEAresults} object.
#' @export
#' @author Jianhong Ou
#' @examples
#'
#' bamExp <- system.file("extdata",
#'                       c("KD.shift.rep1.bam",
#'                         "KD.shift.rep2.bam"),
#'                       package="ATACseqTFEA")
#' bamCtl <- system.file("extdata",
#'                       c("WT.shift.rep1.bam",
#'                         "WT.shift.rep2.bam"),
#'                       package="ATACseqTFEA")
#' bsl <- system.file("extdata", "bindingSites.rds",
#'                    package="ATACseqTFEA")
#' bindingSites <- readRDS(bsl)
#' res <- TFEA(bamExp, bamCtl, bindingSites=bindingSites,
#'             positive=0, negative=0)
#' res
TFEA <- function(bamExp, bamCtl,
                 indexExp=bamExp, indexCtl=bamCtl,
                 positive=4L, negative=5L,
                 bindingSites,
                 proximal=40L, distal=proximal, gap=10L,
                 filter="proximalRegion>0",
                 openscoreZcutoff=0,
                 bindingScoreLog2FCcutoff=0,
                 bindingScorePvalCutoff=1){
  stopifnot("bindingSites must be an GRanges object"=
              is(bindingSites, "GRanges"))
  stopifnot("bindingSites must contain mcols 'motif'"=
              "motif" %in% colnames(mcols(bindingSites)))
  proximal <- checkInteger(proximal)
  distal <- checkInteger(distal)
  gap <- checkInteger(gap)
  stopifnot(isSingleNumber(openscoreZcutoff))
  stopifnot(isSingleNumber(bindingScoreLog2FCcutoff))
  stopifnot(isSingleNumber(bindingScorePvalCutoff))
  stopifnot(proximal>10 && distal>10)
  if(missing(bamExp)){
    stop("bamExp is required.")
  }
  if(missing(bamCtl)){
    stop("bamCtl is required.")
  }

  ## get design table
  bams <- data.frame(bamfiles = bamExp, index = indexExp,
                     group="Exp", stringsAsFactors = FALSE)
  bams <- rbind(bams,
                data.frame(bamfiles = bamCtl, index = indexCtl,
                           group="Ctl", stringsAsFactors = FALSE))
  bams$sample <- sub(".bam", "", make.names(basename(bams$bamfiles)))

  ## get the count regions
  exbs <- expandBindingSites(bindingSites=bindingSites,
                             proximal=proximal,
                             distal=distal,
                             gap=gap)
  ## count reads by 5'ends
  counts <- count5ends(bam=bams$bamfiles, index= bams$index,
                       positive=positive, negative=negative,
                       bindingSites = bindingSites,
                       bindingSitesWithGap=exbs$bindingSitesWithGap,
                       bindingSitesWithProximal=exbs$bindingSitesWithProximal,
                       bindingSitesWithProximalAndGap=
                         exbs$bindingSitesWithProximalAndGap,
                       bindingSitesWithDistal=exbs$bindingSitesWithDistal)
  names(assays(counts)) <- bams$sample

  ## filter 0 counts in proximal
  counts <- eventsFilter(counts, filter = filter)

  ## normalize counts by width of count region
  ## normalize by width
  counts <- countsNormalization(counts, proximal = proximal, distal = distal)

  ## get the weighted binding scores
  counts <- getWeightedBindingScore(counts, openscoreZcutoff = openscoreZcutoff)

  ## get the differential analysis by limma for the binding score
  design <- cbind(CTL=1, EXPvsCTL=bams$group=="Exp")
  counts <- DBscore(counts, design=design, coef="EXPvsCTL")

  ## remove the log2 foldchange smaller than bindingScoreLog2FCcutoff
  ## and pvalue greater than bindingScorePvalCutoff
  ## otherwise dataset is too large
  keep <- rowRanges(counts)$P.Value<=bindingScorePvalCutoff &
    abs(rowRanges(counts)$logFC)>=bindingScoreLog2FCcutoff
  counts <- eventsFilter(counts, keep)
  stopifnot(
    "bindingScorePvalCutoff & bindingScoreLog2FCcutoff is too stringency."=
      length(counts)>0)

  ## do TFEA
  doTFEA(counts)
}



