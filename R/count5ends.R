#' Prepare counts matrix for enrichment analysis
#' @description Prepare the counts matrix by 5'end of reads.
#' @param bam A character vector indicates the file names of the
#' bams or an object of \link[Rsamtools:BamFile]{BamFile}.
#' @param index The names of the index file of the 'BAM' file
#' being processed; This is given without the '.bai' extension.
#' @param yieldSize Number of records to yield each time the file is read.
#' See \link[Rsamtools:BamFile]{BamFile} for details.
#' @param positive,negative integer(1). the size to be shift for
#' positive/negative strand.
#' If the bam file is 5'end shifed files, please set the parameter to 0.
#' @param bindingSites A object of
#' \link[GenomicRanges:GRanges-class]{GenomicRanges} indicates
#' candidate binding sites. The \link{prepareBindingSites} function
#' is a helper function to generate the binding sites.
#' Users can also use other software for example fimo to generate the list.
#' @param bindingSitesWithGap bindingSites with gaps and in both ends,
#' @param bindingSitesWithProximal bindingSites with gaps and proximal region
#' in both ends,
#' @param bindingSitesWithProximalAndGap bindingSites with gaps, and then
#'  proximal and gaps in both ends,
#' @param bindingSitesWithDistal bindingSites with gap, proximal, gap and
#' distal regions.
#' @importFrom GenomicAlignments readGAlignments summarizeOverlaps
#' cigar cigarNarrow cigarQNarrow qwidth
#' @importFrom Rsamtools ScanBamParam BamFile asMates
#' @importFrom BiocGenerics start end width start<- end<- strand
#' @importFrom SummarizedExperiment SummarizedExperiment assays
#' @importFrom S4Vectors isSingleNumber
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
#' res <- count5ends(bam, positive=0L, negative=0L,
#'                   bindingSites=bindingSites,
#'                   bindingSitesWithGap=bsEx$bindingSitesWithGap,
#'                   bindingSitesWithProximal=bsEx$bindingSitesWithProximal,
#'                   bindingSitesWithProximalAndGap=
#'                       bsEx$bindingSitesWithProximalAndGap,
#'                   bindingSitesWithDistal=bsEx$bindingSitesWithDistal)
#' head(res)
count5ends <- function(bam, index=bam,
                       yieldSize=100000,
                       positive=4L, negative=5L,
                       bindingSites,
                       bindingSitesWithGap,
                       bindingSitesWithProximal,
                       bindingSitesWithProximalAndGap,
                       bindingSitesWithDistal){
  if(missing(bam)){
    stop("bam is requred!")
  }
  stopifnot(isSingleNumber(positive))
  stopifnot(isSingleNumber(negative))
  positive <- as.integer(positive)
  negative <- as.integer(negative)
  stopifnot("bindingSites must be an GRanges object"=
              is(bindingSites, "GRanges"))
  stopifnot("bindingSites must contain mcols 'motif'"=
              "motif" %in% colnames(mcols(bindingSites)))
  ## check rowRanges order
  stopifnot(all(distance(bindingSites, bindingSitesWithProximal)==0))
  stopifnot(all(distance(bindingSitesWithProximal, bindingSitesWithDistal)==0))

  if(length(bam)>1){
    stopifnot(length(index)==length(bam))
    counts <- mapply(function(a, b) {
      count5ends(bam=a, index=b,
                 positive=positive, negative=negative,
                 bindingSites = bindingSites,
                 bindingSitesWithGap=bindingSitesWithGap,
                 bindingSitesWithProximal=bindingSitesWithProximal,
                 bindingSitesWithProximalAndGap=
                   bindingSitesWithProximalAndGap,
                 bindingSitesWithDistal=bindingSitesWithDistal)
    }, bam, index, SIMPLIFY = FALSE)
    counts <- lapply(counts, function(.ele) assays(.ele)[[1]])
    names(counts) <- sub(".bam", "", make.names(basename(bam)))
    return(SummarizedExperiment(assays=counts,
                                rowRanges = unname(bindingSites)))
  }
  if(is.character(bam)){
    stopifnot(file.exists(bam))
  }
  if(is(bam, "BamFile")){
    stopifnot("BamFile must be treated as single ends"=!asMates(bam))
    bamfile <- bam
  }else{
    bamfile <- BamFile(bam, index=index, yieldSize=yieldSize, asMates = FALSE)
  }

  open(bamfile)
  on.exit(close(bamfile))
  counts <- data.frame(bs=rep(0, length(bindingSites)),
                       pro.gap=rep(0, length(bindingSitesWithGap)),
                       pro=rep(0, length(bindingSitesWithProximal)),
                       dis.gap=rep(0, length(bindingSitesWithProximalAndGap)),
                       dis=rep(0, length(bindingSitesWithDistal)))
  while (length(chunk0 <- readGAlignments(bamfile))) {
    ## shift 5' end
    if(positive!=0 && negative!=0){
      strands <- as.character(strand(chunk0)) == "-"
      ns <- ifelse(strands, negative, positive)
      cigars <- cigar(chunk0)
      cigars <- as.character(cigarNarrow(cigars))
      cigars <- cigarQNarrow(cigars,
                             start=ifelse(strands, 1, positive+1),
                             end=ifelse(strands, -negative-1, -1))
      chunk0@start <- start(chunk0) + attributes(cigars)$rshift
      chunk0@cigar <- as.character(cigars)
    }

    ## keep 5' end only
    chunk0 <- as(chunk0, "GRanges")
    chunk0 <- promoters(chunk0, upstream = 0, downstream = 1)
    ## counts of binding sites
    so.bs <- countOverlaps(bindingSites, chunk0, ignore.strand=TRUE)
    ## counts of binding sites + gap
    so.pro.gap <- countOverlaps(bindingSitesWithGap, chunk0,
                                ignore.strand=TRUE)
    ## counts of binding sites + proximal
    so.pro <- countOverlaps(bindingSitesWithProximal, chunk0,
                            ignore.strand=TRUE)
    ## counts of binding sites + gap + proximal + gap
    so.dis.gap <- countOverlaps(bindingSitesWithProximalAndGap,
                                chunk0, ignore.strand=TRUE)
    ## counts of binding sites + proximal + distal
    so.dis <- countOverlaps(bindingSitesWithDistal, chunk0, ignore.strand=TRUE)
    counts <- counts + data.frame(bs=so.bs,
                                  pro.gap=so.pro.gap,
                                  pro=so.pro,
                                  dis.gap=so.dis.gap,
                                  dis=so.dis)
  }
  close(bamfile)
  on.exit()
  counts$dis <- counts$dis - counts$dis.gap
  counts$pro <- counts$pro - counts$pro.gap
  counts$pro.gap <- NULL
  counts$dis.gap <- NULL
  stopifnot(all(counts$dis>=0))
  stopifnot(all(counts$pro>=0))
  colnames(counts) <- c("bindingSites", "proximalRegion", "distalRegion")
  SummarizedExperiment(assays=list(counts=counts),
                       rowRanges = unname(bindingSites))
}
