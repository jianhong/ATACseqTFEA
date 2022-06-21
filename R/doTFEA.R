#' Transcription factor enrichment analysis
#' @description Transcription factor enrichment analysis for the filtered
#'  output of \link{DBscore}
#' @param se An \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#'  object. Filtered outputs of \link{DBscore}.
#' @param ... Not used.
#' @importFrom stats p.adjust
#' @importFrom pracma erf
#' @importFrom BiocGenerics start end width start<- end<-
#' @importFrom S4Vectors mcols
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom Matrix Matrix rowSums
#' @return A \link{TFEAresults} object.
#' @export
#' @author Jianhong Ou
#' @examples
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
#' ## get the count regions
#' bsEx <- expandBindingSites(bindingSites)
#' ## count reads by 5'ends
#' res <- count5ends(c(bamExp, bamCtl),
#'                   positive=0L, negative=0L,
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
#' se <- getWeightedBindingScore(se)
#' design <- cbind(CTL=1, EXPvsCTL=c(1, 1, 0, 0))
#' rownames(design) <- colnames(se)
#' counts <- DBscore(se, design=design, coef="EXPvsCTL")
#' doTFEA(counts)
doTFEA <- function(se, ...){
  stopifnot(is(se, "RangedSummarizedExperiment"))
  bindingSites <- rowRanges(se)
  stopifnot("Columne 'motif' and 't' is not in rowRanges of inputs"=
              all(c("motif", "t") %in% colnames(mcols(bindingSites))))

  ## sort the bindingSites matrix by t value
  bindingSites <- bindingSites[order(bindingSites$t, decreasing = TRUE)]

  ## Calculate enrichment
  ## walk down the list L,
  ## incrementing the running sum statistic by sqrt((N-Nh)/Nh)
  ## when we encounter a gene in S and decrementing by sqrt(Nh/(N-Nh))
  ## if the gene is not in S, where N is the number of genes in the list L,
  ## and Nh is the number of gene in the gene set S.

  ## create a motif matrix.
  motifs <- unique(unlist(bindingSites$motif))
  hits <- Matrix(0, nrow = length(motifs), ncol = length(bindingSites))
  rownames(hits) <- motifs
  motifID <- data.frame(i=rep(seq_along(bindingSites),
                              lengths(bindingSites$motif)),
                        j=unlist(bindingSites$motif))
  motifID <- split(motifID$i, motifID$j)
  for(id in names(motifID)){
    hits[id, motifID[[id]]] <- 1
  }

  cal_ES <- function(hits, N, Nh){
    miss <- t(apply(hits==0, 1, cumsum))
    hits <- t(apply(hits, 1, cumsum))
    hits <- hits/Nh
    miss <- miss/(N-Nh)
    ES <- hits - miss
  }

  N <- ncol(hits)
  Nh <- rowSums(hits)

  ES <- cal_ES(hits, N, Nh)

  n <- (N-Nh)*Nh/N
  k <- seq(from=1, to=10, by = 1) ## 10 is enough.

  ESmax <- apply(ES, 1, max, na.rm=TRUE)
  ESmin <- apply(ES, 1, min, na.rm=TRUE)
  ESm <- ifelse(abs(ESmax)>abs(ESmin), ESmax, ESmin)

  p <- mapply(ESm, n, FUN = function(.esm, .n){
    -2*sum((-1)^k*exp(-2*k^2*.esm^2*.n))
  })

  ## expect ES
  k <- seq(1, 10000)
  EES <- vapply(n, FUN = function(.n){
    8*sum((-1)^(k+1)*(exp(-2*k^2*.n)/4 -
                        (sqrt(2*pi)*erf(k*sqrt(2*.n)))/
                        (16*k*sqrt(.n))))
  }, FUN.VALUE = 0.0)
  ## normalized ES = ES/EES
  NES <- ESm/abs(EES)

  res <- data.frame(TF = names(ESm),
                    enrichmentScore=ESm,
                    normalizedEnrichmentScore=NES,
                    p_value = p,
                    adjPval = p.adjust(p, method = "BH"))

  new("TFEAresults",
      enrichmentScore=ES,
      bindingSites = bindingSites,
      motifID = motifID,
      resultsTable = res)
}
