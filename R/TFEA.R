#' Transcription factor enrichment analysis
#' @description Transcription factor enrichment analysis for
#' ATAC-seq (Assay for Transposase-Accessible Chromatin using sequencing).
#' @param bamExp A vector of characters indicates the file names of
#' experiment bams. The bam file must be the one with shifted reads.
#' @param bamCtl A vector of characters indicates the file names of
#' control bams. The bam file must be the one with shifted reads.
#' @param indexExp,indexCtl The names of the index file of the 'BAM' file
#' being processed; This is given without the '.bai' extension.
#' @param bindingSites A object of
#' \link[GenomicRanges:GRanges-class]{GenomicRanges} indicates
#' candidate binding sites. The \link{prepareBindingSites} function
#' is a helper function to generate the binding sites.
#' Users can also use other software for example fimo to generate the list.
#' @param proximal,distal numeric(1) or integer(1).
#'        basepair for open region from binding sites (proximal) and
#'        extented region for background (distal)
#'        of the binding region for aggregate ATAC-seq footprint.
#' @param gap numeric(1) or integer(1). basepair for gaps among binding sites,
#'            proximal, and distal. default is 10L.
#' @param openscoreZcutoff Open score Z value cutoff value. Default is 0.
#'   Open score is calculated by the count ratio of
#'   proximal site and distal site.
#' @param bindingScorePvalCutoff,bindingScoreLog2FCcutoff Binding score cutoff
#'   values. Default is 1 and 0. Binding score is calculated by the count ratio of
#'   proximal site and binding site. The cutoff values are used to decrease
#'   the total number of binding site for ranking. Increasing the log2FCcutoff
#'   value and decreasing the P-value cutoff value can greatly decrease the
#'   memory cost and computing time by decreasing the total binding sites.
#' @importFrom stats p.adjust pnorm
#' @importFrom GenomicAlignments readGAlignments summarizeOverlaps
#' @importFrom Rsamtools ScanBamParam BamFile
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats sd
#' @importFrom pracma erf
#' @import Matrix
#' @import GenomicRanges
#' @return A \link{TFEAresults} object.
#' @export
#' @author Jianhong Ou
#' @examples
#'
#'bamExp <- system.file("extdata",
#'                      c("KD.shift.rep1.bam",
#'                        "KD.shift.rep2.bam"),
#'                      package="ATACseqTFEA")
#'bamCtl <- system.file("extdata",
#'                      c("WT.shift.rep1.bam",
#'                        "WT.shift.rep2.bam"),
#'                      package="ATACseqTFEA")
#'bsl <- system.file("extdata", "bindingSites.rds",
#'                   package="ATACseqTFEA")
#'bindingSites <- readRDS(bsl)
#'res <- TFEA(bamExp, bamCtl, bindingSites=bindingSites)
#'res
TFEA <- function(bamExp, bamCtl, indexExp=bamExp, indexCtl=bamCtl,
                 bindingSites,
                 proximal=40L, distal=proximal, gap=10L,
                 openscoreZcutoff=0,
                 bindingScoreLog2FCcutoff=0,
                 bindingScorePvalCutoff=1){
  stopifnot("bindingSites must be an GRanges object"=
              is(bindingSites, "GRanges"))
  stopifnot("bindingSites must contain mcols 'motif'"=
              "motif" %in% colnames(mcols(bindingSites)))
  stopifnot(is.numeric(proximal))
  stopifnot(is.numeric(distal))
  stopifnot(is.numeric(gap))
  proximal <- as.integer(proximal)
  distal <- as.integer(distal)
  gap <- as.integer(gap)
  stopifnot(proximal>10 && distal>10)
  if(missing(bamExp)){
    stop("bamExp is required.")
  }
  if(missing(bamCtl)){
    stop("bamCtl is required.")
  }
  bindingSites.with.proximal <- bindingSites.with.distal <- bindingSites
  bindingSites.with.proximal.gap <- bindingSites.with.distal.gap <- bindingSites
  start(bindingSites.with.proximal) <- start(bindingSites) - proximal - gap
  end(bindingSites.with.proximal) <- end(bindingSites) + proximal + gap
  start(bindingSites.with.distal) <-
    start(bindingSites) - proximal - 2*gap - distal
  end(bindingSites.with.distal) <-
    end(bindingSites) + proximal + distal + 2*gap
  start(bindingSites.with.proximal.gap) <- start(bindingSites) - gap
  end(bindingSites.with.proximal.gap) <- end(bindingSites) + gap
  start(bindingSites.with.distal.gap) <-
    start(bindingSites) - proximal - 2*gap
  end(bindingSites.with.distal.gap) <- end(bindingSites) + proximal + 2*gap

  ## get total reads counts
  bams <- data.frame(bamfiles = bamExp, index = indexExp,
                     group="Exp", stringsAsFactors = FALSE)
  bams <- rbind(bams,
                data.frame(bamfiles = bamCtl, index = indexCtl,
                           group="Ctl", stringsAsFactors = FALSE))
  bams$sample <- sub(".bam", "", make.names(basename(bams$bamfiles)))
  ## count reads by 5'ends
  count5ends <- function(bam, index, chunk=100000, binding,
                         pro, dis, pro.gap, dis.gap){
    ## check rowRanges order
    stopifnot(all(distance(binding, pro)==0))
    stopifnot(all(distance(pro, dis)==0))
    bamfile <- BamFile(bam, index=index, yieldSize=chunk, asMates = FALSE)
    open(bamfile)
    on.exit(close(bamfile))
    counts <- data.frame(bs=rep(0, length(binding)),
                         pro.gap=rep(0, length(pro.gap)),
                         pro=rep(0, length(pro)),
                         dis.gap=rep(0, length(dis.gap)),
                         dis=rep(0, length(dis)))
    while (length(chunk0 <- readGAlignments(bamfile))) {
      ## keep 5' end only
      chunk0 <- as(chunk0, "GRanges")
      chunk0 <- promoters(chunk0, upstream = 0, downstream = 1)
      ## counts of binding sites
      so.bs <- countOverlaps(binding, chunk0, ignore.strand=TRUE)
      ## counts of binding sites + gap
      so.pro.gap <- countOverlaps(pro.gap, chunk0, ignore.strand=TRUE)
      ## counts of binding sites + proximal
      so.pro <- countOverlaps(pro, chunk0, ignore.strand=TRUE)
      ## counts of binding sites + gap + proximal + gap
      so.dis.gap <- countOverlaps(dis.gap, chunk0, ignore.strand=TRUE)
      ## counts of binding sites + proximal + distal
      so.dis <- countOverlaps(dis, chunk0, ignore.strand=TRUE)
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
    #rownames(counts) <- names(binding)
    counts ## data frame with colnames bs, pro and dis
  }
  counts <- mapply(function(a, b) {
    count5ends(bam=a, index=b, binding = bindingSites,
               pro=bindingSites.with.proximal,
               dis=bindingSites.with.distal,
               pro.gap=bindingSites.with.proximal.gap,
               dis.gap=bindingSites.with.distal.gap)
  }, bams$bamfiles, bams$index, SIMPLIFY = FALSE)
  names(counts) <- bams$sample

  ## filter 0 counts in proximal
  keep <- lapply(counts, function(.ele) .ele$pro>0)
  keep <- do.call(cbind, keep)
  keep <- rowSums(keep) > 0
  counts <- lapply(counts, function(.ele) .ele[keep, , drop=FALSE])
  bindingSites <- bindingSites[keep]

  ## normalize counts by width of count region
  ## normalize by width
  wid <- width(bindingSites)
  names(wid) <- names(bindingSites)
  norm.counts <- lapply(counts, function(.ele){
    round(.ele*max(c(wid, proximal, distal))/data.frame(bs=wid/2,
                                                        pro=proximal,
                                                        dis=distal))
  })

  ## open score = proximal/distal
  openscore <-
    do.call(cbind, lapply(norm.counts, function(.ele)
      log2(.ele[, "pro"]+1) - log2(.ele[, "dis"]+1)))

  ## binding score = proximal/binding
  bindingscore <-
    do.call(cbind, lapply(norm.counts, function(.ele)
      log2(.ele[, "pro"]+1) - log2(.ele[, "bs"]+1)))

  ## weight = openscore>0?1-p:0
  openscoreZ <- apply(openscore, 2, function(.ele){
    mu <- mean(.ele, na.rm = TRUE)
    std <- sd(.ele, na.rm = TRUE)
    (.ele - mu)/std
  })
  openscoreP <- apply(openscoreZ, 2, function(.ele){
    ifelse(.ele>0, 1-2*pnorm(-abs(.ele)), 0)
  })
  ## remove the values that Z-score <openscoreZcutoff
  keep <- rowSums(openscoreZ>=openscoreZcutoff) > 0
  openscoreP <- openscoreP[keep, , drop=FALSE]
  bindingscore <- bindingscore[keep, , drop=FALSE]
  bindingSites <- bindingSites[keep]

  ## following gaussion distribution, should we use Skewness to value it?
  adj.bindingscore <- bindingscore * openscoreP
  ## remove zero sample variances
  keep <- apply(adj.bindingscore, 1, unique)
  keep <- lengths(keep)>1
  adj.bindingscore <- adj.bindingscore[keep, , drop=FALSE]
  bindingSites <- bindingSites[keep]

  design <- cbind(CTL=1, EXPvsCTL=bams$group=="Exp")
  fit <- lmFit(adj.bindingscore, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef="EXPvsCTL", adjust.method="BH",
                 sort.by = "none", number=nrow(fit))

  ## sort the results by t value
  mcols(bindingSites) <- cbind(mcols(bindingSites), tt)
  ## remove the log2 foldchange smaller than bindingScoreLog2FCcutoff
  ## and pvalue greater than bindingScorePvalCutoff
  ## otherwise dataset is too large
  bindingSites <-
    bindingSites[bindingSites$P.Value<=bindingScorePvalCutoff &
                   abs(bindingSites$logFC)>=bindingScoreLog2FCcutoff]
  stopifnot(
    "bindingScorePvalCutoff & bindingScoreLog2FCcutoff is too stringency."=
      length(bindingSites)>0)
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

  p <- vapply(seq_along(ESm),
              function(i) {
                -2*sum((-1)^k*exp(-2*k^2*ESm[i]^2*n[i]))
                },
              FUN.VALUE = 0.0)
  ## expect ES
  k <- seq(1, 10000)
  EES <-
    vapply(seq_along(n),
           function(i){
             8*sum((-1)^(k+1)*(exp(-2*k^2*n[i])/4 -
                                 (sqrt(2*pi)*erf(k*sqrt(2*n[i])))/
                                 (16*k*sqrt(n[i]))))
             },
           FUN.VALUE = 0.0)
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



