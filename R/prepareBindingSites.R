#' Prepare binding site for TFEA
#' @description Prepare binding sites by ginven position weight matrix and
#'  genome.
#' @param pwms either \code{\link[TFBSTools]{PFMatrix}},
#' \code{\link[TFBSTools]{PFMatrixList}}, \code{\link[TFBSTools]{PWMatrix}},
#' \code{\link[TFBSTools]{PWMatrixList}}
#' @param genome \code{\link[BSgenome:BSgenome-class]{BSgenome}} object.
#' @param seqlev A character vector. Sequence levels to be searched.
#' @param p.cutoff p-value cutoff for returning motifs; default is 5e-05
#' @param w parameter controlling size of window for filtration; default is 7
#' @param grange GRanges for motif search. If it is set, function will only
#' search the binding site within the grange.
#' @param maximalBindingWidth A numeric vector(length=1).
#' Maximal binding site width. Default is 40.
#' @param mergeBindingSitesByPercentage A numeric vector (length=1).
#' The percentage of overlapping region of binding sites to merge as one
#' binding site.
#' @return A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} with
#' all the positions of matches.
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqinfo seqlevels `seqlevels<-` Seqinfo seqlengths
#' @importFrom motifmatchr matchMotifs
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics `%in%`
#' @export
#' @author Jianhong Ou
#' @examples
#' library(TFBSTools)
#' motifs <- readRDS(system.file("extdata", "PWMatrixList.rds",
#'                               package="ATACseqTFEA"))
#' library(BSgenome.Drerio.UCSC.danRer10)
#' seqlev <- "chr1" #paste0("chr", 1:25)
#' mts <- prepareBindingSites(motifs, Drerio, seqlev,
#'                            grange=GRanges("chr1",
#'                                           IRanges(5000, 100000)))
prepareBindingSites <- function(pwms, genome, seqlev=seqlevels(genome),
                                p.cutoff = 5e-05, w = 7, grange,
                                maximalBindingWidth = 40L,
                                mergeBindingSitesByPercentage = 0.8){
  stopifnot("genome must be a BSgenome object."=is(genome, "BSgenome"))
  if(!is.na(maximalBindingWidth[1])){
    maximalBindingWidth <- maximalBindingWidth[1]
    stopifnot(is.numeric(maximalBindingWidth))
  }
  si <- seqinfo(genome)
  stopifnot("All seqlev must be in genome."=all(seqlev %in% names(si)))
  if(missing(grange)){
    p <- as(si[seqlev], "GRanges")
  }else{
    stopifnot("grange must be a GRanges Object"=is(grange, "GRanges"))
    p <- grange
  }
  motif_pos <- matchMotifs(pwms=pwms, subject=p,
                           genome = genome, out = "positions",
                           p.cutoff = p.cutoff, w = w)
  mts.unlist <- unlist(motif_pos, use.names = FALSE)
  mts.unlist$score <- NULL
  mts.unlist$motif <- rep(names(motif_pos), lengths(motif_pos))
  seqlev <- intersect(seqlevels(mts.unlist), seqlev)
  mts.unlist <- mts.unlist[seqnames(mts.unlist) %in% seqlev]
  seqlevels(mts.unlist) <-
    seqlevels(mts.unlist)[seqlevels(mts.unlist) %in% seqlev]
  seqinfo(mts.unlist) <- si[seqlev]
  rm(motif_pos)
  ## reduce is not good.
  ## try a way to reduce by percentage of overlaps
  reduceByPercentage <- function(query, percentage){
    ol <- findOverlaps(query = query, drop.self=TRUE, drop.redundant=TRUE)
    if(length(ol)==0) return(query)
    qh <- query[queryHits(ol)]
    sh <- query[subjectHits(ol)]
    oh <- GRanges(seqnames(qh),
                  IRanges(start=ifelse(start(qh)>start(sh),
                                       start(qh), start(sh)),
                          end=ifelse(end(qh)<end(sh),
                                     end(qh), end(sh))),
                  strand = strand(qh))
    wqh <- width(qh)
    wsh <- width(sh)
    woh <- width(oh)
    pqh <- woh/wqh
    psh <- woh/wsh
    keep <- pqh>=percentage | psh>=percentage
    if(sum(keep)==0) return(query)
    ol <- ol[keep]

    li <- as.data.frame(ol)
    qHid <- unique(queryHits(ol))
    li <- rbind(li, data.frame(queryHits=qHid, subjectHits=qHid))
    li$queryHits <- formatC(li$queryHits,
                            width = nchar(as.character(max(li$queryHits))),
                            flag = "0")
    q2 <- split(query[li$subjectHits], li$queryHits)
    q2 <- reduce(GRangesList(q2))
    q2 <- unlist(q2)
    q2$motif <- split(query$motif[li$subjectHits], li$queryHits)[names(q2)]

    c(query[-unique(c(queryHits(ol), subjectHits(ol)))], q2)
  }
  mts.unlist <- split(mts.unlist, seqnames(mts.unlist)) ## memory
  mts.reduce <- list()
  for(i in seq_along(mts.unlist)){
    mts.reduce[[i]] <-
      reduceByPercentage(mts.unlist[[i]],
                         percentage = mergeBindingSitesByPercentage)
  }
  rm(mts.unlist)
  mts.reduce <- mts.reduce[lengths(mts.reduce)>0]
  mts.reduce <- unlist(GRangesList(mts.reduce))

  ## remove the binding site that greater than maximalBindingWidth.
  if(!is.na(maximalBindingWidth)){
    mts.reduce <- mts.reduce[width(mts.reduce)<=maximalBindingWidth]
  }

#  names(mts.reduce) <-
#    paste0("p", formatC(seq.int(length(mts.reduce)),
#                          width = nchar(as.character(length(mts.reduce))),
#                          flag = "0"))
  return(mts.reduce)
}
