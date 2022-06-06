#' Prepare binding site for TFEA
#' @description Prepare binding sites by given position weight matrix and
#'  genome.
#' @param pwms either \code{\link[TFBSTools]{PFMatrix}},
#' \code{\link[TFBSTools]{PFMatrixList}}, \code{\link[TFBSTools]{PWMatrix}},
#' \code{\link[TFBSTools]{PWMatrixList}}
#' @param genome \code{\link[BSgenome:BSgenome-class]{BSgenome}} object.
#' @param seqlev A character vector. Sequence levels to be searched.
#' @param p.cutoff p-value cutoff for returning motifs; default is 1e-05
#' @param w parameter controlling size of window for filtration; default is 7
#' @param grange GRanges for motif search. If it is set, function will only
#' search the binding site within the grange. Usually a peak list should be
#' supplied.
#' @param maximalBindingWidth A numeric vector(length=1).
#' Maximal binding site width. Default is 40.
#' @param mergeBindingSitesByPercentage A numeric vector (length=1).
#' The percentage of overlapping region of binding sites to merge as one
#' binding site.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#' the calculations.
#' @return A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} with
#' all the positions of matches.
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqinfo seqlevels `seqlevels<-` Seqinfo seqlengths
#' @importFrom motifmatchr matchMotifs
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom BiocGenerics `%in%`
#' @importFrom TFBSTools ID
#' @importFrom utils combn
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
                                p.cutoff = 1e-05, w = 7, grange,
                                maximalBindingWidth = 40L,
                                mergeBindingSitesByPercentage = 0.8,
                                ignore.strand = TRUE){
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
  stopifnot((inherits(pwms, c("PFMatrixList", "PWMatrixList"))))
  if(length(names(pwms))==0){
    n <- vapply(pwms, FUN=ID, FUN.VALUE = character(1))
    stopifnot('Can not identify the name of pwms'=length(n)==length(pwms))
    names(pwms) <- n
  }
  motif_pos <- matchMotifs(pwms=pwms, subject=p,
                           genome = genome, out = "positions",
                           p.cutoff = p.cutoff, w = w)
  mts.unlist <- unlist(motif_pos, use.names = FALSE)
  mts.unlist$motif <- rep(names(motif_pos), lengths(motif_pos))
  seqlev <- intersect(seqlevels(mts.unlist), seqlev)
  mts.unlist <- mts.unlist[seqnames(mts.unlist) %in% seqlev]
  seqlevels(mts.unlist) <-
    seqlevels(mts.unlist)[seqlevels(mts.unlist) %in% seqlev]
  seqinfo(mts.unlist) <- si[seqlev]
  rm(motif_pos)

  mts.unlist <- split(mts.unlist, seqnames(mts.unlist)) ## memory
  mts.reduce <- lapply(mts.unlist, reduceByPercentage,
                       percentage = mergeBindingSitesByPercentage,
                       ignore.strand = ignore.strand)
  rm(mts.unlist)
  mts.reduce <- mts.reduce[lengths(mts.reduce)>0]
  mts.reduce <- unlist(GRangesList(mts.reduce))

  ## remove the binding site that greater than maximalBindingWidth.
  if(!is.na(maximalBindingWidth)){
    mts.reduce <- mts.reduce[width(mts.reduce)<=maximalBindingWidth]
  }

  return(mts.reduce)
}


findOverlaps1 <- function(query, percentage, ignore.strand, ...){
  stopifnot(percentage[1]>0 && percentage[1]<1)
  hits <- findOverlaps(query = query,
                       minoverlap = 0L,
                       maxgap = -1L,
                       ...)
  overlaps <- pintersect(query[queryHits(hits)],
                         query[subjectHits(hits)],
                         ignore.strand = ignore.strand)
  mcols(overlaps) <- NULL
  percentOverlap0 <- width(overlaps)/width(query[queryHits(hits)])
  percentOverlap1 <- width(overlaps)/width(query[subjectHits(hits)])
  percentOverlap <- ifelse(
    percentOverlap0>percentOverlap1,
    percentOverlap0, percentOverlap1
  )
  keep <- percentOverlap>=percentage
  list(hits=hits[keep], overlaps=overlaps[keep])
}

reduceList <- function(query_new){
  l <- query_new$qid
  if(any(duplicated(l))){ ## remove the duplicated items
    query_new <- split(query_new,
               vapply(l, FUN=paste, collapse=",",FUN.VALUE = character(1L)))
    query_new <- reduce(query_new)
    query_new <- unlist(query_new)
    query_new$qid <- strsplit(names(query_new), split=",")
    names(query_new) <- NULL
    l <- query_new$qid
  }
  ## remove the subset items
  ol <- findOverlaps(query_new,
                     drop.self=TRUE,
                     drop.redundant=FALSE,
                     minoverlap = 1L)
  ol <- ol[lengths(l[queryHits(ol)])<lengths(l[subjectHits(ol)])]
  sub <- mapply(setdiff, l[queryHits(ol)], l[subjectHits(ol)],
                SIMPLIFY = FALSE)
  sub <- lengths(sub)==0
  sub <- split(sub, queryHits(ol))
  sub <- vapply(sub, FUN = any, FUN.VALUE = logical(1L))
  sub <- sub[as.character(seq_along(query_new))]
  sub[is.na(sub)] <- FALSE
  query_new[!sub]
}

#' Reduce by percentage of overlaps of GRanges object
#' @description Merge the ranges by percentage of overlaps to avoid
#' broad ranges of continues ranges overlapped with limit bases.
#' @param query An object of GRanges
#' @param percentage A numeric vector (length=1).
#' The percentage of overlapping region of binding sites to merge as one
#' range.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#' the calculations.
#' @param colnToKeep The metadata colnums should be kept for reduced GRanges
#' @return An object of GRanges.
#' @export
#' @import GenomicRanges
#' @importFrom S4Vectors SimpleList
#' @examples
#' gr <- GRanges("chr1", IRanges(c(1, 5, 10), width=c(10, 5, 2)))
#' reduceByPercentage(gr, 0.5, colnToKeep=NULL)

reduceByPercentage <- function(query, percentage, ignore.strand=TRUE,
                               colnToKeep=c("score", "motif")){
  stopifnot(is(query, "GRanges"))
  stopifnot(percentage[1]>0 && percentage[1]<1)
  stopifnot(is.logical(ignore.strand))
  stopifnot(length(ignore.strand)==1)
  if(length(colnToKeep)) {
    stopifnot(all(colnToKeep %in% colnames(mcols(query))))
  }
  query_new <- query
  mcols(query_new) <- NULL
  query_new$qid <- seq_along(query)
  while(TRUE){
    len1 <- length(query_new)
    ol <- findOverlaps1(query_new,
                        percentage = percentage,
                        ignore.strand=ignore.strand,
                        select="all",
                        drop.self=FALSE,
                        drop.redundant=FALSE)
    if(length(ol$hits)==length(query_new)) break
    overlaps_key <- ol$overlaps
    ol <- ol$hits
    stopifnot(all(subjectHits(ol) %in% queryHits(ol)))
    stopifnot(all(queryHits(ol) %in% subjectHits(ol)))
    stopifnot(sum(queryHits(ol)==subjectHits(ol))==length(query_new))
    query_new <- query_new[subjectHits(ol)]
    q2_qid <- split(query_new$qid,
                            queryHits(ol))
    q2_qid <- lapply(q2_qid, unlist, use.names = FALSE)
    q2_qid <- lapply(q2_qid, unique)
    stopifnot(all(seq_along(query) %in% unlist(q2_qid)))
    query_new <- split(overlaps_key, queryHits(ol))
    query_new <- reduce(query_new, ignore.strand=ignore.strand)
    query_new <- unlist(query_new)
    query_new$qid <- q2_qid
    query_new <- sort(query_new)
    q2_qid <- paste(
      seqnames(query_new),
      start(query_new),
      vapply(query_new$qid,
             FUN=function(.ele) paste(.ele, collapse=";"),
             FUN.VALUE = character(1)))
    query_new <- query_new[!duplicated(q2_qid)]
    #check qid correlation
    query_new <- reduceList(query_new)
    len2 <- length(query_new)
    if(len2==len1) break
  }
  stopifnot('The reduce of binding sites does not work correctly.'=
              all(seq_along(query) %in% unlist(query_new$qid)))
  query_new <- sort(query_new)
  ## reduce binding range
  q2 <- lapply(query_new$qid, function(.ele){
    query[as.numeric(.ele)]
  })
  query_new$qid <- NULL
  if(length(colnToKeep)){
    for(i in colnToKeep){
      mcols(query_new)[i] <-
        SimpleList(lapply(q2, function(.ele) mcols(.ele)[, i]))
    }
  }
  query_new
}

