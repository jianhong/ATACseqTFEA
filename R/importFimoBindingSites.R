#' Prepare binding site by fimo results
#' @description Prepare binding sites by given fimo gff files
#' @param fimoGFFfiles Filenames of gff files of fimo output.
#' @param maximalBindingWidth A numeric vector(length=1).
#' Maximal binding site width. Default is 40.
#' @param mergeBindingSitesByPercentage A numeric vector (length=1).
#' The percentage of overlapping region of binding sites to merge as one
#' binding site.
#' @param ignore.strand When set to TRUE, the strand information is ignored in
#' the calculations.
#' @param ... Parameter to be passed to \link[rtracklayer:GFFFile-class]{import.gff}
#' @return A \code{\link[GenomicRanges:GRanges-class]{GenomicRanges}} with
#' all the positions of matches.
#' @importFrom GenomeInfoDb seqinfo seqlevels `seqlevels<-` Seqinfo seqlengths
#' seqnames seqinfo<-
#' @importFrom rtracklayer import
#' @importFrom IRanges IRanges reduce findOverlaps pintersect countOverlaps
#' distance promoters
#' @importFrom S4Vectors queryHits subjectHits split mcols mcols<-
#' @importFrom BiocGenerics `%in%`
#' @importFrom GenomicRanges GRangesList GRanges
#' @export
#' @author Jianhong Ou
#' @examples
#' extdata <- system.file('extdata', package='ATACseqTFEA')
#' fimoGFFfiles <- dir(extdata, 'fimo.*.gff', full.names=TRUE)
#' mts <- importFimoBindingSites(fimoGFFfiles)
importFimoBindingSites <- function(fimoGFFfiles,
                                   maximalBindingWidth = 40L,
                                   mergeBindingSitesByPercentage = 0.8,
                                   ignore.strand = TRUE,
                                   ...){
  stopifnot(is.character(fimoGFFfiles))
  if(!is.na(maximalBindingWidth[1])){
    maximalBindingWidth <- maximalBindingWidth[1]
    stopifnot(is.numeric(maximalBindingWidth))
  }
  motif_pos <- lapply(fimoGFFfiles, import, ...)

  mts.unlist <- unlist(GRangesList(motif_pos), use.names = FALSE)
  mts.unlist$motif <- mts.unlist$Name
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
