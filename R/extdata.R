#' Data in extdata
#' @description The list of data saved in extdata folder.
#' @aliases cluster_PWMs
#' @aliases PWMatrixList
#' @name extdata
#' @rdname extdata
#' @docType data
#' @details
#' The `PWMatrixList` is a collection of jasper2018, jolma2013 and
#' cisbp_1.02 from package motifDB (v 1.28.0) and merged by distance smaller
#' than 1e-9 calculated by MotIV::motifDistances function (v 1.42.0).
#' The merged motifs were exported by motifStack (v 1.30.0).
#'
#' The `cluster_PWMs` is a list of non-redundant TF motifs downloaded from
#' [DeepSTARR](https://github.com/bernardo-de-almeida/motif-clustering).
#' There are 6502 motifs in the data set.
#'
#' The `best_curated_Human` is a list of TF motifs downloaded from
#' [TFEA github](https://github.com/Dowell-Lab/TFEA).
#' There are 1279 human motifs in the data set.
#' @examples
#' motifs <- readRDS(system.file("extdata", "PWMatrixList.rds",
#'                   package="ATACseqTFEA"))
#' motifs2 <- readRDS(system.file("extdata", "cluster_PWMs.rds",
#'                    package="ATACseqTFEA"))
#' motifs3 <- readRDS(system.file("extdata", "best_curated_Human.rds",
#'                    package="ATACseqTFEA"))
NULL
