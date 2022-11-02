#!/usr/bin/env Rscript
# 2022 Jianhong Ou
# report issue at https://github.com/jianhong/ATACseqTFEA/issue
suppressPackageStartupMessages(library(optparse))
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-m", "--motifpath"), type="character", default='PWMatrixList.rds',
              help="The file path 'PWMatrixList' or 'PFMatrixList' saved in R object (read by function readRDS) [default %default]. The seaching path including the extdata of ATACseqTFEA package"),
  make_option(c("-g", "--bsgenome"), type="character", default="BSgenome.Hsapiens.UCSC.hg38",
              help = "BSgenome package name [default \"%default\"]"),
  make_option(c("-t", "--txdb"), type="character", default="TxDb.Hsapiens.UCSC.hg38.knownGene",
              help="TxDb package name [default %default]"),
  make_option(c("-s", "--seqlevels"), type="character", default='paste0("chr", 1:25)',
              help="seqlevels for binding sites [default %default]"),
  make_option(c("-p", "--peakfiles"), type="character",
              help="Peak files for ATAC-seq (required)."),
  make_option(c("-f", "--format"), type="character", default='narrowPeak',
              help="Format of ATAC-seq peak files [default %default]"),
  make_option(c("-e", "--experiment"), type="character",
              help="Bam file of experiments/treatments"),
  make_option(c("-c", "--control"), type="character",
              help="Bam file of controls"),
  make_option(c("-o", "--outfolder"), type="character", default=".",
              help="Output folder [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$peakfiles)){
  stop("--peakfiles is required.")
}

console_log <- function(...){
  if ( opt$verbose ) {
    message(...)
  }
}

bamExp <- eval(parse(text=opt$experiment))
bamCtl <- eval(parse(text=opt$control))

for(x in c(bamExp, bamCtl)){
  if(!file.exists(x)){
    stop(x, ' is not available.')
  }
}

if(file.exists(system.file("extdata", opt$motifpath, package = "ATACseqTFEA"))){
  motifpath <- system.file("extdata", opt$motifpath, package = "ATACseqTFEA")
}else{
  motifpath <- opt$motifpath
}

suppressPackageStartupMessages(library(ATACseqTFEA))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(ATACseqQC))
if(!require(package=opt$bsgenome, character.only = TRUE)){
  stop("can not load bsgenome.")
}
if(!require(package=opt$txdb, character.only = TRUE)){
  stop("can not load txdb.")
}

console_log("loading motifs from ", motifpath)
motifs <- readRDS(motifpath)

genome <- get(opt$bsgenome)

txdb <- get(opt$txdb)
genes <- genes(txdb)

## only search the ranges for peaks.
seqlev <- eval(parse(text=opt$seqlevels))
peaks <- lapply(opt$peakfiles, import, format=opt$format)
peaks <- reduce(unlist(GRangesList(peaks)))
console_log("prepare binding sites.")
mts <- prepareBindingSites(motifs, genome, seqlev,
                           grange=peaks,
                           p.cutoff = 1e-05)
console_log("doing TFEA")
res <- TFEA(bamExp, bamCtl, positive=0, negative=0, bindingSites=mts,
            filter="proximalRegion>1")

dir.create(opts$outfolder, recursive=TRUE)
console_log("save the results as a csv file")
## export the results into a csv file
write.csv(res$resultsTable, file.path(opts$outfolder, "enrichment.csv"),
          row.names=FALSE)
## export the results into a rds file
saveRDS(res, file.path(opts$outfolder, "res.rds"))

console_log("plot Enrichment scores")
## get all the enrichment score plots
dir.create(file.path(opts$outfolder, "ESplots"), recursive=TRUE)
g <- plotES(res, outfolder=file.path(opts$outfolder, "ESplots"))
## volcanoplot
pdf(file.path(opts$outfolder, "volcanoplot.pdf"))
ESvolcanoplot(TFEAresults=res)
dev.off()

console_log("plot footprints")
## plot the footprint
outfolder <- file.path(opts$outfolder, "footprint")
dir.create(outfolder)
for(TF in res$resultsTable$TF[res$resultsTable$adjPval<0.05]){
  bs <- getBindingSites(res)
  bs <- bs[vapply(bs$motif,
                  FUN=function(.ele) TF %in% .ele,
                  FUN.VALUE = logical(1L))]
  bs$score <- mapply(bs$score, bs$motif, FUN = function(score, motif){
    max(score[motif==TF])
  })
  pdf(file.path(outfolder, paste(TF, ".pdf")))
  factorFootprints(c(bamCtl, bamExp),
                   pfm = as.matrix(motifs[[TF]]),
                   bindingSites = bs,
                   seqlev = seqlev, genome = Drerio,
                   upstream = 100, downstream = 100,
                   group = rep(c("WT", "KD"),
                               c(length(bamCtl), length(bamExp))))
  dev.off()
}

