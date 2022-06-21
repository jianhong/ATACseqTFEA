#' Plot enrichment score for one transcription factor
#' @description Plot GSEA style enrichment score curve.
#' @param TFEAresults A \link{TFEAresults} object. Output of \link{TFEA}.
#' @param TF A character vector. The transcription factor names.
#' @param outfolder character(1). Output file path.
#' @param xlab,ylab character string giving label for x-axis/y-axis.
#' @param resolution integer(1). The number of bars plotted in the bottom of
#' figure to show the density of occurrence of events.
#' @param ... parameter passed to pdf.
#' @return NULL if outfolder is set or ggplot object.
#' @importFrom ggplot2 ggplot aes_string geom_line geom_rug xlab ylab
#' theme_classic geom_hline ggtitle
#' @importFrom dplyr sample_n
#' @importFrom grDevices dev.off pdf
#' @export
#' @examples
#' res <- system.file("extdata", "res.rds", package="ATACseqTFEA")
#' res <- readRDS(res)
#' g <- plotES(res, TF="KLF9", outfolder=NA)
#' print(g)
plotES <- function(TFEAresults,
                   TF,
                   outfolder=".",
                   xlab="rank",
                   ylab="Enrichment",
                   resolution=500L,
                   ...){
  stopifnot("TFEAresults must be output of TFEA function"=
              is(TFEAresults, "TFEAresults"))
  ES <- t(getEnrichmentScore(TFEAresults))
  ESplot <- function(ES, i, xlab, ylab, resolution){
    dat <- data.frame(cbind(x=seq.int(nrow(ES)), y=ES[, i]))
    p <- ggplot(dat, aes_string(x="x", y="y")) +
      geom_line() +
      geom_rug(data=
                 sample_n(subset(dat, dat$x %in%
                                   getMotifID(TFEAresults)[[i]]),
                          size=min(resolution,
                                   nrow(subset(dat, dat$x %in%
                                                 getMotifID(TFEAresults)[[i]]
                                               )))),
               sides = "b", position = "jitter") +
      xlab(xlab) + ylab(ylab) + theme_classic() +
      geom_hline(yintercept = 0) + ggtitle(i)
    p
  }
  if(missing(TF)){
    if(!is.na(outfolder)&&!is.null(outfolder)){
      if(!file.exists(outfolder)){
        dir.create(outfolder)
      }
      null <- lapply(colnames(ES), function(i){
        pdf(file.path(outfolder, paste0(make.names(i), ".pdf")), ...)
        print(ESplot(ES, i, xlab, ylab, resolution))
        dev.off()
      })
    }else{
      stop("outfolder can not be NA if TF is not set.")
    }
  }else{
    if(!is.na(outfolder)&&!is.null(outfolder)){
      if(!file.exists(outfolder)){
        dir.create(outfolder)
      }
      null <- lapply(TF, function(i){
        if(i %in% colnames(ES)){
          pdf(file.path(outfolder, paste0(make.names(i), ".pdf")), ...)
          print(ESplot(ES, i, xlab, ylab, resolution))
          dev.off()
        }else{
          warning(i, "is not a valid TF name.")
        }
      })
    }else{
      if(length(TF)==1){
        i <- TF
        if(i %in% colnames(ES)){
          p <- ESplot(ES, i, xlab, ylab, resolution)
          p
        }else{
          warning(i, "is not a valid TF name.")
        }
      }else{
        stop("If multiple TFs are given, outfolder can not be NA.")
      }
    }
  }
}
