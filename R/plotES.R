#' Plot enrichment score for one transcription factor
#' @description Plot GSEA style enrichment score curve.
#' @param TFEAresults A \link{TFEAresults} object. Output of \link{TFEA}.
#' @param TF A character vector. The transcription factor names.
#' @param outfolder character(1). Output file path.
#' @param xlab,ylab character string giving label for x-axis/y-axis.
#' @param ... parameter passed to pdf.
#' @return NULL if outfolder is set or ggplot object.
#' @importFrom ggplot2 ggplot aes_string geom_line geom_rug xlab ylab
#' theme_classic geom_hline
#' @importFrom grDevices dev.off pdf
#' @export
#' @examples
#' res <- system.file("extdata", "res.rds", package="ATACseqTFEA")
#' res <- readRDS(res)
#' g <- plotES(res, TF="KLF9", outfolder=NA)
plotES <- function(TFEAresults, TF, outfolder=".",
                   xlab="rank", ylab="Enrichment", ...){
  stopifnot("TFEAresults must be output of TFEA function"=
              is(TFEAresults, "TFEAresults"))
  ES <- t(TFEAresults@enrichmentScore)
  ESplot <- function(ES, i, xlab, ylab){
    dat <- data.frame(cbind(x=seq.int(nrow(ES)), ES))
    p <- ggplot(dat, aes_string(x="x", y=i)) +
      geom_line() +
      geom_rug(data=subset(dat, dat$x %in% TFEAresults@motifID[[i]]),
               sides = "b", position = "jitter") +
      xlab(xlab) + ylab(ylab) + theme_classic() +
      geom_hline(yintercept = 0)
  }
  if(missing(TF)){
    if(!is.na(outfolder)&&!is.null(outfolder)){
      if(!file.exists(outfolder)){
        dir.create(outfolder)
      }
      for(i in colnames(ES)){
        pdf(file.path(outfolder, paste0(make.names(i), ".pdf")), ...)
        ESplot(ES, i, xlab, ylab)
        dev.off()
      }
    }else{
      stop("outfolder can not be NA if TF is not set.")
    }
  }else{
    if(!is.na(outfolder)&&!is.null(outfolder)){
      if(!file.exists(outfolder)){
        dir.create(outfolder)
      }
      for(i in TF){
        if(i %in% colnames(ES)){
          pdf(file.path(outfolder, paste0(make.names(i), ".pdf")), ...)
          ESplot(ES, i, xlab, ylab)
          dev.off()
        }else{
          warning(i, "is not a valid TF name.")
        }
      }
    }else{
      if(length(TF)==1){
        i <- TF
        if(i %in% colnames(ES)){
          ESplot(ES, i, xlab, ylab)
        }else{
          warning(i, "is not a valid TF name.")
        }
      }else{
        stop("If multiple TFs are given, outfolder can not be NA.")
      }
    }
  }
}
