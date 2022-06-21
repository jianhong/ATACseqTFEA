#' Plot enrichment score for one transcription factor
#' @description Plot GSEA style enrichment score curve.
#' @param TFEAresults A \link{TFEAresults} object. Output of \link{TFEA}.
#' @param xlab,ylab character string giving label for x-axis/y-axis.
#' @param TFnameToShow Transcription factor names to be drawn.
#' @param significantCutoff Cutoff value for significant.
#' @param col Color sets for the points.
#' @param ... parameter passed to pdf.
#' @return ggplot object.
#' @importFrom ggplot2 ggplot aes_string theme_bw labs geom_point
#' scale_color_manual unit
#' @importFrom ggrepel geom_text_repel
#' @export
#' @examples
#' res <- system.file("extdata", "res.rds", package="ATACseqTFEA")
#' res <- readRDS(res)
#' ESvolcanoplot(TFEAresults=res)
ESvolcanoplot <- function(TFEAresults,
                          xlab="Enrichment Score",
                          ylab="-log10(p value)",
                          TFnameToShow = 20,
                          significantCutoff = 0.05,
                          col = c("red", "blue", "gray"),
                          ...){
  stopifnot("TFEAresults must be output of TFEA function"=
              is(TFEAresults, "TFEAresults"))
  res <- as(TFEAresults, "data.frame")
  res$qvalue <- -log10(res$p_value)
  res <- res[!is.na(res$enrichmentScore), , drop=FALSE]
  res <- res[!is.na(res$qvalue), , drop=FALSE]
  stopifnot("No valid results."=nrow(res)>0)
  res$significant <- res$adjPval < significantCutoff
  res$Color <- col[3]
  res$Color[res$significant & res$enrichmentScore>0] <- col[1]
  res$Color[res$significant & res$enrichmentScore<0] <- col[2]
  res$symbol <- as.character(res$TF)
  if(!is.numeric(TFnameToShow)[1]){
    res$symbol[!res$symbol %in% TFnameToShow] <- NA
  }else{
    res$symbol[res$Color=="gray"] <- NA
    res$symbol[res$symbol=="-"] <- NA
    if(sum(!is.na(res$symbol))>TFnameToShow[1]){
      res <- res[order(abs(res$enrichmentScore),
                                 decreasing = TRUE), ]
      res$counts <- 0
      res$counts[!is.na(res$symbol)] <- seq.int(sum(!is.na(res$symbol)))
      res$symbol[res$counts>TFnameToShow[1]] <- NA
    }
  }
  g <- ggplot(data=res, aes_string(x="enrichmentScore", y="qvalue")) +
    theme_bw() +
    labs(x=xlab, y=ylab) +
    geom_point(aes_string(color="Color")) +
    geom_text_repel(data=subset(res, !is.na(res$symbol)),
                    aes_string(label="symbol"),
                    box.padding = unit(0.45, "lines")) +
    scale_color_manual(values = levels(factor(res$Color)), guide="none")
  g
}
