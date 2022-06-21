#' Differential binding analysis
#' @description Use limma to do differential binding analysis for binding
#'  scores.
#' @param se An \link[SummarizedExperiment:RangedSummarizedExperiment-class]{RangedSummarizedExperiment}
#'  object. Outputs of \link{getWeightedBindingScore}.
#' @param design Design table for \link[limma:lmFit]{lmFit}.
#' @param ... Parameters can be used by \link[limma:lmFit]{lmFit}.
#' @param coef column number or column name specifying which coefficient or
#' contrast of the linear model is of interest.
#' See \link[limma:topTable]{topTable}.
#' @return A RangedSummarizedExperiment object with the dataframe returned by
#'  \link[limma:topTable]{topTable} as appendence of the origin rowData.
#' @importFrom limma lmFit eBayes topTable
#' @importFrom SummarizedExperiment assays rowData rowData<-
#' @export
#' @author Jianhong Ou
#' @examples
#' library(SummarizedExperiment)
#' set.seed(1)
#' sigma2 <- 0.05 / rchisq(100, df=10) * 10
#' y <- matrix(rnorm(100*6,sd=sqrt(sigma2)),100,6)
#' design <- cbind(Intercept=1,Group=c(0,0,0,1,1,1))
#' y[1,4:6] <- y[1,4:6] + 1
#' se <- SummarizedExperiment(assays=list(counts=y))
#' DBscore(se, design, coef=1)
DBscore <- function(se, design, coef, ...){
  stopifnot(is(se, "SummarizedExperiment"))
  adj.bindingscore <- assays(se)[[1]]
  fit <- lmFit(adj.bindingscore, design, ...)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef=coef, adjust.method="BH",
                 sort.by = "none", number=nrow(fit))
  rowData(se) <- cbind(rowData(se), adj.bindingscore, tt)
  se
}
