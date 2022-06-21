#' Prepare the genomic ranges for proximal and distal regions for counting
#' @description Create multiple GRanges objects for downstream counting.
#' The GRanges objects including
#' bindingSitesWithGap: bindingSites with gaps and in both ends,
#' bindingSitesWithProximal: bindingSites with gaps and proximal region
#' in both ends,
#' bindingSitesWithProximalAndGap: bindingSites with gaps, and then proximal
#' and gaps in both ends,
#' and bindingSitesWithDistal: bindingSites with gaps, proximal, gaps and
#' distal regions.
#' @param bindingSites A object of
#' \link[GenomicRanges:GRanges-class]{GenomicRanges} indicates
#' candidate binding sites. The \link{prepareBindingSites} function
#' is a helper function to generate the binding sites.
#' Users can also use other software for example fimo to generate the list.
#' @param proximal,distal numeric(1) or integer(1).
#'        bases for open region from binding sites (proximal) and
#'        extended region for background (distal)
#'        of the binding region for aggregate ATAC-seq footprint.
#' @param gap numeric(1) or integer(1). bases for gaps among binding sites,
#'            proximal, and distal. default is 10L.
#' @return an \link[GenomicRanges:GRangesList-class]{GRangesList} object with
#' elements bindingSitesWithGap, bindingSitesWithProximal,
#' bindingSitesWithProximalAndGap, and bindingSitesWithDistal
#' for \link{count5ends}
#' @importFrom BiocGenerics start end width start<- end<-
#' @export
#' @author Jianhong Ou
#' @examples
#' bsl <- system.file("extdata", "bindingSites.rds",
#'                    package="ATACseqTFEA")
#' bindingSites <- readRDS(bsl)
#' bs <- expandBindingSites(bindingSites)
#' length(bs)
#' names(bs)
#' lengths(bs)
expandBindingSites <- function(bindingSites,
                               proximal=40L,
                               distal=proximal,
                               gap=10L){
  stopifnot("bindingSites must be an GRanges object"=
              is(bindingSites, "GRanges"))
  proximal <- checkInteger(proximal)
  distal <- checkInteger(distal)
  gap <- checkInteger(gap)
  bindingSites.with.proximal <-
    bindingSites.with.distal <-
    bindingSites
  bindingSites.with.proximal.gap <-
    bindingSites.with.distal.gap <-
    bindingSites
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
  return(GRangesList(
    bindingSitesWithGap=bindingSites.with.proximal.gap,
    bindingSitesWithProximal=bindingSites.with.proximal,
    bindingSitesWithProximalAndGap=bindingSites.with.distal.gap,
    bindingSitesWithDistal=bindingSites.with.distal
  ))
}
