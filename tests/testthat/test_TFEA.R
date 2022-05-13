test_that("TFEA works not correct", {
  bamExp <- system.file("extdata",
                        c("KD.shift.rep1.bam",
                          "KD.shift.rep2.bam"),
                        package="ATACseqTFEA")
  bamCtl <- system.file("extdata",
                        c("WT.shift.rep1.bam",
                          "WT.shift.rep2.bam"),
                        package="ATACseqTFEA")
  bsl <- system.file("extdata", "bindingSites.rds",
                     package="ATACseqTFEA")
  bindingSites <- readRDS(bsl)
  res <- TFEA(bamExp, bamCtl, bindingSites=bindingSites)
  es <- getEnrichmentScore(res)
  es_end <- es[, ncol(es)]
  expect_true(all(es_end==0))
})
