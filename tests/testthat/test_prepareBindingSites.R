test_that("prepareBindingsites works not correct", {
  motifs <- readRDS(system.file("extdata", "PWMatrixList.rds",
                                package="ATACseqTFEA"))
  seqlev <- "chr1"
  mts <- prepareBindingSites(motifs, Drerio, seqlev="chr1",
                             grange=GRanges("chr1",
                                            IRanges(5000, 100000)),
                             p.cutoff = 5e-05)
  expect_false(any(duplicated(mts)))
  expect_true(all(width(mts)<40L))
  ## test ignore.strand = FALSE
  mts <- prepareBindingSites(motifs, Drerio, seqlev="chr1",
                             grange=GRanges("chr1",
                                            IRanges(5000, 100000)),
                             p.cutoff = 1e-05,
                             ignore.strand = FALSE)
})
