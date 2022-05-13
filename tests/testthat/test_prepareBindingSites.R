test_that("prepareBindingsites works not correct", {
  motifs <- readRDS(system.file("extdata", "PWMatrixList.rds",
                                package="ATACseqTFEA"))
  seqlev <- "chr1"
  mts <- prepareBindingSites(motifs, Drerio, seqlev="chr1",
                             grange=GRanges("chr1",
                                            IRanges(5000, 100000)),
                             p.cutoff = 5e-05)
  res <- readRDS((system.file("extdata", "bindingSites.rds",
                              package="ATACseqTFEA")))
  expect_true(all(start(res)==start(mts)))
  expect_true(all(end(res)==end(mts)))
})
