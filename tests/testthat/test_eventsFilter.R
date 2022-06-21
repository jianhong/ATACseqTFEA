test_that("eventsFilter works not correct", {
  se <- SummarizedExperiment(assays=list(counts=data.frame(A=rep(1, 10))),
                             rowRanges = GRanges('1', IRanges(seq.int(10),
                                                              width=1)))
  expect_true(length(eventsFilter(se, rep(c(TRUE, FALSE), 5)))==5)
  expect_true(length(eventsFilter(se, start(se)>8))==2)
  expect_true(length(eventsFilter(se, A==1))==10)
  expect_true(length(eventsFilter(se, "A==1"))==10)
  expect_true(length(eventsFilter(se, A!=1))==0)
  expect_true(length(eventsFilter(se, "A!=1"))==0)
  expect_true(length(eventsFilter(se, seqnames(se)=='1'))==10)
  expect_true(length(eventsFilter(se, seqnames(se)=='chr1'))==0)
  filter <- "A>0"
  expect_true(length(eventsFilter(se, filter))==10)
})
