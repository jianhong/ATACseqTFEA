test_that("expandBindingSites works not correct", {
  bindingSites <- GRanges('chr1', IRanges(1000, width = 10))
  bs <- expandBindingSites(bindingSites = bindingSites,
                           proximal=40L,
                           distal=40L,
                           gap=10L)
  expect_true(width(bs$bindingSitesWithGap)==width(bindingSites) + 20)
  expect_true(width(bs$bindingSitesWithProximal)==width(bindingSites) + 100)
  expect_true(width(bs$bindingSitesWithProximalAndGap)==
                width(bindingSites) + 120)
  expect_true(width(bs$bindingSitesWithDistal)==width(bindingSites) + 200)
})
