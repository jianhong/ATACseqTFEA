require("ATACseqTFEA") || stop("unable to load Package:ATACseqTFEA")
require("BSgenome.Drerio.UCSC.danRer10") ||
  stop("unable to load Package:BSgenome.Drerio.UCSC.danRer10")
require("testthat") || stop("unable to load testthat")
test_check("ATACseqTFEA")
