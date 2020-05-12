require("ATACseqTFEA") || stop("unable to load Package:ATACseqTFEA")
require("GenomicRanges") || stop("unable to load Package::GenomicRanges")
require("GenomicAlignments") || stop("unalbe to load Package:GenomicAlignments")
require("utils") || stop("unable to load Package:utils")
require("BSgenome.Drerio.UCSC.danRer10") ||
  stop("unable to load Package:BSgenome.Drerio.UCSC.danRer10")
require("Rsamtools") || stop("unable to load Rsamtools")
require("testthat") || stop("unable to load testthat")
test_check("ATACseqTFEA")
