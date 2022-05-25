# ATACseqTFEA

[![platforms](http://bioconductor.org/shields/availability/devel/ATACseqTFEA.svg)](http://bioconductor.org/packages/devel/bioc/html/ATACseqTFEA.html)
[![build](http://bioconductor.org/shields/build/devel/bioc/ATACseqTFEA.svg)](http://bioconductor.org/packages/devel/bioc/html/ATACseqTFEA.html)

Transcription Factor Enrichment Analysis for ATAC-seq

Assay for Transpose-Accessible Chromatin using sequencing (ATAC-seq) is a
technique to assess genome-wide chromatin accessibility by probing open
chromatin with hyperactive mutant Tn5 Transposase that inserts sequencing
adapters into open regions of the genome. ATACseqTFEA is a improvement of
current computational method that detects differential activity of transcription
factors (TFs). ATACseqTFEA not only use the difference of open region
information, but also (or emphasize) the difference of TFs
footprints (cutting sites or insertion sites).

ATACseqTFEA provides an easy, rigorous way to broadly assess TF activity changes
between two conditions.

## Installation

To install this package, start R and enter:

```r
library(BiocManager)
BiocManager::install("ATACseqTFEA")
```

## Documentation

To view documentation of ATACseqTFEA, start R and enter:
```r
browseVignettes("ATACseqTFEA")
```
