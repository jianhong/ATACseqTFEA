# File list

## Motifs saved in R objects.

- `best_curated_Human.rds`: The `best_curated_Human` is a list of TF motifs
downloaded from [TFEA github](https://github.com/Dowell-Lab/TFEA).
There are 1279 human motifs in the data set.

- `cluster_PWMs.rds`: The `cluster_PWMs` is a list of non-redundant TF 
motifs downloaded from
[DeepSTARR](https://github.com/bernardo-de-almeida/motif-clustering).
There are 6502 motifs in the data set.

- `PWMatrixList.rds`: The `PWMatrixList` is a collection of jasper2018,
jolma2013 and cisbp_1.02 from package motifDB (v 1.28.0) and merged by distance
smaller than 1e-9 calculated by MotIV::motifDistances function (v 1.42.0).
The merged motifs were exported by motifStack (v 1.30.0).

## Command line scripts

- `sample_scripts.R`: This script will run the TFEA step by step and save the
intermediate results in R objects. 

## Unit test or documentation used R objects.

- `bindingSites.rds`: The sample outputs of the function `prepareBindingSites`
- `res.rds`: The sample outputs of the function `TFEA`

## Sample bam files and their index files.

The bam files are insertion site shifted files.
KD means Knockdown; WT means wild type.

- `KD.shift.rep1.bam`
- `KD.shift.rep1.bam.bai`
- `KD.shift.rep2.bam`
- `KD.shift.rep2.bam.bai`
- `WT.shift.rep1.bam`
- `WT.shift.rep1.bam.bai`
- `WT.shift.rep2.bam`
- `WT.shift.rep2.bam.bai`
