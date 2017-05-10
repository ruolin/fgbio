---
title: EstimatePoolingFractions
---

# EstimatePoolingFractions

## Overview
Group: SAM/BAM

Examines a pooled sample and estimates the fraction of each constituent sample.
Uses a VCF of known genotypes for the samples within the mixture along with a
BAM of sequencing data derived from the pool.  Performs a multiple regression
for the alt allele fractions at each SNP locus, using as inputs the individual
sample's genotypes.  Only SNPs that are bi-allelic within the pooled samples are
used.

Various filtering parameters can be used to control which loci are used:
 --intervals will restrict analysis to variants within the described intervals
 --min-genotype-quality will filter out any site with any genotype with GQ < n
 --min-mean-sample-coverage requires that the coverage of a site in the BAM be
     >= min-mean-sample-coverage * n_samples
 --min-mapping-quality filters out reads in the BAM with MQ < n
 --min-base-quality filters out bases in the BAM with Q < n

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Values|
|----|----|----|-----------|---------|----------|--------------||vcf|v|PathToVcf|VCF of individual sample genotypes.|Required|1||
|bam|b|PathToBam|Path to BAM file of sequencing data.|Required|1||
|output|o|FilePath|Output file to write with pooling fractions.|Required|1||
|intervals|l|PathToIntervals|Zero or more set of regions to restrict analysis to.|Optional|Unlimited||
|samples|s|String|Optional subset of samples from VCF to use.|Optional|Unlimited||
|non-autosomes|n|String|Non-autosomal chromosomes to avoid.|Required|Unlimited|M, chrM, MT, X, chrX, Y, chrY|
|min-genotype-quality|g|Int|Minimum genotype quality. Use -1 to disable.|Optional|1|30|
|min-mean-sample-coverage|c|Int|Minimum (sequencing coverage @ SNP site / n_samples).|Optional|1|6|
|min-mapping-quality|m|Int|Minimum mapping quality.|Optional|1|20|
|min-base-quality|q|Int|Minimum base quality.|Optional|1|5|

