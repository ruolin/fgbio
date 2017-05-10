---
title: TrimFastq
---

# TrimFastq

## Overview
Group: FASTQ

Trims reads in one or more line-matched fastq files to a specific read length,
and optionally supports dropping reads from all files when one or more reads
is shorter than the desired length.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Values|
|----|----|----|-----------|---------|----------|--------------||input|i|PathToFastq|One or more input fastq files.|Required|Unlimited||
|output|o|PathToFastq|A matching number of output files.|Required|Unlimited||
|length|l|Int|Length to trim reads to.|Required|1||
|exclude|x|Boolean|Exclude reads below the trim length.|Optional|1|false|

