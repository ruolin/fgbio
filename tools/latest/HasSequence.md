---
title: HasSequence
---

# HasSequence

## Overview
Group: Personal

Searches for DNA sequences in the read pairs.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Values|
|----|----|----|-----------|---------|----------|--------------|
|input|i|PathToBam|Input SAM or BAM.|Required|1||
|output|o|PathToBam|Output SAM or BAM.|Required|1||
|sequences|m|Path|File containing the DNA sequences to search for.|Required|1||
|metrics|s|Path|Sequence match histogram written to this file.|Required|1||
|max-mismatches||Int|Maximum mismatches for matching a sequence.|Optional|1|1|
|minimum-base-quality|q|Int|Minimum base quality. Any bases falling below this quality will be considered a mismatch even in the bases match.|Optional|1|0|
|include-no-calls||Boolean|Count no calls as mismatches unless both bases are no calls.|Optional|1|false|
|start-only||Boolean|Search for mismatches at the start only|Optional|1|true|
|has-sequence-tag||String|The tag to store the result.|Optional|1|XW|

