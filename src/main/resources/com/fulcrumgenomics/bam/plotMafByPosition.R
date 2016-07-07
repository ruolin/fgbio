################################################################################
# The MIT License
#
# Copyright (c) 2016 Fulcrum Genomics LLC
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
################################################################################

# Requirements
require(ggplot2)

# Read in some arguments
args  = commandArgs(trailing=T)
input  = args[1]
output = args[2]
name   = args[3]

# Load up and filter the data
data   = read.table(input, header=TRUE, sep="\t")
rows   = length(data$id)
maxId  = max(data$id)
data   = subset(data, data$id == 0 | data$id == maxId | data$group == 'Expected' | data$non_ref_fraction > 0)
ref    = subset(data, data$group != 'Expected')

# Build up the plot title
meanCoverage = formatC(mean(data$depth), digits=2, format="f")
title1 = paste("Non-Ref Fraction By Position for ", name, "\nMean Coverage=", meanCoverage, ";  Covered Loci=", rows, sep="")
title2 = paste("Non-Ref Fraction Distribution At Expected Ref Sites for ", name, "\nMean Coverage=", meanCoverage, ";  Covered Loci=", rows, sep="")

pdf(output, width=11, height=8.5)
ggplot(data) + aes(x=id, y=non_ref_fraction, color=group) + geom_point(size=1, shape=1) +
  scale_y_log10(limits=c(0.0001, 1), breaks=c(0.0001, 0.001, .01, .1, 1)) +
  labs(x="Position along Target Regions", y="Fraction of Non-Reference Observations", title=title1)

ggplot(ref) + aes(x=non_ref_fraction) + geom_histogram(aes(fill=..count..)) +
  scale_y_log10() + scale_x_log10(limits=c(0.0001, 1), breaks=c(0.0001, 0.001, .01, .1, 1)) +
  labs(x="Non-Ref Fraction", y="Sites with Non-Ref Fraction", title=title2)

dev.off()
