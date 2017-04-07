#################################################################################
# The MIT License                                                               #
# Copyright (c) 2017 Fulcrum Genomics LLC                                       #
# Permission is hereby granted, free of charge, to any person obtaining a copy  #
# of this software and associated documentation files (the "Software"), to deal #
# in the Software without restriction, including without limitation the rights  #
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     #
# copies of the Software, and to permit persons to whom the Software is         #
# furnished to do so, subject to the following conditions:                      #
# The above copyright notice and this permission notice shall be included in    #
# all copies or substantial portions of the Software.                           #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN     #
# THE SOFTWARE.                                                                 #
#################################################################################

# R script to generate plots for expected tag family collision rates.

options(warn = -1) # Don't emit warnings, only errors
library(ggplot2)

args    = commandArgs(trailingOnly=T)
metrics = args[1]
output  = args[2]
data = read.table(metrics, header=T)
data$input_dna = factor(data$input_dna)


isize = unique(data$mean_insert_size)
sd    = unique(data$sd_insert_size)

title1 = paste("Family Sizes w/start-stop w/Insert Size mean =", isize, "and sd =", sd)
title2 = paste("Cumulative Family Sizes w/start-stop w/Insert Size mean =", isize, "and sd =", sd)

pdf(output)
  # regular plot
  ggplot(data) + aes(x=family_size, y=fraction_of_families, color=input_dna) + geom_line() +
    scale_y_log10() +
    labs(x="Family size using start/stop alone", y="Fraction of Families at Size", title=title1)

  # cumulative plot
  ggplot(data) + aes(x=family_size, y=cumulative_fraction, color=input_dna) + geom_line() +
    scale_y_continuous(limits=c(0, 1)) +
    labs(x="Family size using start/stop alone", y="Fraction of Families <= Size", title=title2)

dev.off()
