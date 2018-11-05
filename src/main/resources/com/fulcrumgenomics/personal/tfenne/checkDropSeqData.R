require(ggplot2)

options(warn = -1) # Don't emit warnings, only errors
library(ggplot2)

args   = commandArgs(trailingOnly=T)
input  = args[1]
output = args[2]
name   = args[3]
cells  = as.integer(args[4])

data = read.table(input, sep="\t", header=T)
data$cell_number = seq(1, nrow(data))

png(output, width=8, height=8, units='in', res=300)

ggplot(data) +
  aes(x=cell_number, y=cumsum(count)) +
  geom_line() +
  geom_vline(xintercept=cells, color="lightgrey") +
  scale_x_sqrt() + 
  labs(x="Number of cells", y="Cumulative UMIs", title=name) +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()
