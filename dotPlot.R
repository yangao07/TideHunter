library(seqinr)
args <- commandArgs(TRUE)


seq=read.fasta(args[1])[[1]]
len_start=as.integer(args[2])
len_end=as.integer(args[3])
w = as.integer(args[4])
s = as.integer(args[5])
m = as.integer(args[6])

pdf("dot.pdf")

dotPlot(seq[len_start:len_end], seq[len_start:len_end], wsize = w, wstep = s, nmatch = m, pch=100,
        main = paste("wsize = ",w, ", wstep = ", s, ", nmatch = ", m, "\n"))

dev.off()
