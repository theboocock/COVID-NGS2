args = commandArgs(trailingOnly=T)

qc_table = read.delim(args[1], sep="\t", header=T)
print(qc_table)
nextseq_table = read.delim(args[2], sep="\t", header=T)
print(nextseq_table)
