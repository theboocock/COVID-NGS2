
args = commandArgs(trailingOnly=T)
qc = read.delim(args[1],sep="\t", header=T)
out_dir = args[2] 
in_dir = args[3] 
dir.create(out_dir)
qc = qc[qc$sample_type == "amp",]
for (num in  qc$sample){
    print(num)
    files_to_copy = list.files(in_dir,paste(num,".*",sep=""),full.names=T) 
    for (file_in in files_to_copy){
        print(file_in)
        file.copy(file_in, out_dir)
    }
}

write.table(qc, file=paste(out_dir,"/qc_report.tsv",sep=""), quote=F, row.names=F,col.names=T,sep="\t")
