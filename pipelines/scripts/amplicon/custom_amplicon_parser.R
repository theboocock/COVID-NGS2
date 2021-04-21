#!/usr/bin/env Rscript

print("HERE")
library(S4Vectors)
library(ShortRead)
library(Biostrings)
library(data.table)

#James see:
#/u/project/kruglyak/jsbloom/COVID_genomes/reference/COVID_384rtPrimers.srt.bed



#Function to trim reads
trimReadsAmp=function(    #read1
                       in.file1,
                        #read2
                        in.file2, 
                        #trimmed read1
                        out.file1,
                        #trimmed read 2
                        out.file2,
                        #DNAStringSet of RT primers
                        oPool.primers,
                        #Variables for Read 1 (RT Read)
                             # UMI length
                             r1.N=3,
                             # Additional trimming after primer
                             r1.Lpad=2,
                             # additional trimming at end of read
                            r1.Rpad=2,

                        # Variables for Read 2 (SSS Read) 
                              # N length (or NSR length) typically random hex
                              r2.N=6,
                              # additional trimming after random seq
                              r2.Lpad=2,
                              # additional trimming at end of read
                              r2.Rpad=2,

                        #total number of records to read at a time into memory
                        nbuffer=1e7,
                        # sometimes just read 2 is sometimes all dark cycles , remove these reads
                        polyG="GGGGGGGGGG",
                        #both read 1 and read 2 should be longer than this amount$
                        min.length.filter=30
                    ){
    fi1 = FastqStreamer(in.file1, nbuffer, readerBlockSize = 1e9, verbose = T)
    fi2 = FastqStreamer(in.file2, nbuffer, readerBlockSize = 1e9, verbose = T)
    idx_iter = 0
    repeat {
        idx_iter= idx_iter + 1 
        rfq1 <- yield(fi1) 
        rfq2 <- yield(fi2) 
         if(length(rfq1) == 0 ) { break }
                            
        cread1 = sread(rfq1)
        q1=quality(rfq1)
        cread2 = sread(rfq2)
        q2=quality(rfq2)
        
        w1=width(cread1)
        w2=width(cread2)
        cread1=cread1[w1>min.length.filter & w2>min.length.filter]
        q1=q1[w1>min.length.filter & w2>min.length.filter]
        cread2=cread2[w1>min.length.filter & w2>min.length.filter]
        q2=q2[w1>min.length.filter & w2>min.length.filter]
    
        for(oname in names(oPool.primers)[1:length(oPool.primers)]){
					r1.primer.N = length(oPool.primers[[oname]])
					m1=which.isMatchingStartingAt(oPool.primers[[oname]], cread1, 
															starting.at = (r1.N+1):(r1.N+3), 
															max.mismatch = 1, 
															follow.index = T, 
															with.indels = F)
				  to.keep=!is.na(m1)
				  ### Let's seee if we can filter the reads
				  cread1s=cread1[to.keep]
				  q1s=q1[to.keep]
				  m1k=m1[to.keep]
				  

				  umi1=narrow(cread1s, start=m1k-r1.N, end=m1k-1)
				 
				  r1w=width(cread1s) 

				  r1.start=m1k+r1.primer.N+r1.Lpad
				  r1.end=r1w-r1.Rpad

				  r1.trimmed=narrow(cread1s, start=r1.start, end=r1.end)
				  q1.trimmed=narrow(q1s, start=r1.start, end=r1.end)

			   
				  cread2s= cread2[to.keep]
				  q2s= q2[to.keep]


				  r2w=width(cread2s) 


				  umi2=narrow(cread2s, start=1, end=r2.N)
				
				  #clip a few bases on r2 to help with systematic errors there 
							  
				  r2.start=r2.N+r2.Lpad
				  r2.end=r2w-r2.Rpad

				  r2.trimmed=narrow(cread2s, start=r2.start,end=r2.end)
				  q2.trimmed=narrow(q2s, start=r2.start,end=r2.end)
				  
				  #dump reads where read2 is  polyG
				  r2.to.keep=vcountPattern(polyG, r2.trimmed)==0
				  #could comment out next 2 lines, mostly for debugging
				  #print('polyG')
				  #print(sum(!r2.to.keep))
				  umi1=umi1[r2.to.keep]
				  r1.trimmed=r1.trimmed[r2.to.keep]
				  q1.trimmed=q1.trimmed[r2.to.keep]
				  umi2=umi2[r2.to.keep]
				  r2.trimmed=r2.trimmed[r2.to.keep]
				  q2.trimmed=q2.trimmed[r2.to.keep]

				  #could comment out next 6 lines, mostly for debugging
				  u.umis=rle(sort(as.character(xscat(umi1))))
				  bdf=data.frame(barcode=u.umis$values, count=u.umis$lengths, stringsAsFactors=F)
				  #omatch[[oname]]=bdf
				  #print(oname)
				  #print(nrow(bdf))
				  #print(sum(bdf$count))
				  if(length(r1.trimmed)>0){
				  cname=xscat(BStringSet(rep(oname,length(r1.trimmed))), ':', BStringSet(1:length(r1.trimmed)), ':', BStringSet(umi1), ':', BStringSet(umi2),":",BStringSet(rep(idx_iter,length(r1.trimmed))))
				  writeFastq( ShortReadQ(sread=r1.trimmed, quality=q1.trimmed,  id= cname),  file=out.file1, mode='a', full=FALSE, compress=TRUE)
				  writeFastq( ShortReadQ(sread=r2.trimmed, quality=q2.trimmed, id= cname),  file=out.file2, mode='a', full=FALSE, compress=TRUE)
              }
    }
    }

    #print("HERE")
    close(fi1)
    close(fi2)
}
trimReadsArctic=function(    #read1
                       in.file1,
                        #read2
                        in.file2, 
                        #trimmed read1
                        out.file1,
                        #trimmed read 2
                        out.file2,
                        #DNAStringSet of RT primers
                        primers,
                        #Variables for Read 1 (RT Read)
                             # UMI length
                             r1.N=8,
                             # Additional trimming after primer
                             r1.Lpad=3,
                             # additional trimming at end of read
                            r1.Rpad=3,
                        # Variables for Read 2 (SSS Read) 
                              # N length (or NSR length) typically random hex
                              r2.N=6,
                              # additional trimming after random seq
                              r2.Lpad=3,
                              # additional trimming at end of read
                              r2.Rpad=3,
                        #total number of records to read at a time into memory
                        nbuffer=1e7,
                        # sometimes just read 2 is sometimes all dark cycles , remove these reads
                        polyG="GGGGGGGGGG",
                        #both read 1 and read 2 should be longer than this amount
                        min.length.filter=30,
                        RTprimerExtraL='CTCGATGGAG'
                        
                    ){
    fi1 = FastqStreamer(in.file1, nbuffer, readerBlockSize = 1e9, verbose = T)
    fi2 = FastqStreamer(in.file2, nbuffer, readerBlockSize = 1e9, verbose = T)
    idx_iter = 0
      repeat {
        idx_iter = idx_iter + 1

        rfq1 <- yield(fi1) 
        rfq2 <- yield(fi2) 
      
         if(length(rfq1) == 0 ) { break }
                            
        cread1 = sread(rfq1)
        q1=quality(rfq1)
        cread2 = sread(rfq2)
        q2=quality(rfq2)
        
        w1=width(cread1)
        w2=width(cread2)
        cread1=cread1[w1>min.length.filter & w2>min.length.filter]
        q1=q1[w1>min.length.filter & w2>min.length.filter]
        cread2=cread2[w1>min.length.filter & w2>min.length.filter]
        q2=q2[w1>min.length.filter & w2>min.length.filter]
        for(oname in names(primers)) { #oPool.primers)[1:length(oPool.primers)]){
                
               r1.primer.N = primers[[oname]]$length #length(oPool.primers[[oname]])
               rseqplus=paste0(RTprimerExtraL, primers[[oname]]$seq)
               m1=which.isMatchingStartingAt(rseqplus, cread1, 
                                                        starting.at = r1.N+1, #(r1.N-1):(r1.N+1), 
                                                        max.mismatch = 1, 
                                                        follow.index = T, 
                                                        with.indels = F)
              to.keep=!is.na(m1)
              ### Let's seee if we can filter the reads
              cread1s=cread1[to.keep]
              q1s=q1[to.keep]
              m1k=m1[to.keep]
              
              umi1=narrow(cread1s, start=1, end=r1.N)
             
              r1w=width(cread1s) 
              r1.start=m1k+r1.primer.N+r1.Lpad
              r1.end=r1w-r1.Rpad
              r1.trimmed=narrow(cread1s, start=r1.start, end=r1.end)
              q1.trimmed=narrow(q1s, start=r1.start, end=r1.end)
           
              cread2s= cread2[to.keep]
              q2s= q2[to.keep]
              r2w=width(cread2s)
             # writeXStringSet(cread2s, filepath='~/Desktop/amp13.txt')
             
              umi2=narrow(cread2s, start=1, end=r2.N)
            
              #clip a few bases on r2 to help with systematic errors there 
                          
              r2.start=r2.N+r2.Lpad
              r2.end=r2w-r2.Rpad
              r2.trimmed=narrow(cread2s, start=r2.start,end=r2.end)
              q2.trimmed=narrow(q2s, start=r2.start,end=r2.end)
              
              #dump reads where read2 is  polyG
              r2.to.keep=vcountPattern(polyG, r2.trimmed)==0
              #could comment out next 2 lines, mostly for debugging
              #print('polyG')
              #print(sum(!r2.to.keep))
              umi1=umi1[r2.to.keep]
              r1.trimmed=r1.trimmed[r2.to.keep]
              q1.trimmed=q1.trimmed[r2.to.keep]
              umi2=umi2[r2.to.keep]
              r2.trimmed=r2.trimmed[r2.to.keep]
              q2.trimmed=q2.trimmed[r2.to.keep]
              #could comment out next 6 lines, mostly for debugging
              u.umis=rle(sort(as.character(xscat(umi1))))
              bdf=data.frame(barcode=u.umis$values, count=u.umis$lengths, stringsAsFactors=F)
              #omatch[[oname]]=bdf
              #print(oname)
              #print(nrow(bdf))
              #print(sum(bdf$count))
              
              if(length(r1.trimmed)>0){
	          cname=xscat(BStringSet(rep(oname,length(r1.trimmed))), ':', BStringSet(1:length(r1.trimmed)), ':', BStringSet(umi1), ':', BStringSet(umi2),":",BStringSet(rep(idx_iter,length(r1.trimmed))))
              writeFastq( ShortReadQ(sread=r1.trimmed, quality=q1.trimmed,  id= cname),  file=out.file1, mode='a', full=FALSE, compress=TRUE)
              writeFastq( ShortReadQ(sread=r2.trimmed, quality=q2.trimmed, id= cname),  file=out.file2, mode='a', full=FALSE, compress=TRUE)
              }
        }
    }
    close(fi1)
    close(fi2)
} 



# some input stuff for running in R
args = commandArgs(trailingOnly=T)
fastq_one= args[1] 
print(fastq_one)
fastq_two = args[2] 
fastq_out_one = args[3]
fastq_out_two = args[4]
opool.primers = args[5]
type = args[6]
arctic.primers = args[7]
all.primers = args[8]
sample_list = args[9] 
oPool.primers.table=read.delim(opool.primers, header=F, stringsAsFactors=F)
oPool.primers=DNAStringSet(oPool.primers.table[,5])
names(oPool.primers) = oPool.primers.table[,4]
#James see /u/project/kruglyak/jsbloom/COVID_genomes/outAMP/
#experiment name
#clinic name
#sequence name
# stuff starts happening here
#loop through fastq files 
primers = read.table(arctic.primers, header=T, sep="\t", stringsAsFactors=F)
#print(head(primers))
primers = primers[grep('RIGHT', primers$name),]
if (type == "arctic"){
	primers = read.table(arctic.primers, header=T, sep="\t", stringsAsFactors=F)
	primers = primers[grep('RIGHT', primers$name),]
	primers=split(primers, primers$name)
	trimReadsArctic(#read1
                   fastq_one,
                    #read2
                    fastq_two, 
                    #trimmed read1
                    fastq_out_one,
                    #trimmed read 2
                    fastq_out_two,
                    #DNAStringSet of RT primers
                    primers)
}else if (type == "amp"){
	oPool.primers.table=read.delim(opool.primers, header=F, stringsAsFactors=F)
	oPool.primers=DNAStringSet(oPool.primers.table[,5])
	names(oPool.primers) = oPool.primers.table[,4]

	trimReadsAmp(#read1
                   fastq_one,
                    #read2
                    fastq_two, 
                    #trimmed read1
                    fastq_out_one,
                    #trimmed read 2
                    fastq_out_two,
                    #DNAStringSet of RT primers
                    oPool.primers)

}else if (type == "amp_complicated"){
    library(stringr)
    x = str_split(sample_list,",")
    x = as.numeric(x[[1]])
    #print(x)
    #oPool.primers.table = oPool.primers.table[which(oPool.primers.table[,3] %in% x),]
    #print(nrow(oPool.primers.table))
    #oPool.primers.table = oPool.primers.table[!duplicated(oPool.primers.table[,2]),]
    #print(nrow(oPool.primers.table))
    #oPool.primers = DNAStringSet(oPool.primers.table[,2])
    #names(oPool.primers) = oPool.primers.table[,1]
	
    ### Send the arctic information through. 
    # 2 is arctic, 1,3,4 are amplicon. 
    print(x)
    if (any(c(1,3,4) %in% x)){
        y  = x[x!=2]        
        print(y)
        oPool.primers.table = read.delim(all.primers, header=F, sep="\t", stringsAsFactors=F)
        oPool.primers.table = oPool.primers.table[which(oPool.primers.table[,3] %in% y),]
        oPool.primers.table = oPool.primers.table[!duplicated(oPool.primers.table[,2]),]
        #write.table(oPool.primers.table, file="test_data1.txt", quote=F, row.names=F,col.names=T,sep="\t")
        #oPool.primers.table = oPool.primers.table[oPool.primers.table[,1] == "c19_1452_1477",]
        #write.table(oPool.primers.table, file="test_data.txt", quote=F, row.names=F,col.names=T,sep="\t")
        oPool.primers = DNAStringSet(oPool.primers.table[,2])
        #print(oPool.primers.table)
        names(oPool.primers) = oPool.primers.table[,1]
        trimReadsAmp(#read1
                       fastq_one,
                        #read2
                        fastq_two, 
                        #trimmed read1
                        fastq_out_one,
                        #trimmed read 2
                        fastq_out_two,
                        #DNAStringSet of RT primers
                        oPool.primers)
    }
    if (2 %in% x){
        y  = x[x==2]        
        oPool.primers.table = read.delim(all.primers, header=F, sep="\t", stringsAsFactors=F)
        oPool.primers.table = oPool.primers.table[which(oPool.primers.table[,3] %in% y),]
        oPool.primers.table = oPool.primers.table[!duplicated(oPool.primers.table[,2]),]
        oPool.primers.arctic = DNAStringSet(oPool.primers.table[,2])
        write.table(oPool.primers.table, file="out1.txt", quote=F,row.names=F,col.names=T)
        primers.arctic = data.frame(seq=oPool.primers.table[,2],name=oPool.primers.table[,1],length=sapply(oPool.primers.table[,2],nchar))
        primers.arctic = split(primers.arctic,primers.arctic$name)
        trimReadsArctic(#read1
                       fastq_one,
                        #read2
                        fastq_two, 
                        #trimmed read1
                        fastq_out_one,
                        #trimmed read 2
                        fastq_out_two,
                        #DNAStringSet of RT primers
                        primers.arctic)
    }
    ## REMOVE DUPLICETASE
    #print(oPool.primers.table)

} 

# following code will report a per-amplicon summary from the generated bam files above
#library(rbamtools)
#bam.files=list.files(out.dir, pattern='.bam$') 
#for(bam.file in bam.files) {
#
#reader = bamReader(paste0(out.dir,bam.file),idx=TRUE)
## get total expected oligos
#contigs = (refSeqDict(getHeaderText(reader)))@SN
#contigs.length=(refSeqDict(getHeaderText(reader)))@LN
#br=bamRange(reader, coords=c(0, 1, contigs.length[1]))
#
#adepth=alignDepth(br)
#png(filename=paste0(out.dir,bam.file,'.png'), width=1024, height=512)
#plotAlignDepth(adepth, main=paste(bam.file, 'median depth=', median(adepth@depth)))
#dev.off()
#
#dbr=data.frame(br)
#
#primer.count=rle(sort(tstrsplit(dbr$name, ':')[[1]]))
#primer.count=data.frame(primer=primer.count$values, count=primer.count$lengths)
#
#opOrig=oPool.primers.table
#names(opOrig)=c('genome', 'start','stop', 'name', 'seq')
#oPmerge=merge(opOrig, primer.count, by.x='name', by.y='primer',all.x=T, sort=F)
#oPmerge$count[is.na(oPmerge$count)]=0
#oPmerge=oPmerge[order(oPmerge$start),]
