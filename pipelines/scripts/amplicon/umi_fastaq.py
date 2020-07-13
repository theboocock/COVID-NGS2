#!/usr/bin/env python

# Copyright (c) 2020 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import argparse
import pysam

### fgbio ### collapse and remap reads ###

def umi_tag(bam,out_bam):
    """

    """
    ### TODO: create arctic UMIs and map them back to the original primer sequence or something.
    samfile = pysam.AlignmentFile(bam)
    samfile_out = pysam.AlignmentFile(out_bam, "wb", template=samfile)
    only_pairs= {} 
    for read in samfile:
            # TOOD: umi just  becomes the primer name
            # TODO: umi becomes the 3-bp UMI and the NSR.

        split_name = (read.query_name.split(":"))
        primer_name = split_name[0]
        read.set_tag(tag="PN",value=primer_name, value_type="Z")
        umi = split_name[2]
        rh = split_name[3]
        umi_string = umi + rh 
        read.set_tag(tag="UB",value=umi_string,value_type="Z")
        try:
            only_pairs[read.query_name].append(read)
        except:
            only_pairs[read.query_name] = [read]
        if len(only_pairs[read.query_name]) == 2:
            read_one=False
            read_two=False
            for value in only_pairs[read.query_name]:
                if value.is_read1:
                    read_one = True
                if value.is_read2:
                    read_two = True
            if read_one and read_two: 
                for value in only_pairs[read.query_name]:
                    samfile_out.write(value)
    

sort_bam=""" fgbio SortBam -s TemplateCoordinate -i {0} -o {1}""" 
set_mate_info=""" fgbio SetMateInformation -i {0} -o {1} """ 
group_read_by_umi =""" fgbio GroupReadsByUmi -i {0} -o {1} -t UB -s edit --edits 0 """
molecular_consensus=""" fgbio  CallMolecularConsensusReads -i {0} -o {1} -t MI -M 1 """  
bam_to_fastq=""" bedtools bamtofastq -i {0} -fq {1} -fq2 {2} """ 
quick_align="""bwa mem -t 1 -R "{0}" {1} {2} {3} | samtools sort > {4} && samtools index {4}"""
import subprocess
import shutil
def sort_and_correct_reads(in_bam,out_bam, fastq_one, fastq_two,reference,read_group, summarise_molecular_info):
    """
        Sort,get molecular consensus and finally call variants. 
    """
    sort_bam_cmd = sort_bam.format(in_bam,"test1.bam") 
    subprocess.check_call(sort_bam_cmd,shell=True)
    set_mate_info_cmd = set_mate_info.format("test1.bam","test2.bam")
    subprocess.check_call(set_mate_info_cmd,shell=True)
    group_read_by_umi_cmd = group_read_by_umi.format("test2.bam","test3.bam")
    subprocess.check_call(group_read_by_umi_cmd,shell=True)
    ## TODO: damn order
    shutil.copy("test3.bam", summarise_molecular_info)
    #subprocess.check_call("samtools index {0}".format(summarise_molecular_info),shell=True)
    molecular_consensus_cmd = molecular_consensus.format("test3.bam","test4.bam")
    subprocess.check_call(molecular_consensus_cmd,shell=True)
    bam_to_fastq_cmd = bam_to_fastq.format("test4.bam", fastq_one,fastq_two)
    subprocess.check_call(bam_to_fastq_cmd,shell=True)
    ## ADD BACK READ GROUP
    quick_align_cmd = quick_align.format(read_group, reference, fastq_one, fastq_two, out_bam) 
    subprocess.check_call(quick_align_cmd, shell=True)

def main():
    parser = argparse.ArgumentParser(description="Generate fastq containing UMIS")
    parser.add_argument("-i","--input-bam",dest="input_bam",help="input bam")
    parser.add_argument("-o","--output-bam",dest="output_bam",help="output bam")
    parser.add_argument("--output-summarise-molecular-info",dest="summarise_molecular_info",help="Summarise molecular info")
    parser.add_argument("--fastq-one-out",dest="fastq_one_out")
    parser.add_argument("--fastq-two-out",dest="fastq_two_out")
    parser.add_argument("--read-group", dest="read_group")
    parser.add_argument("--reference",dest="reference_genome")
    args = parser.parse_args()
    bam = args.input_bam
    out_bam = args.output_bam
    print(args.input_bam)
    umi_tag(args.input_bam,"test.bam")
    sort_and_correct_reads("test.bam", args.output_bam, args.fastq_one_out,args.fastq_two_out,args.reference_genome,args.read_group,args.summarise_molecular_info)



if __name__ == "__main__":
    main()
