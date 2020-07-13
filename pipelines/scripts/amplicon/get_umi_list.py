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

import pysam
import argparse

###
###
sort_bam=""" fgbio SortBam -s TemplateCoordinate -i {0} -o {1}"""
set_mate_info=""" fgbio SetMateInformation -i {0} -o {1} """
group_read_by_umi =""" fgbio GroupReadsByUmi -i {0} -o {1} -t UB -s edit --edits 2 """
def read_bam_input(bam,sample, output_file):
    """
        Read bam input
    """
    samfile = pysam.AlignmentFile(bam,"rb")
    with open(output_file,"w") as out_f:
        for read in samfile:
            umi = read.get_tag("UB") 
            md = read.get_tag("MI")
            pos = read.reference_start
            out_f.write(sample + "\t" +str(pos) + "\t" + str(umi) + "\t" + str(md) +"\n") 

def main():
    parser= argparse.ArgumentParser(description="Calculate the read-depth UMI")
    parser.add_argument("--bam",dest="bam",help="BAM input file")
    parser.add_argument("--sample",dest="sample", help="Sample input")
    parser.add_argument("--output-df",dest="output_df",help="Output dataframe")
    args = parser.parse_args()
    bam_input = args.bam
    read_bam_input(bam_input,args.sample,args.output_df)

if __name__=="__main__":
    main()
