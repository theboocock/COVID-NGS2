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

def filter_reads(bams, output_bam, summary_file):
    """
        Filter reads
    """
    pysam.sort("-n","-o","test.bam", bams)
    samfile = pysam.AlignmentFile("test.bam","rb")
    out_bam = pysam.AlignmentFile("test1.bam", "wb", template=samfile)
    last_read = ""
    read_one = False
    count  = [0,0]
    total_count = 0
    first = True
    with open(summary_file,"w") as summary_out:
        for read in samfile:
            query_name = read.query_name
            current_read = query_name 
            if first:
                last_read = current_read 
                first = False
            if current_read != last_read:
                if read_one:
                    start = int(tmp_read1.query_name.split("_")[1])
                    position = tmp_read1.reference_end
                    if position != None:
                        position = position + 2 
                        if start == position:
                            out_bam.write(tmp_read1)
                            if read_two:
                                out_bam.write(tmp_read2)
                            count[0] += 1 
                        else:
                            count[1] +=1
                        total_count +=1
                        last_read = current_read
                read_one = False
                read_two = False
            if read.is_read1:
                read_one = True
                tmp_read1 = read
            if read.is_read2:
                read_two = True
                tmp_read2 = read
        count[0] = count[0] / float(total_count)
        count[1] = count[1] / float(total_count)
        summary_out.write(str(count[0]) + " " + str(count[1]) + "\n")
    out_bam.close()
    pysam.sort("-o",output_bam, "test1.bam")


def main():
    parser = argparse.ArgumentParser(description="Read lines")
    parser.add_argument("-i","--input-bam",dest="bam",help="Bam input")
    parser.add_argument("-o","--output-bam",dest="output_bam",help="BAM output")
    parser.add_argument("-s","--summary",dest="summary_file",help="summary")
    args = parser.parse_args()
    bams = args.bam
    out_bam = args.output_bam
    summary_file = args.summary_file
    filter_reads(bams, out_bam,summary_file)

if __name__ == "__main__":
    main()
