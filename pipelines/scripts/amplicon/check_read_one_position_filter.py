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

def filter_reads(bams, output_bam, summary_file,arctic_hash):
    """
        Filter reads
    """
    pysam.sort("-n","-o","test.bam", bams)
    samfile = pysam.AlignmentFile("test.bam","rb")
    out_bam = pysam.AlignmentFile("test1.bam", "wb", template=samfile)
    last_read = ""
    read_one = False
    read_two = False
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
                # can't match if there isn't a read one
                if read_one:
                    if "nCoV" in tmp_read1.query_name:
                        read_short = tmp_read1.query_name.split(":")[0]
                        arctic_match_end = (arctic_hash[read_short][1])
                        position = tmp_read1.reference_end
                        if position is not None:
                            #print(position)
                            #print(arctic_match_end + 6)
                            if position == (arctic_match_end + 6):
                                out_bam.write(tmp_read1)
                                if read_two:
                                    out_bam.write(tmp_read2)
                                    count[0] += 1 
                                count[0] += 1 
                    #must be amplicon
                    else:
                        start = int(tmp_read1.query_name.split("_")[1])
                        position = tmp_read1.reference_end
                        if position != None:
                            position = position + 2 
                            if start == position:
                                out_bam.write(tmp_read1)
                                if read_two:
                                    out_bam.write(tmp_read2)
                                    count[0] += 1 
                                count[0] += 1 
                last_read = current_read
                read_one = False
                read_two = False

            if read.is_read1:
                read_one = True
                tmp_read1 = read
            if read.is_read2:
                read_two = True
                tmp_read2 = read
            total_count +=1
        if (total_count != 0):
            count[0] = count[0] / float(total_count)
            count[1] = 1 - count[0] 
        summary_out.write(str(count[0]) + " " + str(count[1]) + "\n")
    out_bam.close()
    pysam.sort("-o",output_bam, "test1.bam")


def get_arctic_hash(arctic_bed):
    arctic_hash = {}
    with open(arctic_bed) as arctic_in:
        for line in arctic_in:
            line_s=line.split("\t")
            start = line_s[8]
            end = line_s[9]
            name = line_s[0]
            arctic_hash[name]=[int(start),int(end)]

    return arctic_hash

def main():
    parser = argparse.ArgumentParser(description="Read lines")
    parser.add_argument("-i","--input-bam",dest="bam",help="Bam input")
    parser.add_argument("-o","--output-bam",dest="output_bam",help="BAM output")
    parser.add_argument("-s","--summary",dest="summary_file",help="summary")
    parser.add_argument("-a","--arctic-positions",dest="arctic_positions",help="arctic")
    args = parser.parse_args()
    bams = args.bam
    out_bam = args.output_bam
    summary_file = args.summary_file
    arctic_hash = get_arctic_hash(args.arctic_positions)
    filter_reads(bams, out_bam,summary_file,arctic_hash)

if __name__ == "__main__":
    main()
