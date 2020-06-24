# Copyright (c) 2020 Boocock James <james.boocock@ucla.edu>
# Author: Boocock James <james.boocock@ucla.edu>
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
import os
import pysam



def filter_reads(bam_in, min_length, output_bam):
    """
        filter reads by lenth
    """
    samfile = pysam.AlignmentFile(bam_in,"rb")
    filtered_reads = pysam.AlignmentFile(output_bam, "wb", template=samfile)
    for read in samfile:
        if  read.reference_length is None: 
            filtered_reads.write(read)
        elif read.reference_length >= min_length:
            filtered_reads.write(read)
    samfile.close()
    filtered_reads.close()

def main():
    """
        Filter reads by mapping length 
    """
    parser = argparse.ArgumentParser(description="Filter reads by mapping length") 
    parser.add_argument("-b","--bam-in",dest="bam_in", help="Bam input in")
    parser.add_argument("-q","--min-aligned-length",dest="min_length", help="Minimum length for aligned read")
    parser.add_argument("-o","--output",dest="output_bam", help="output bam")
    args = parser.parse_args()
    bam_in = args.bam_in
    min_length= int(args.min_length)
    bam_out = args.output_bam
    filter_reads(bam_in, min_length, bam_out)

if __name__=="__main__":
    main()
