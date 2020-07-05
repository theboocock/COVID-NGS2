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
import subprocess


def create_umi_fasta(bam,umi_primer_site):
    """

    """`


def main():
    parser = argparse.ArgumenntParser(description="Generate fastq containing UMIS")
    parser.add_argument("-i","--input-bam",dest="input_bam",help="input bam")
    parser.add_argument("-o","--output-bam",dest="output_bam",help="output bam")
    parser.add_argument("--umi-primer-site",dest="umi_primer_site",help="Collapse all reads across their primer site",action="store_true",default=False)
    args = parser.parse_args()
    bam = args.input_bam
    out_bam = args.output_bam
    umi_primer_site = args.umi_primer_site
    create_umi_fasta(args.input_bam,args.umi_primer_site)
