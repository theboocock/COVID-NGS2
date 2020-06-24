#!/usr/bin/env python
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
import subprocess
import glob
KRAKEN_PATH=""
KRAKEN_DB=""

KRAKEN_CMD="{exe} --db {kraken_db} --classified-out {out_prefix}.cseqs_#.fastq --unclassified-out {out_prefix}.ucseqs_#.fastq  --threads {threads} --report {report} --output {kraken_output} --paired --gzip-compressed {fastq_in_one} {fastq_in_two}" 
gzip_parallel="""parallel --xapply -j {threads} "gzip -c {{1}} > {{2}}" ::: {in_fastq} ::: {out_fastq_gz}
"""

def run_kraken(fastq_in_one, fastq_in_two, kraken_report, kraken_output, fastq_out_prefix, threads,  KRAKEN_PATH, KRAKEN_DB):
    """
        Run kraken2 software to classify all the reads in each library

    """
    kraken_cmd = KRAKEN_CMD.format(exe=KRAKEN_PATH, kraken_db=KRAKEN_DB, out_prefix=fastq_out_prefix, report=kraken_report, 
            kraken_output=kraken_output, fastq_in_one=fastq_in_one,fastq_in_two=fastq_in_two, threads=threads)
    subprocess.check_call(kraken_cmd, shell=True)     

def main():
    parser = argparse.ArgumentParser(description="Process kraken and filter reads mappig SARS2-COV")
    parser.add_argument("--fastq-in-one",dest="fastq_in_one",help="Read one input")
    parser.add_argument("--fastq-in-two",dest="fastq_in_two",help="Read two input")
    parser.add_argument("--fastq-out-one", dest="fastq_out_r1",help="Read one input")
    parser.add_argument("--fastq-out-two",dest="fastq_out_r2", help="Read two input")
    parser.add_argument("--kraken-report",dest="kraken_report",help="Kraken report")
    parser.add_argument("--kraken-output",dest="kraken_output",help="Kraken output")
    parser.add_argument("--kraken-db",dest="kraken_db", help="Kraken db")
    parser.add_argument("--kraken-path",dest="kraken_path", help="Kraken db")
    parser.add_argument("--fastq-out-prefix",dest="fastq_out_prefix", help="Fastq out prefix")
    parser.add_argument("--threads",dest="threads",default=1)
    args = parser.parse_args()

    fastq_in_one = args.fastq_in_one
    fastq_in_two = args.fastq_in_two
    kraken_report = args.kraken_report
    fastq_out_r1 = args.fastq_out_r1
    fastq_out_r2 = args.fastq_out_r2
    KRAKEN_PATH=args.kraken_path 
    KRAKEN_DB=args.kraken_db
    kraken_output = args.kraken_output 
    fastq_out_prefix = args.fastq_out_prefix 
    threads = args.threads
    run_kraken(fastq_in_one, fastq_in_two, kraken_report, kraken_output,
            fastq_out_prefix,threads,
             KRAKEN_PATH, KRAKEN_DB)

if __name__ == "__main__":
    main()
