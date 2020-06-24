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
import subprocess
import yaml

"""
samtools view -c -b $INFILE -L | samtools view -c > ${INFILE}_human_count.txt
samtools view -s $FRAC -b $INFILE -L hg38_rRNA.bed |samtools view -c > ${INFILE}_rrna_count.txt
samtools view -s $FRAC -b $INFILE sars2 | samtools view -c > ${INFILE}_sars2_count.txt
samtools view -f 4 -s $FRAC -c $INFILE > ${INFILE}_unmapped_count.txt
"""

SAMTOOLS_REGION_COUNT= """
    samtools view -c {bam_in} -L {region_in} 
"""
SAMTOOLS_ONE_REGION_COUNT= """
    samtools view -c {bam_in} {region_in} 
"""
SAMTOOLS_UNMAPPED_COUNT = """
    samtools view -f 4 -c {bam_in}
"""

def _load_cfg(config_file):

    with open(config_file) as in_f:
        data = yaml.load(in_f)
    return(data)

def calculate_alignment_statistics(bam_in, config, output):
    """
        Calculate alignment statistics 
    """
    rrna_hg38 = config["PATHS"]["RRNA_HG38"]
    hg38_genome = config["PATHS"]["GENOME_HG38"]
    sars_genome_name = config["PARAMS"]["SARS_GENOME_NAME"]

    sam_cmd  = SAMTOOLS_REGION_COUNT.format(bam_in=bam_in, region_in=hg38_genome)
    count_one = int(subprocess.check_output(sam_cmd, shell=True).strip())
    sam_cmd  = SAMTOOLS_REGION_COUNT.format(bam_in=bam_in, region_in=rrna_hg38)
    count_two= int(subprocess.check_output(sam_cmd, shell=True).strip())
    sam_cmd  = SAMTOOLS_ONE_REGION_COUNT.format(bam_in=bam_in, region_in=sars_genome_name)
    count_three = int(subprocess.check_output(sam_cmd, shell=True).strip())
    sam_cmd  = SAMTOOLS_UNMAPPED_COUNT.format(bam_in=bam_in) 
    count_four = int(subprocess.check_output(sam_cmd, shell=True).strip())
    with open(output, "w") as out_f:
        out_f.write(str(count_one)  + "\t" + str(count_two) + "\t" + str(count_three) + "\t" + str(count_four) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Calculate COVID19 mapping statistics")
    parser.add_argument("-b","--bam-in", dest="bam_input", required=True)
    parser.add_argument("-c","--config-file",dest="config_file", required=True)
    parser.add_argument("-o","--output-file",dest="output_file", required=True)
    
    args = parser.parse_args() 
    print(args.config_file)
    config = _load_cfg(args.config_file)
    calculate_alignment_statistics(args.bam_input, config, args.output_file)
if __name__=="__main__":
    main()
