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

__BCFTOOLS_FILT_TEMPLATE__ = """
    bcftools view {0} -e " INFO/DP {1} | QUAL < {2}" > filt.vcf" 
"""


## TODO: remove SNPS ONLY and biallelic calls
__BCFTOOLS_SAMPLE_TEMPLATE__ = """
    bcftools view -m2 -M2 -v snps {0} -s {1} | bgzip -c > sample.vcf.gz && tabix -p vcf sample.vcf.gz 
"""

__BCFTOOLS_QUERY_TEMPLATE__= """
    bcftools view {0} -i " FORMAT/DP < {1} | QUAL < {2}" | bcftools query -f '%CHROM\\t%POS0\\t%END\\t%ID\\n' > mask.bed 
"""

__BCFTOOLS_CONSENSUS_TEMPLATE__="""
    bcftools consensus -m merged_masked.bed -s {0} -f {1} sample.vcf.gz > {2} 
"""
__BEDTOOLS_MASKING_STRAND__="""
    bedtools genomecov -ibam {0} -d -strand {1} > {2} 
"""
### SET GT ###
#bcftools +setGT sample.vcf.gz -- -t q --new-gt . -i 'FMT/DP < 5'
#

__BIOAWK_COV__="""
        cat {0} | bioawk -c 'fastx' '{{a=$seq; print 1 - gsub(/N/,"",$seq)/length(a)}}' > {1} 
    """

__BEDTOOLS_MERGE__= """
    cat {} {} | sort -k 1,1 -k2,2n | uniq | bedtools merge -i stdin > merged_masked.bed
"""


def get_masking_bed(bam_input, min_depth, max_strand_prop):

    plus_cmd = __BEDTOOLS_MASKING_STRAND__.format(bam_input ,"+", "plus.bed")
    subprocess.check_call(plus_cmd, shell=True)
    minus_cmd = __BEDTOOLS_MASKING_STRAND__.format(bam_input,"-","minus.bed")
    subprocess.check_call(minus_cmd,shell=True)
    with open("plus.bed") as plus_in:
        with open("minus.bed") as minus_in:
            with open("mask.bed","w") as mask_out:
                for line1, line2 in zip(plus_in, minus_in):
                    line1s = line1.split()
                    line2s = line2.split()
                    chrom = line1s[0]
                    end = int(line1s[1])
                    start = int(line1s[1]) - 1
                    plusc = int(line1s[2])
                    minc = int(line2s[2])
                    depth = plusc + minc
                    if depth < min_depth:
                        mask_out.write(chrom + "\t" + str(start) + "\t" + str(end)+"\n")
                        continue
                    max_c = max(minc, plusc)
                    if max_c/(depth*1.0) > max_strand_prop:
                        mask_out.write(chrom + "\t" + str(start) + "\t" + str(end) + "\n")

def get_sample_vcf(vcf_gz, sample):
    sample_cmd = __BCFTOOLS_SAMPLE_TEMPLATE__.format(vcf_gz, sample)
    print(sample_cmd)
    subprocess.check_call(sample_cmd,shell=True)

def get_consensus_fasta(reference_genome, sample, output_consensus, output_cov, snps_to_exclude_bed):

    merged_bed_cmd = __BEDTOOLS_MERGE__.format("mask.bed",snps_to_exclude_bed)
    print(merged_bed_cmd)
    subprocess.check_call(merged_bed_cmd,shell=True)

    con_cmd = __BCFTOOLS_CONSENSUS_TEMPLATE__.format(sample, reference_genome, "con.fasta")
    subprocess.check_call(con_cmd,shell=True)
    with open("con.fasta") as con:
        with open(output_consensus,'w') as out_f:
            for line in con:
                if ">" in line:
                    out_f.write(">" + sample + "\n")
                else:
                    out_f.write(line)
    cov_cmd =__BIOAWK_COV__.format(output_consensus, output_cov)
    subprocess.check_call(cov_cmd, shell=True)

def main():
    parser = argparse.ArgumentParser(description="Process joint VCF and generate consensus fastas from biallelic sites only")
    print("Welcome to the sars-cov-2 variannt caller")
    parser.add_argument("-s","--sample",dest="sample_id", help="Sample id")
    parser.add_argument("-v","--vcf", dest="vcf_gz", help="VCF input GZ")
    parser.add_argument("-d","--min-depth", dest="min_depth", help="Minimum depth per sample", default=5)
    parser.add_argument("-p","--max-strand-prop",dest="max_strand_prop", help="Minimum quality",default=1.0)
    parser.add_argument("-b","--bcftools-path", dest="bcftools_path", help="Bcftools path")
    parser.add_argument("-r","--reference-genome",dest="reference_genome", help="Reference genome",required=True)
    parser.add_argument("-m","--masked-bed",dest="masked_bed", help="Masked bed file", required=True)
    parser.add_argument("-i","--bam",dest="bam_input", help="Bam input file") 
    parser.add_argument("-o","--output-fasta",dest="output_fasta", help="output fasta", required=True)
    parser.add_argument("-c","--coverage-output",dest="coverage_output", help="Coverage output", required=True)
    args = parser.parse_args()
    vcf_gz = args.vcf_gz
    min_depth = int(args.min_depth)
    max_strand_prop = float(args.max_strand_prop)
    bcftools_path = args.bcftools_path
    sample = args.sample_id
    reference = args.reference_genome
    bam_input = args.bam_input
    output_fasta = args.output_fasta
    coverage_output = args.coverage_output
    snps_to_exclude = args.masked_bed
    print(sample)
    get_sample_vcf(vcf_gz=vcf_gz, sample=sample)
    get_masking_bed(bam_input, min_depth, max_strand_prop) 
    get_consensus_fasta(reference,sample, output_fasta, coverage_output, snps_to_exclude)
if __name__=="__main__":
    main()
