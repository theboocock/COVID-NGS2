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
import os
import pysam

def vcf_input_filter(vcf_file,lower, upper, vcf_outfile):
    """

    """
    vcf_in = pysam.VariantFile(vcf_file)
    with open(vcf_outfile,"w") as out_f:
        vcf_in.header.add_line("""##FORMAT=<ID=AR,Number=1,Type=Float,Description="Allelic ratio">""")
        out_f.write(str(vcf_in.header))
        for rec in vcf_in.fetch():
            out_f.write(str(rec.chrom)+ "\t")
            out_f.write(str(rec.pos) + "\t")
            identification = rec.id
            if rec.id == None:
                record_id = "."
            else:
                record_id = rec.id

            out_f.write(str(record_id) + "\t")
            out_f.write(str(rec.ref) + "\t")
            out_f.write(str(",".join(rec.alts)) + "\t")
            if rec.qual is None:
                record_qual = "."
            else:
                record_qual = rec.qual
            out_f.write(str(int(record_qual)) + "\t")
            if rec.filter is None:
                record_filter = "." 
            else:
                filter_keys_len = len(rec.filter.keys())
                if filter_keys_len  == 0:
                    record_filter = "." 
                else:
                    record_filter = rec.filter.keys()[0]
            out_f.write(str(record_filter) + "\t")
            info_keys=(list(rec.info))
            info_str = ";".join([key  +"=" + str(rec.info[key]) for key in info_keys])
            out_f.write(str(info_str) + "\t")
            out_f.write(str(":".join(rec.format)) + ":AR" + "\t")
            for sample_idx, keys in enumerate(rec.samples.items()):
                key = keys[0]
                sample = keys[1]
                allelic_depth = (sample["AD"])
                dp = (sample["DP"])
                ref_allele = sample["AD"][0]
                if ref_allele != None:
                    alt_allele = 0
                    other = ""
                    for i, alt in enumerate(sample["AD"][1:]):
                        if alt != None:
                            alt_allele = alt 
                            idx = i
                    if ref_allele == 0 and alt_allele == 0:
                        ratio = 0.0
                    else:
                        ratio  = alt_allele/(alt_allele  + ref_allele)
                    gt = sample["GT"] 
                    if int(dp) < 3:
                        gt = "./."
                    else:
                        if ratio < lower:
                            gt = "0/0"
                        elif ratio > upper:
                            gt = str(sample["GT"][1]) + "/" + str(sample["GT"][1])
                        else:
                            other = str(idx+1)
                            gt = "0" + "/" + other 
                            if int(dp) < 10:
                                gt = './.'
                                gt = './.'

                    keys = sample.keys()[1:]
                    keys_out = [sample[key] for key in keys]
                    out_str_list = [] 
                    for item in keys_out:
                        try:
                            out_str_list.append(",".join([str(it) if it is not None else "." for it in item ]))
                        except:
                            out_str_list.append(str(item))
                    out_f.write(gt +":" +":".join(out_str_list))
                    out_f.write(":"+str(ratio))
                else:
                    gt = sample["GT"] 
                    # TODO: make sure to extract the relevant information from the reference
                    if gt[0] == None:
                        gt = "./."
                        ratio = "."
                    else:
                        gt = "0/0"
                        ratio = "0"
                    keys = sample.keys()[1:]
                    out_f.write(gt +":.:.:.:"+ratio)
                if sample_idx != (len(rec.samples.keys()) -1):
                    out_f.write("\t")
            out_f.write("\n")


def main():
    parser = argparse.ArgumentParser(description="Allelic balance filter")
    parser.add_argument("--vcf",dest="vcf_input", help="VCF input")
    parser.add_argument("-l","--lower", dest="lower_threshold", default=.2)
    parser.add_argument("-u","--upper", dest="upper_threshold", default=.8)
    parser.add_argument("-o","--vcf-out", dest="vcf_output") 
    args = parser.parse_args()
    lower = float(args.lower_threshold)
    upper = float(args.upper_threshold)
    vcf_outfile = args.vcf_output 
    vcf_file = args.vcf_input
    vcf_input_filter(vcf_file, lower, upper, vcf_outfile)

if __name__ == "__main__":
    main()
