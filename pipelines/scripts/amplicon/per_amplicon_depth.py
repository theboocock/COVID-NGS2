#!/usr/bin/env python

import argparse
import pysam


def count_read_per_amplicon(bam_input, quality_filter):

    count_hash= {}
    bam_reads = pysam.AlignmentFile(bam_input, "rb")
    for read in bam_reads:
        if read.mapping_quality > quality_filter:
            read_name = (read.query_name).split(":")
            amplicon = read_name[0]
            umi = read_name[2]
            try:
                count_hash[amplicon] += 1
            except:
                count_hash[amplicon] = 1
    return count_hash

def output_dataframe(count_hash, oligo_pool, output_df, library_type):

    with open(oligo_pool) as oligo_in:
        with open(output_df,"w") as out_df:
            out_df.write("primer_name\tcount\tlibrary_type\n")
            for line in oligo_in:
                split_line = line.strip().split("\t")
                primer_name = split_line[3]
                try:
                    count = count_hash[primer_name]
                except:
                    count = 0
                out_df.write(primer_name + "\t" + str(count) + "\t" + library_type+ "\n")


def output_complicated_dataframe(count_hash, oligo_pool, output_df, oligos_used):
    
    with open(oligo_pool) as oligo_in:
        with open(output_df,"w") as out_df:
            out_df.write("primer_name\tcount\toligo_set\toligos_used\n")
            for line in oligo_in:
                split_line = line.strip().split("\t")
                primer_name = split_line[0]
                set_of_oligos = split_line[2]
                try:
                    count = count_hash[primer_name]
                except:
                    count = 0
                out_df.write(primer_name + "\t" + str(count) + "\t"  + set_of_oligos + "\t"+ oligos_used+ "\n")



def main():
    parser = argparse.ArgumentParser(description="Per get read depth per amplicon")
    parser.add_argument("-b","--bam-input", dest="bam_input", help="BAM input",required=True)
    parser.add_argument("-p","--oligo-pool", dest="oligo_pool", help="Oligo_pool", required=True)
    parser.add_argument("-q","--quality-filter",dest="mapping_quality",help="Amplicon depth", required=True)
    parser.add_argument("-o","--output-dataframe",dest="output_df", help="output df", required=True)
    parser.add_argument("-s","--sample-type",dest="sample_type", help="sample type", required=True) 
    parser.add_argument("-c","--complicated",dest="complicated", help="Complicated run", default=False, action="store_true")
    parser.add_argument("--oligos-used", dest="oligos_used", help="Oligos used")
    args = parser.parse_args()
    bam_input = args.bam_input
    oligo_pool = args.oligo_pool
    output_df = args.output_df
#    complicated = args.complicated
    min_mapping_quality = int(args.mapping_quality)
    count_hash = count_read_per_amplicon(bam_input,min_mapping_quality)
    if not args.complicated:
        output_dataframe(count_hash, oligo_pool, output_df, args.sample_type) 
    else:
        output_complicated_dataframe(count_hash,oligo_pool, output_df, args.oligos_used)

if __name__ == "__main__":
    main()
