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


## Process picard insert sizes ## 
def process_insert_sizes(insert_sizes):
    with open(insert_sizes) as in_insert:
        next_line = False 
        for i, line in enumerate(in_insert):
            if next_line: 
                split_line= line.split()
                return(split_line[0] + " " + split_line[2])
            if line.startswith("MEDIAN_INSERT_SIZE"):
                next_line = True 

def combine_qc(map_stat_in,coverage_in, insert_size_line,output_file, complete_coverage_in,read_count, mapstat_md,no_dedup_cov):
    with open(output_file,"w") as out_f:
        with open(map_stat_in) as map_in:
            with open(coverage_in) as coverage_f:
                with open(complete_coverage_in) as coverage_full:
                    with open(read_count) as read_count_in:
                        with open(mapstat_md) as map_stat_md_sars:
                            with open(no_dedup_cov) as no_dedup_cov_out:
                                cov_total = 0
                                for i, line in enumerate(coverage_full):
                                    s_line = line.split("\t")
                                    cov = int(s_line[2])
                                    cov_total += cov
                                cov_mean = str(cov_total / ( i * 1.0))
                                cov_total = 0
                                for i, line in enumerate(no_dedup_cov_out):
                                    s_line = line.split("\t")
                                    cov = int(s_line[2])
                                    cov_total += cov
                                cov_mean_no_dedup = str(cov_total / ( i * 1.0))
                                coverage = coverage_f.readline().strip()
                                map_stat_in = map_in.readline().strip()
                                map_stat_in = map_stat_in.replace("\t"," ")
                                read_count = str(int(read_count_in.readline().strip())/2)
                                mapstat_md_sars_number= str(int(map_stat_md_sars.readline().strip()))
                                out_f.write(coverage + " "  + cov_mean+ " " + cov_mean_no_dedup + " " + map_stat_in + " " + read_count + " "+ 
                                        mapstat_md_sars_number  + " " + insert_size_line+ "\n")



def main():
    parser = argparse.ArgumentParser(description="Aggregate QC outputs for COVID19 NGS pipeline")
    parser.add_argument("-m","--map-stat-in",dest="map_stat_in", help="Mapping statistics in", required=True)
    parser.add_argument("-c","--coverage-in",dest="coverage_in", help="Coverage input",required=True)
    parser.add_argument("-a","--complete-coverage-in",dest="complete_coverage_in", help="Coverage input",required=True)
    parser.add_argument("-i","--insert-sizes",dest="insert_sizes",help="Insert sizes",required=True)
    parser.add_argument("-o","--output-file",dest="output_file",help="output_file",required=True)
    parser.add_argument("--read-count",dest="read_count_in",help="Read count", required=True)
    parser.add_argument("--mapstat-md",dest="mapstat_md",help="Map statmp sars2 mapping count",required=True)
    parser.add_argument("--no-dedup-coverage",dest="no_dedup_coverage",help="No deduplications coverae total", required=True)
    args = parser.parse_args()
    map_stat_in = args.map_stat_in
    coverage_in = args.coverage_in
    insert_sizes = args.insert_sizes
    output_file = args.output_file
    complete_coverage = args.complete_coverage_in
    read_count = args.read_count_in
    insert_size_line = process_insert_sizes(insert_sizes)
    mapstat_md = args.mapstat_md
    no_dedup_cov = args.no_dedup_coverage
    if insert_size_line == None:
        insert_size_line = "0 0"
    combine_qc(map_stat_in, coverage_in, insert_size_line, output_file, complete_coverage, read_count, mapstat_md,no_dedup_cov) 



if __name__ == "__main__":
    main()
