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
import gzip

human_taxid="9606"
sars_taxid="2697049"


class centrifuge_info:

    def __init__(self, read_id, score, taxid, genome_name):
        self.read_id = read_id
        self.score = score
        self.taxid= taxid 
        self.genome_name = genome_name
    def get_read_id(self):
        return self.read_id
    def get_score(self):
        return self.score
    def get_taxid(self):
        return self.taxid
    def get_genome_name(self):
        return self.genome_name
    def __str__(self):
        return self.read_id +" " + self.score + " " + self.taxid

def get_info(line_s):
    
    read_in = line_s[0]
    score = line_s[3]
    taxid = line_s[2]
    genome_name = line_s[1]
    info = centrifuge_info(read_in, score, taxid, genome_name)
    return(info)

def write_lines(read_one_f, read_two_f, read_one_sars2, read_two_sars2):
        line_one =read_one_f.readline().decode()
        read_one_sars2.write(line_one)
        line_two =read_two_f.readline().decode()
        read_two_sars2.write(line_two)
        line_one =read_one_f.readline().decode()
        read_one_sars2.write(line_one)
        line_two =read_two_f.readline().decode()
        read_two_sars2.write(line_two)
        line_one =read_one_f.readline().decode()
        read_one_sars2.write(line_one)
        line_two =read_two_f.readline().decode()
        read_two_sars2.write(line_two)
        line_one =read_one_f.readline().decode()
        read_one_sars2.write(line_one)
        line_two =read_two_f.readline().decode()
        read_two_sars2.write(line_two)

def skip_lines(read_one_f, read_two_f):
        line_one =read_one_f.readline().decode()
        line_two =read_two_f.readline().decode()
        line_one =read_one_f.readline().decode()
        line_two =read_two_f.readline().decode()
        line_one =read_one_f.readline().decode()
        line_two =read_two_f.readline().decode()
        line_one =read_one_f.readline().decode()
        line_two =read_two_f.readline().decode()


def extract_reads(read_one_f,read_two_f, read_match,count_hash, sars2_out_read_one, sars2_out_read_two):

    ### count hash ##
    
    if read_match.get_taxid() == sars_taxid:
        count_hash[sars_taxid]  +=1
        write_lines(read_one_f, read_two_f, sars2_out_read_one,sars2_out_read_two) 
    elif read_match.get_taxid() == human_taxid:
        count_hash[human_taxid]  +=1
        skip_lines(read_one_f, read_two_f)
    elif read_match.get_taxid() == "0":
        count_hash["unclassified"] += 1
        skip_lines(read_one_f, read_two_f)
    else:
        count_hash["unknown"] += 1
        skip_lines(read_one_f, read_two_f)
    return(count_hash)

def filter_reads(read_one, read_two, centrifuge_output, read_one_out, read_two_out):
    count_hash = {}
    count_hash[sars_taxid] = 0
    count_hash[human_taxid] = 0
    count_hash["unclassified"]=0 
    count_hash["unknown"]=0 
    with gzip.open(read_one) as read_one_f:
        with  gzip.open(read_two) as read_two_f:
            with open(centrifuge_output) as centrifuge_output_f:
                line = centrifuge_output_f.readline()
                line = centrifuge_output_f.readline()
                line_s = line.split("\t")
                current_read = line_s[0]
                read_match = get_info(line_s)      
                while line: 
                    line_s = line.split("\t")
                    new_read = line_s[0]
                    if new_read != current_read:
                        current_read = new_read
                        count_hash = extract_reads(read_one_f, read_two_f, read_match, count_hash, read_one_out, read_two_out) 
                        read_match = get_info(line_s)
                    line = centrifuge_output_f.readline()

#                    for i, line in enumerate(centrifuge_output_f):
#                        if i > 0:
#                            line_split = line.split("\t")
#                            read_one_line = read_one_f.readline()
    return(count_hash)

def write_count_hash(count_report, count_read_out_report):
    with open(count_read_out_report,"w") as count_out:
        for key, value in count_report.items():
            count_out.write(key + " " + str(value) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Filter reads from centrifuge run") 
    parser.add_argument("--read-one",dest="read_one",help="read_one")
    parser.add_argument("--read-two",dest="read_two",help="read_two")
    parser.add_argument("--read-one-out", dest="read_one_out",help="read one out")
    parser.add_argument("--read-two-out", dest="read_two_out", help="read two out")
    parser.add_argument("--centrifuge-report",dest="centrifuge_output",help="centrifuge output file")
    parser.add_argument("--count-out-table", dest="count_read_out_report",help="count report")
    args = parser.parse_args()
    read_one = args.read_one
    read_two = args.read_two
    read_one_out = open(args.read_one_out,"w")
    read_two_out = open(args.read_two_out,"w")
    centrifuge_output = args.centrifuge_output
    count_report = filter_reads(read_one,read_two, centrifuge_output, read_one_out, read_two_out)
    write_count_hash(count_report, args.count_read_out_report)
    read_one_out.close()
    read_two_out.close()
if __name__=="__main__":
    main()
