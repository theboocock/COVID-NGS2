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

bioawk_cmd="""bioawk -c fastx '{{if($comment ~ "{taxon_to_extract}") {{ print "@"$name" "$comment"\\n"$seq"\\n+\\n"$qual}}}}' {read_in} | bgzip -c > {read_out} """ 
bioawk_deplete_cmd ="""bioawk -c fastx '{{if($comment !~ "{taxon_to_deplete}") {{ print "@"$name" "$comment"\\n"$seq"\\n+\\n"$qual}}}}' {read_in} | bgzip -c > {read_out} """ 
human_taxid = "kraken:taxid\\\\|9606"

def filter_reads(read_one_in, read_two_in, read_one_out, read_two_out, taxon_to_extract, deplete): 
    """
        Filter reads
    """
    
    if deplete: 
        cmd_one = bioawk_deplete_cmd.format(taxon_to_deplete=human_taxid, read_in=read_one_in, read_out=read_one_out)
        cmd_two = bioawk_deplete_cmd.format(taxon_to_deplete=human_taxid, read_in=read_two_in, read_out=read_two_out)
        subprocess.check_call(cmd_one, shell=True)
        subprocess.check_call(cmd_two, shell=True)
    else:
        search_string = "kraken:taxid\\\\|{taxon_id}".format(taxon_id=taxon_to_extract)
        cmd_one = bioawk_cmd.format(taxon_to_extract=search_string, read_in=read_one_in, read_out=read_one_out)
        cmd_two = bioawk_cmd.format(taxon_to_extract=search_string, read_in=read_two_in, read_out=read_two_out)
        subprocess.check_call(cmd_one, shell=True)
        subprocess.check_call(cmd_two, shell=True)

def main():
    parser=argparse.ArgumentParser(description="Filter zipped kraken fastqs to only extract a read unique assigned to a specific taxon")
    parser.add_argument("--taxon-to-extract", dest="taxon_to_extract")
    parser.add_argument("--read-one-in",dest="read_one_in")
    parser.add_argument("--read-two-in",dest="read_two_in")
    parser.add_argument("--read-one-out",dest="read_one_out")
    parser.add_argument("--read-two-out",dest="read_two_out")
    parser.add_argument("--deplete", dest="deplete",action="store_true", default=False)
    args = parser.parse_args()
    filter_reads(args.read_one_in,args.read_two_in, args.read_one_out, args.read_two_out, args.taxon_to_extract,args.deplete)

if __name__=="__main__":
    main()

