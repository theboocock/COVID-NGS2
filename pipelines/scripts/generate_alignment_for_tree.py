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
from Bio import SeqIO
import codecs 

### Coverage filter will let us remove genomes that are dogshit

def combine_fasta(query_fasta, clades_in, sub_sampled_tree, coverage_filter=0.9):
    seq_hash = {}
    for seq in SeqIO.parse(clades_in,"fasta"):
        decrypted_seq = codecs.decode(str(seq.seq), "rot_13")
        if "outgroup" in seq.id:
            seq_id = seq.id
        else:
            seq_id_split = seq.id.split("_")
            seq_id = seq_id_split[1]
        if sub_sampled_tree:
            seq_hash[seq_id] = decrypted_seq
        else:
            seq_hash[seq.id] = decrypted_seq 
    for seq in SeqIO.parse(query_fasta,"fasta"):
        seq_hash[seq.id] = str(seq.seq)
    with open("outputs.fasta","w") as out_f:
        for key, value in seq_hash.items():
            out_f.write(">" + key +"\n" + value + "\n")
            print(key)

MUSCLE_CMD="""mafft  outputs.fasta > outputs.align.fasta"""
RAXML_CMD="""../../../ngs/software/raxml-ng/bin/raxml-ng  --msa outputs.align.fasta --model GTR+G --prefix tree --threads 1 --seed  2"""
def align_sequneces():
    subprocess.check_call(MUSCLE_CMD,shell=True)
    subprocess.check_call(RAXML_CMD,shell=True)
def main():
    parser = argparse.ArgumentParser(description="Take query and database sequences and generate alignments")
    parser.add_argument("-q","--query-fasta",dest="query_fasta",help="Query fasta input")
    parser.add_argument("-c","--clades-in",dest="clades_in", help="Clades in")
    parser.add_argument("-s","--sub-sampled-tree",dest="sub_sampled_tree",action="store_true",default=False)
    args =parser.parse_args()
    query_fasta = args.query_fasta
    clades_in = args.clades_in
    sub_sampled_tree = args.sub_sampled_tree
    combine_fasta(query_fasta, clades_in, sub_sampled_tree)
    align_sequneces()
if __name__== "__main__":
    main()

