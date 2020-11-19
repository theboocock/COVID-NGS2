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
import pandas as pd
import numpy as np
sheet_header="""
[Header],,,,,,,,,
IEMFileVersion,4,,,,,,,,
Investigator_Name,HJ,,,,,,,,
Experiment_Name,{name},,,,,,,,
Date,{date},,,,,,,,
Workflow,GenerateFASTQ,,,,,,,,
Application,HiSeq_FASTQ_Only,,,,,,,,
Assay,TruSeq_LT,,,,,,,,
Description,,,,,,,,,
Chemistry,Default,,,,,,,,
,,,,,,,,,
[Reads],,,,,,,,,
150,,,,,,,,,
,,,,,,,,,
[Settings],,,,,,,,,
ReverseComplement,0,,,,,,,,
Adapter,AGATCGGAAGAGCACACGTCTGAACTCCAGTCA,,,,,,,,
AdapterRead2,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT,,,,,,,,
,,,,,,,,,
[Data],,,,,,,,,
"""
__bcl_to_fastq__="""
    /u/scratch/s/smilefre/covidngs/COVID19-NGS/software/bcl2fastq/bcl2fastq  --runfolder-dir . -o  {run_folder} --reports-dir report/ --sample-sheet {sample_sheet} --ignore-missing-bcl --barcode-mismatches {barcode_mismatches} --no-lane-splitting --minimum-trimmed-read-length 0 --mask-short-adapter-reads 0
"""
def make_sample_sheet(sample_sheet,  experiment_name,date, sample_sheet_neb, sample_sheet_amp):
    """
        Experiment name 
    """
    sheet_line = sheet_header.format(name=experiment_name,date=date)
    print("HERERER")
    print(sample_sheet)
    pd_in = pd.read_csv(sample_sheet,sep=",")
    print(pd_in["index"])
    amp = pd_in["sample"].str.contains("Ms")
    neb = pd_in["sample"].str.contains("NEB")
    pd_in["Lane"] = 1
    pd_in["uid"] = pd_in["uid"].apply(np.char.upper)
    print(pd_in)
    pd_in["index"] = pd_in["index"].apply(np.char.upper)
    pd_in["index2"] = pd_in["index2"].apply(np.char.upper)
    pd_in["Sample_ID"]  = pd_in["sample"] + "_" + pd_in["uid"]
    pd_in["library_type"] = "amp_complicated"
    pd_in["sample_type"] = "amp"
    pd_in["uid_sample_type"] = pd_in["uid"] + "_" + pd_in["sample_type"]
    pd_in["library_type"][pd_in["sample"].str.contains("NEB")] = "neb"
    if "oligos_used" not in pd_in.columns:
        # Hack to get the columns right
        pd_in["oligos_used"] = "3,4"
    neb_df  = pd_in[neb] 
    amp_df =  pd_in[amp]
    amp_df_first = amp_df.copy()
    neb_df_first = neb_df.copy()
    pd_in["sample"] = list(range(1,pd_in.shape[0]+1))
    for i in range(2,5):
        tmp_row = neb_df_first.copy()
        tmp_row["Lane"] = i
        neb_df= neb_df.append(tmp_row)
        tmp_row = amp_df_first.copy()
        tmp_row["Lane"] = i
        amp_df = amp_df.append(tmp_row) 
    with open(sample_sheet_neb, "w") as out_f:
        out_f.write(sheet_line)
        neb_df.to_csv(out_f, mode="a",index=False)
    with open(sample_sheet_amp, "w") as out_f:
        out_f.write(sheet_line)
        amp_df.to_csv(out_f, mode="a",index=False)
    pd_in.to_csv("sample_sheet.csv",index=False, sep="\t")


import subprocess
def run_bcl2fastq(sample_sheet_neb, sample_sheet_amp):
    """
        run bcl2fastq
    """
    neb_cmd = __bcl_to_fastq__.format(run_folder=".",sample_sheet=sample_sheet_neb, barcode_mismatches=1)
    subprocess.check_call(neb_cmd, shell=True)
    amp_cmd = __bcl_to_fastq__.format(run_folder=".",sample_sheet=sample_sheet_amp, barcode_mismatches=0)
    #subprocess.check_call(amp_cmd, shell=True)

def main():

    parser = argparse.ArgumentParser(description="Create sample sheet")
    parser.add_argument("-s","--sheet", dest="sample_sheet", help="Pre sample sheet")
    parser.add_argument("-e","--experiment-name",dest="exp_name", help="Experiment name", required=True)
    parser.add_argument("-d","--date",dest="date", help="experiment date", required=True)
    parser.add_argument("--output-file-neb",dest="output_file_neb", help="Output_file_neb", required=True)
    parser.add_argument("--output-file-amp",dest="output_file_amp", help="Output_file_amp", required=True)
    args = parser.parse_args()
    make_sample_sheet(args.sample_sheet, args.exp_name, args.date, args.output_file_neb, args.output_file_amp)
    run_bcl2fastq(args.output_file_neb, args.output_file_amp)

if __name__=="__main__":
    main()
