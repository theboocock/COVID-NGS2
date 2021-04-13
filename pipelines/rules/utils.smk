### Starting index = 0
#def get_uid_library_type(wildcards):
#    uid = metadata["uid_library_type"][metadata["uid_library_type"] == wildcards.sample]
#    return(uid)

#def getuid_sample_map_r1(wildcards):
#    r1 = metadata["read_one"][metadata["uid_library_type"] == library_uid]
#    return(r1)


def get_mean_coverage(wildcards):
    print(wildcards.sample)
    cov = metdata_dedup[metdata_dedup["merged_id"] == wildcards.sample]
    return cov["coverage_mean_dedup"] 

def get_pango_in_all(wildcards):
    return(f"outputs/consensus/merged/phased/sars2/all_{MIN_DEPTH}_none.fasta")
def get_sample_depth(wildcards):
    return wildcards.sample +  "_"  +  str(wildcards.down)

def get_bams(wildcards):
    samples = mapped_uid[wildcards.sample]
    return(expand("outputs/md/{intervals}/{sample}.sorted.md.bam", intervals=intervals, sample= samples))

def get_bai(wildcards):
    samples = mapped_uid[wildcards.sample]
    return(expand("outputs/md/{intervals}/{sample}.sorted.md.bam.bai", intervals=intervals, sample= samples))
 
def get_sample_name_for_merge(wildcards):
    sample_id = sample_names_hash[wildcards.sample][0]
    return(sample_id)
#TODO: merged bams
#rule merged_bams:
#    shadow: "minimal"
#    group: "align"
#    input:
#        bams = get_bams, 
#        bai = get_bai
 #   output:
##        bam="outputs/md/merged/{intervals}/{sample}.bam",
#        bai="outputs/md/merged/{intervals}/{sample}.bam.bai"
#    params:
#        rg=get_sample_name_for_merge    
#    run:
#        if len(input.bams) != 1:
#            # Skip merge
#            shell("samtools merge tmp.bam {input.bams} && samtools sort tmp.bam > sorted.bam")
#            shell("java -jar {PICARD_PATH} AddOrReplaceReadGroups I=sorted.bam O={output.bam} RGID=1 RGSM={params.rg} RGLB=4 RGPL=ILLUMINA RGPU=unit1")
#        else:
#            shell("java -jar {PICARD_PATH} AddOrReplaceReadGroups I={input.bams} O={output.bam} RGID=1 RGSM={params.rg} RGLB=4 RGPL=ILLUMINA RGPU=unit1")
#        shell("samtools index {output}")
#

def get_rg(wildcards):
    strain_name =get_sample_name(wildcards)
    rg="@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:lib1\\tPU:unit1".format(sample=strain_name)
    return(rg)
def uid_library_to_sample_name(): 
    id_library_type  = metadata["uid"] 
    metadata["id_library_type"] = id_library_type
    unique_id_library_type = list(pd.unique(metadata["id_library_type"]))
    mapped_uids_library = {} 
    sample_names = {}
    mapped_uids_read_one = {}
    mapped_uids_read_two = {}
    start_id = 1  
    #print(len(unique_id_library_type))
    for uid_name in unique_id_library_type: 
        g = list((metadata["sample"][metadata["id_library_type"] == uid_name]))
        r1 =  list((metadata["fastq_one"][metadata["id_library_type"] == uid_name]))
        r2 =  list((metadata["fastq_two"][metadata["id_library_type"] == uid_name]))
        read_one = list((metadata["sample"][metadata["id_library_type"] == uid_name]))
        read_two = list((metadata["sample"][metadata["id_library_type"] == uid_name]))
        sample_names["LA_{:04d}".format(start_id)] = ["LA-{:04d}".format(start_id), uid_name]
        mapped_uids_library["LA_{:04d}".format(start_id)]= g
        mapped_uids_read_one["LA_{:04d}".format(start_id)] = r1
        mapped_uids_read_two["LA_{:04d}".format(start_id)] = r2
        start_id += 1
    return(mapped_uids_library, sample_names, mapped_uids_read_one, mapped_uids_read_two)


def uid_library_to_sample_hash(merge_rule="uid_type"):
    if merge_rule == "uid_type":
        id_library_type = metadata["uid"] + "_" + metadata["sample_type"] 
    elif merge_rule == "uid":
        id_library_type = metadata["uid"]
    elif merge_rule == "uid_oligo":
        id_library_type = metadata["uid"]  + "_" + metadata["oligos_used"]
    elif merge_rule == "column_uid_sample_type":
        #   print("HERE")
        id_library_type = metadata["uid_sample_type"]
    elif merge_rule == "per_sample":
        id_library_type = metadata["sample"]

    metadata["id_library_type"] = id_library_type
    unique_id_library_type = list(pd.unique(metadata["id_library_type"]))
    mapped_uids_library = {} 
    sample_names = {}
    mapped_uids_read_one = {}
    mapped_uids_read_two = {}
    start_id = 1 
    for uid_name in unique_id_library_type: 
        g = list((metadata["sample"][metadata["id_library_type"] == uid_name]))
        r1 =  list((metadata["fastq_one"][metadata["id_library_type"] == uid_name]))
        r2 =  list((metadata["fastq_two"][metadata["id_library_type"] == uid_name]))
        read_one = list((metadata["sample"][metadata["id_library_type"] == uid_name]))
        read_two = list((metadata["sample"][metadata["id_library_type"] == uid_name]))
        sample_names["LA_{:04d}".format(start_id)] = ["LA-{:04d}".format(start_id), uid_name]
        mapped_uids_library["LA_{:04d}".format(start_id)]= g
        mapped_uids_read_one["LA_{:04d}".format(start_id)] = r1
        mapped_uids_read_two["LA_{:04d}".format(start_id)] = r2
        start_id += 1
    return(mapped_uids_library, sample_names, mapped_uids_read_one, mapped_uids_read_two)


def uid_library_to_sample_name_and_readgroup(): 
    id_library_type  = metadata["uid"] + "_" + metadata["sample_type"] 
    metadata["id_library_type"] = id_library_type
    unique_id_library_type = list(pd.unique(metadata["id_library_type"]))
    mapped_uids_library = {} 
    sample_names = {}
    mapped_uids_read_one = {}
    mapped_uids_read_two = {}
    start_id = 1 
    for uid_name in unique_id_library_type: 
        g = list((metadata["sample"][metadata["id_library_type"] == uid_name]))
        r1 =  list((metadata["fastq_one"][metadata["id_library_type"] == uid_name]))
        r2 =  list((metadata["fastq_two"][metadata["id_library_type"] == uid_name]))
        read_one = list((metadata["sample"][metadata["id_library_type"] == uid_name]))
        read_two = list((metadata["sample"][metadata["id_library_type"] == uid_name]))
        sample_names["LA_{:04d}".format(start_id)] = ["LA-{:04d}".format(start_id), uid_name]
        mapped_uids_library["LA_{:04d}".format(start_id)]= g
        mapped_uids_read_one["LA_{:04d}".format(start_id)] = r1
        mapped_uids_read_two["LA_{:04d}".format(start_id)] = r2
        start_id += 1
    return(mapped_uids_library, sample_names, mapped_uids_read_one, mapped_uids_read_two)

def getuid_sample_map_r2(wildcards):
    r2 = metadata["read_two"][metadata["uid_library_type"] == library_uidi]
    return(r2)

def get_sample_type(wildcards):
    sample_type = metadata["sample_type"][metadata["sample"] == int(wildcards.sample)]
    sample_type = str(sample_type).split()[1]
    return(sample_type)

def get_oligos_used(wildcards):
    if "oligos_used" in metadata.columns:
        oligo_used = metadata["oligos_used"][metadata["sample"] == int(wildcards.sample)]
        oligo_used = str(oligo_used).split()[1]
        return oligo_used 
    else:
        return "NA"

def get_fq1(wildcards):
    out_fastq = metadata["fastq_one"][metadata["sample"] == int(wildcards.sample)]
    return out_fastq 

def get_fq2(wildcards):
    out_fastq = metadata["fastq_two"][metadata["sample"] == int(wildcards.sample)]
    return out_fastq

def get_sample_name(wildcards):
    sample_name = metadata["sample"][metadata["sample"] == int(wildcards.sample)]
    sample_name = str(sample_name).split()[1]
    return sample_name

def read_taxid(taxid_file):
    taxid_dict = {}
    with open(taxid_file) as out_f:
        for line in out_f:
            refseq = line.strip().split()[0]
            taxid = line.strip().split()[1]
            taxid_dict[refseq] = taxid 
    return(taxid_dict)

import pandas as pd
def fastq_match(config):
    try:
        config["sample_list"]
    except KeyError:
        sys.stderr.write("Please specifiy the sample input file --sample in on the command line")
        sys.exit(1)
    meta_in=pd.read_csv(config["sample_list"],sep="\t")
    return meta_in 

