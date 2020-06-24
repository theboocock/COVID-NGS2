
KRAKEN2_PATH=config["PATHS"]["KRAKEN2_PATH"]
KRAKEN2_DB=config["PATHS"]["KRAKEN2_DB"]

import os



rule count_reads:
    input:
        read_one=get_fq1,
    output:
        "outputs/read_counts/{sample}.txt"
    shell:
        "zcat {input} | wc -l > {output}"

rule trim_reads:
    conda: workflow.basedir + "/envs/r-covid.yaml"
    input:
        read_one=get_fq1,
        read_two=get_fq2
    output:
        out_fastq_one="outputs/fastqs/{sample}.R1.fastq.gz",
        out_fastq_two="outputs/fastqs/{sample}.R2.fastq.gz"
    params:
        sample_type =  get_sample_type
    log:
        trim_read_out="log/fastqs/{sample}.log"
    shell:
        "{SCRIPTS_DIR}/amplicon/run_custom_amplicon.sh {params.sample_type} {input.read_one} {input.read_two} {output.out_fastq_one} {output.out_fastq_two} {OLIGO_POOL_BED} {ARCTIC_FILE} {log}"

rule kraken_assignment:
    group: "align"
    input:
        read_one="outputs/fastqs/{sample}.R1.fastq.gz",
        read_two="outputs/fastqs/{sample}.R2.fastq.gz"
    output:
        read_one = temp("outputs/kraken2/{sample}.cseqs__1.fastq"),
        read_two = temp("outputs/kraken2/{sample}.cseqs__2.fastq"),
        unclassified_one = temp("outputs/kraken2/{sample}.ucseqs__1.fastq"),
        unclassified_two = temp("outputs/kraken2/{sample}.ucseqs__2.fastq"),
        kraken_report = "outputs/kraken2/{sample}.kraken_report.txt",
        kraken_output = temp("outputs/kraken2/{sample}_output.txt")
    resources:
        mem_mb=4000,
        walltime=43200
    threads: 4
    params:
        out_prefix= "outputs/kraken2/{sample}",
        taxon_to_extract= "9606"
    log:
        "logs/kraken2/{sample}.log"
    run:
        try:
            shell("touch {output.read_one}")
            shell("touch {output.read_two}")
            shell("touch {output.unclassified_one}")
            shell("touch {output.unclassified_two}")
            shell("touch {output.kraken_report}")
            shell("touch {output.kraken_output}")
            run_cmd = "{SCRIPTS_DIR}/kraken_run_and_filter.py --kraken-path {KRAKEN2_PATH} --kraken-db {KRAKEN2_DB} --fastq-in-one  {input_read_one} --fastq-in-two {input_read_two} --fastq-out-prefix {params_out_prefix}  --fastq-out-one {output_read_one} --fastq-out-two {output_read_two} --kraken-report {output_kraken_report} --kraken-output {output_kraken_output}  --threads {threads} 2> {log}".format(SCRIPTS_DIR=SCRIPTS_DIR,KRAKEN2_PATH=KRAKEN2_PATH, KRAKEN2_DB=KRAKEN2_DB, input_read_one=input.read_one, input_read_two=input.read_two, params_out_prefix=params.out_prefix,output_read_one=output.read_one, output_read_two=output.read_two,output_kraken_report=output.kraken_report, output_kraken_output=output.kraken_output,threads=threads,log=log)
            subprocess.check_call(run_cmd,shell=True)
        except subprocess.CalledProcessError as exc:
            shell("touch {output.read_one}")
            shell("touch {output.read_two}")
            shell("touch {output.unclassified_one}")
            shell("touch {output.unclassified_two}")
            shell("touch {output.kraken_report}")
            shell("touch {output.kraken_output}")
            pass

   
def get_all_uid_read_one(wildcards):
    samples = expand("outputs/fastqs/{sample}.R1.fastq.gz", sample=mapped_uid[wildcards.sample])
    return(samples)

def get_all_uid_read_two(wildcards):
    samples = expand("outputs/fastqs/{sample}.R2.fastq.gz", sample=mapped_uid[wildcards.sample])
    return(samples)
   
def get_all_uid_kraken_report(wildcards):
    samples = expand("outputs/kraken2/{sample}.kraken_report.txt" , sample=mapped_uid[wildcards.sample])
    return(samples)

rule kraken_assignment_merge:
    group: "align"
    shadow: "minimal"
    input:
        read_one=get_all_uid_read_one,
        read_two=get_all_uid_read_two,
        kraken_inputs= get_all_uid_kraken_report
    output:
        read_one = temp("outputs/kraken2/merged/{sample}.cseqs__1.fastq"),
        read_two = temp("outputs/kraken2/merged/{sample}.cseqs__2.fastq"),
        unclassified_one = temp("outputs/merged/kraken2/{sample}.ucseqs__1.fastq"),
        unclassified_two = temp("outputs/merged/kraken2/{sample}.ucseqs__2.fastq"),
        kraken_report = "outputs/kraken2/merged/{sample}.kraken_report.txt",
        kraken_output = temp("outputs/kraken2/merged/{sample}_output.txt")
    resources:
        mem_mb=4000,
        walltime=43200
    threads: 4
    params:
        out_prefix= "outputs/kraken2/merge/{sample}",
        taxon_to_extract= "9606"
    log:
        "logs/kraken2/merge/{sample}.log"
    run:
        try:
            shell("touch {output.read_one}")
            shell("touch {output.read_two}")
            shell("touch {output.unclassified_one}")
            shell("touch {output.unclassified_two}")
            shell("touch {output.kraken_report}")
            shell("touch {output.kraken_output}")
            shell("touch {log}")
            if(len(input.kraken_inputs) == 1):
                import shutil
                shutil.copy(input.kraken_inputs[0], output.kraken_report)
            else:
                for r1,r2 in zip(input.read_one, input.read_two):
                    shell("cat {r1} >> r1.fastq.gz")
                    shell("cat {r2} >> r2.fastq.gz")
                run_cmd = "{SCRIPTS_DIR}/kraken_run_and_filter.py --kraken-path {KRAKEN2_PATH} --kraken-db {KRAKEN2_DB} --fastq-in-one  r1.fastq.gz --fastq-in-two r1.fastq.gz --fastq-out-prefix {params_out_prefix}  --fastq-out-one {output_read_one} --fastq-out-two {output_read_two} --kraken-report {output_kraken_report} --kraken-output {output_kraken_output}  --threads {threads} 2> {log}".format(SCRIPTS_DIR=SCRIPTS_DIR,KRAKEN2_PATH=KRAKEN2_PATH, KRAKEN2_DB=KRAKEN2_DB, input_read_one=input.read_one, input_read_two=input.read_two, params_out_prefix=params.out_prefix,output_read_one=output.read_one, output_read_two=output.read_two,output_kraken_report=output.kraken_report, output_kraken_output=output.kraken_output,threads=threads,log=log)
                subprocess.check_call(run_cmd,shell=True)
        except subprocess.CalledProcessError as exc:
            pass


rule braken_assignment_merged:
    group: "align"
    input:
        kraken_report = "outputs/kraken2/merged/{sample}.kraken_report.txt",
    output:
        abundances = "outputs/abundances/merged/{sample}.bracken"
    log:
        "logs/abundances/{sample}.log"
    run:
        try:
            run_cmd="{BRACKEN_PATH} -d {KRAKEN2_DB} -i {input_kraken_report} -o {output_abundances} 2> {log}".format(BRACKEN_PATH=BRACKEN_PATH, KRAKEN2_DB=KRAKEN2_DB, input_kraken_report=input.kraken_report, output_abundances=output.abundances, log=log)
            subprocess.check_call(run_cmd,shell=True)
        except subprocess.CalledProcessError as exc:
            shell("touch {output.abundances}")
            pass


#### TODO: make an agg bracken that is filtered and unfiltered ####
### For that I need the final output file which will contain all the samples that pass QC
rule ag_braken:
    input:
        expand("outputs/abundances/merged/{sample}.bracken",sample=mapped_uid.keys())
    output:
        "outputs/final/abundances.bracken"
    run:
        with open(output[0],"w") as out_f:
            first_file = True 
            for in_f in input:
                sample = os.path.basename(in_f).split(".bracken")[0]
                with open(in_f) as read_in:
                    for line in read_in:
                        if first_file:
                            string_add = "sample\t"
                            first_file = False
                        else: 
                            string_add = sample  + "\t"
                        out_f.write(string_add + line)
 

rule kraken_deplete_host:
    group: "align"
    input:
        read_one = "outputs/kraken2/{sample}.cseqs__1.fastq",
        read_two = "outputs/kraken2/{sample}.cseqs__2.fastq"
    output:
        read_one= "outputs/depleted_fastqs/human/{sample}.R1.fastq.gz",
        read_two= "outputs/depleted_fastqs/human/{sample}.R2.fastq.gz"
    shell:
        "{SCRIPTS_DIR}/kraken_filter_fastqs.py --read-one-in {input.read_one} --read-two-in {input.read_two} --read-one-out {output.read_one} --read-two-out {output.read_two} --deplete" 

rule kraken_filter_reads:
    group: "align"
    input:
        read_one = "outputs/kraken2/{sample}.cseqs__1.fastq",
        read_two = "outputs/kraken2/{sample}.cseqs__2.fastq"
    output:
        read_one = "outputs/filtered_fastqs/{intervals}/{sample}.R1.fastq.gz",
        read_two = "outputs/filtered_fastqs/{intervals}/{sample}.R2.fastq.gz"
    resources:
        mem_mb=16000,
        walltime=43200
    params:
        taxon_to_extract=lambda wildcards: taxon_dict["{intervals}".format(intervals=wildcards.intervals)]
    shell:
        "{SCRIPTS_DIR}/kraken_filter_fastqs.py --read-one-in {input.read_one} --read-two-in {input.read_two} --read-one-out {output.read_one} --read-two-out {output.read_two} --taxon-to-extract {params.taxon_to_extract}" 


## Merge libraries that are the sample library type + sample id then rerun the pipeline for these samples. 
def get_bams(wildcards):
    samples = mapped_uid[wildcards.sample]
    return(expand("outputs/md/{intervals}/{sample}.sorted.md.bam", intervals=intervals, sample= samples))

def get_bai(wildcards):
    samples = mapped_uid[wildcards.sample]
    return(expand("outputs/md/{intervals}/{sample}.sorted.md.bam.bai", intervals=intervals, sample= samples))
 
def get_sample_name_for_merge(wildcards):
	sample_id = sample_names_hash[wildcards.sample][0]
	# TODO REMOVE
	return(sample_id)
rule merged_bams:
    shadow: "minimal"
    group: "align"
    input:
        bams = get_bams, 
        bai = get_bai
    output:
        bam="outputs/md/merged/{intervals}/{sample}.bam",
        bai="outputs/md/merged/{intervals}/{sample}.bam.bai"
    params:
        rg=get_sample_name_for_merge    
    run:
        if len(input.bams) != 1:
            # Skip merge
            shell("samtools merge tmp.bam {input.bams} && samtools sort tmp.bam > sorted.bam")
            shell("java -jar {PICARD_PATH} AddOrReplaceReadGroups I=sorted.bam O={output.bam} RGID=1 RGSM={params.rg} RGLB=4 RGPL=ILLUMINA RGPU=unit1")
        else:
            shell("java -jar {PICARD_PATH} AddOrReplaceReadGroups I={input.bams} O={output.bam} RGID=1 RGSM={params.rg} RGLB=4 RGPL=ILLUMINA RGPU=unit1")
        shell("samtools index {output}")


def get_rg(wildcards):
    strain_name =get_sample_name(wildcards)
    rg="@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:lib1\\tPU:unit1".format(sample=strain_name)
    return(rg)

rule bwa_map2: 
    input: 
        read_one="outputs/filtered_fastqs/{intervals}/{sample}.R1.fastq.gz",
        read_two="outputs/filtered_fastqs/{intervals}/{sample}.R2.fastq.gz"
    output:
        "outputs/kraken_bams/{intervals}/{sample}.kraken.bam"
    params:
        rg=get_rg
        ### Could get the reference seq here in a function base on the the taxid extraction
    resources:
        mem_mb=16000,
        walltime=43200
    threads: 1
    log:
        "logs/align/{intervals}/{sample}.log"
    shell:
        "{SCRIPTS_DIR}/align.sh '{params.rg}' {SARS_REF} {input.read_one} {input.read_two} {output} {log} {threads}" 

rule picard_sort:
    group: "align"
    input:
        ("outputs/kraken_bams/{intervals}/{sample}.kraken.bam")
    output:
        temp("outputs/sorted/{intervals}/{sample}.sorted.bam")
    resources:
        mem_mb=16000,
        walltime=93600
    params:
        order="SORT_ORDER=coordinate",
        tmpdir="-Djava.io.tmpdir=tmpdir",
    log:
        "logs/picard_sort/{intervals}/{sample}.log"
    shell:
        "{SCRIPTS_DIR}/picard_sort.sh {params.tmpdir} {PICARD_PATH} {input} {output} {params.order} {log} {resources.mem_mb}"


rule picard_mark_duplicates:
    input:
        ("outputs/sorted/{intervals}/{sample}.sorted.bam")
    output:
        bam_out=("outputs/md/{intervals}/{sample}.sorted.md.bam"),
        metrics_out="outputs/{intervals}/md/{sample}.metrics.txt"
    group: "align"
    log:
        "logs/picard_mark_duplicates/{intervals}/{sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir",
    shell:
        "{SCRIPTS_DIR}/picard_md.sh {params.tmpdir} {PICARD_PATH} {input} {output.bam_out} {output.metrics_out} {log}"

rule picard_index:
    input:
        "outputs/md/{intervals}/{sample}.sorted.md.bam"
    output:
        "outputs/md/{intervals}/{sample}.sorted.md.bam.bai"
    group: "align"
    log:
        "logs/picard_index/{intervals}/{sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir"
    shell:
        "{SCRIPTS_DIR}/picard_index.sh {params.tmpdir} {PICARD_PATH} {input} {output} {log}"

rule bwa_map:
    input:
        read_one="outputs/fastqs/{sample}.R1.fastq.gz",
        read_two="outputs/fastqs/{sample}.R2.fastq.gz"
    output:
        temp("outputs/mapping_stats/align/{sample}.bam")
    params:
        rg=get_rg
    group: "align2"
    resources:
        mem_mb=16000,
        walltime=43200
    threads: 4
    log:
        "logs/align2/{sample}.log"
    shell:
        #echo {input.read_one} {input.read_two} 2> {log}"
        "{SCRIPTS_DIR}/align.sh '{params.rg}' {REFERENCE} {input.read_one} {input.read_two} {output} {log} {threads} "
rule picard_sort2:
    input:
        ("outputs/mapping_stats/align/{sample}.bam")
    output:
        temp("outputs/mapping_stats/sorted/{sample}.sorted.bam")
    group: "align2"
    resources:
        mem_mb=16000,
        walltime=93600
    params:
        order="SORT_ORDER=coordinate",
        tmpdir="-Djava.io.tmpdir=tmpdir",
    log:
        "logs/picard_sort2/{sample}.log"
    shell:
        "{SCRIPTS_DIR}/picard_sort.sh {params.tmpdir} {PICARD_PATH} {input} {output} {params.order} {log} {resources.mem_mb}"


rule picard_mark_duplicates2:
    input:
        ("outputs/mapping_stats/sorted/{sample}.sorted.bam")
        #("outputs/sorted_mapping_stats/{sample}.sorted.bam")
    output:
        bam_out=("outputs/mapping_stats/md/{sample}.sorted.md.bam"),
        metrics_out="outputs/mapping_stats/md/{sample}.metrics.txt"
    group: "align2"
    resources:
        mem_mb=16000,
        walltime=93600
    log:
        "logs/picard_mark_duplicates2/{sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir",
    shell:
        "{SCRIPTS_DIR}/picard_md.sh {params.tmpdir} {PICARD_PATH} {input} {output.bam_out} {output.metrics_out} {log}"

rule picard_index2:
    input:
        ("outputs/mapping_stats/md/{sample}.sorted.md.bam")
    output:
        "outputs/mapping_stats/md/{sample}.sorted.md.bam.bai"
    group: "align2"
    resources:
        mem_mb=16000,
        walltime=93600
    log:
        "logs/picard_index/{sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir",
    shell:
        "{SCRIPTS_DIR}/picard_index.sh {params.tmpdir} {PICARD_PATH} {input} {output} {log}"
#rule trim_reads:
#    shadow: "minimal"
#    output:
#        read_one="outputs/trimmed_fastqs/{sample}.R1.fastq.gz",
#        read_two="outputs/trimmed_fastqs/{sample}.R2.fastq.gz"
#    log:
#        "outputs/trim_reads/{sample}.log#"
#    shell:
#        "{SCRIPTS_DIR}/trimmomatic.sh {input.read_one} {input.read_two} {output.read_one} {output.read_two} {TRIMMOMATIC_PATH} {TRIMMOMATIC_ADAPTERS} {log}"
    
#rule bwa_map:
#    input:
#        read_one = "outputs/filtered_fastqs/{sample}.R1.fastq.gz",
#        read_two = "outputs/filtered_fastqs/{sample}.R2.fastq.gz"
#        #read_one=get_fq1,
#        #read_two=get_fq2
#    #    read_one="outputs/trimmed_fastqs/{sample}.R1.fastq.gz",
#    #    read_two="outputs/trimmed_fastqs/{sample}.R2.fastq.gz"
#    group: "align"
#    output:
#        temp("outputs/mapped/{sample}.unsorted.bam")
#    resources:
#        mem_mb=16000,
#        walltime=43200
#    threads: 8
#    params:
#        rg="@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:lib1\\tPU:unit1",
#    log:
#        "logs/align/{sample}.log"
#    shell:
#        #echo {input.read_one} {input.read_two} 2> {log}"
#        "{SCRIPTS_DIR}/align.sh '{params.rg}' {SARS_REF} {input.read_one} {input.read_two} {output} {log}  "
#
#
#### Extract 
#rule extract_reads_mapping_to_sars2:
#    shadow: "minimal"
#    group: "variant_calling"
#    input:
#        bam="outputs/md/{sample}.sorted.md.bam",
#        bam_index="outputs/md/{sample}.sorted.md.bam.bai"
#    output:
#        bam="outputs/bam_subset/{sample}.{intervals}.bam",
#        bai="outputs/bam_subset/{sample}.{intervals}.bam.bai"
#    log:
#        "logs/picard_index/{sample}.{intervals}.log"
#    params:
#        tmpdir="-Djava.io.tmpdir=tmpdir",
#        cluster="-l h_rt=24:00:00 -l h_data=16G -cwd",
#        rg="@RG\\tID:{sample}\\tSM:{sample}\\tPL:illumina\\tLB:lib1\\tPU:unit1"
#    shell:
#        "{SCRIPTS_DIR}/subset_and_index_bam.sh  {input.bam} {output.bam} {intervals} {MAPQ_FILT} {log} {PICARD_PATH} {params.tmpdir} {SARS_REF} {params.rg}"
