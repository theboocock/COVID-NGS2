rule haplotype_caller:
    """
        Run the GATK haplotype caller to get variants called against the SARS2-COV genome
    """
    group: "variant_calling"
    input:
        bam="outputs/sorted/{sample}.sorted.bam",
        bam_index="outputs/sorted/{sample}.sorted.bam.bai"
    output:
        "outputs/gvcfs/{sample}.g.vcf"
    threads: 8
    log:
        "logs/haplotype_caller/{sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir",
        cluster="-l time=24:00:00 -l h_data=1.5G -pe shared 12 -cwd"
    shell:
        "{SCRIPTS_DIR}/haplotype_caller.sh {GATK_PATH} {REFERENCE} {input.bam} {output} {threads}"

def _gatk_multi_arg(flag, files):
    cmd_string = flag + flag.join(files)
    return cmd_string

rule genomics_db_import:
    group: "vcf_output"
    input:
        gvcfs = expand("outputs/gvcfs/{sample}.g.vcf", sample=samples)
    output:
        "outputs/genomesdb/{intervals}_complete.txt",
    log:
        "logs/genomicsdb/{intervals}.log"
    run:
        shell("""{SCRIPTS_DIR}/make_genomicsdb.sh {GATK_PATH} {REFERENCE} {output} {intervals} {input.gvcfs}""")

rule output_vcf:
    group: "vcf_output"
    input:
        expand("outputs/genomesdb/{intervals}_complete.txt",intervals=intervals)
    log:
        "logs/vcfs/all_{intervals}.log"
    output:
        "outputs/vcfs/all_{intervals}.vcf.gz"
    run:
        shell("""{SCRIPTS_DIR}/extract_variants_from_genomes_db.sh {GATK_PATH} {input} {output} {log} {REFERENCE}""")

rule filter_vcf:
    group: "vcf_output"
    input:
        expand("outputs/vcfs/all_{intervals}.vcf.gz",intervals=intervals)
    output:
        "outputs/vcfs_filt/all_{intervals}.filt.vcf.gz"        
    log:
        "outputs/vcfs_filt/all_{intervals}.filt.log"
    shell:
        "{SCRIPTS_DIR}/filter_vcf.sh  {input} {output} {log} {MIN_DEPTH} {MIN_QUAL}"
