rule star_align:
    shadow: "minimal"
    input:
        read_one=get_fq1,
        read_two=get_fq2
    resources:
        mem_mb=16000,
        walltime=43200
    output:
        "outputs/star/{sample}.bam"
    threads: 8
    params:
        rg="ID:{sample} SM:{sample} PL:illumina LB:lib1 PU:unit1",
    log:
        "logs/star/{sample}.log"
    shell:
        "{SCRIPTS_DIR}/star_align.sh {threads} {input.read_one} {input.read_two} {STAR_REFERENCE} '{params.rg}' {output}  {log}"

rule picard_sort:
    input:
        ("outputs/star/{sample}.bam")
    output:
        ("outputs/sorted/{sample}.sorted.bam")
    resources:
        mem_mb=16000,
        walltime=93600
    params:
        order="SORT_ORDER=coordinate",
        tmpdir="-Djava.io.tmpdir=tmpdir",
    log:
        "logs/picard_sort/{sample}.log"
    shell:
        "{SCRIPTS_DIR}/picard_sort.sh {params.tmpdir} {PICARD_PATH} {input} {output} {params.order} {log} {resources.mem_mb}"


rule picard_index:
    input:
        "outputs/sorted/{sample}.sorted.bam"
    output:
        "outputs/sorted/{sample}.sorted.bam.bai"
    log:
        "logs/picard_index/{sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir",
    shell:
        "{SCRIPTS_DIR}/picard_index.sh {params.tmpdir} {PICARD_PATH} {input} {output} {log}"

