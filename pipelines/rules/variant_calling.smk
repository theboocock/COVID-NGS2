rule haplotype_caller:
    """
        Run the GATK haplotype caller to get variants called against the SARS2-COV genome
    """
    group: "variant_calling"
    input:
        bam="outputs/md/{intervals}/{sample}.sorted.md.bam",
        bam_index="outputs/md/{intervals}/{sample}.sorted.md.bam.bai"
    group: "haplotype"
    output:
        "outputs/gvcfs/{intervals}/{sample}.g.vcf"
    threads: 1
    log:
        "logs/haplotype_caller/{intervals}/{sample}.log"
    shell:
        "{SCRIPTS_DIR}/haplotype_caller.sh {GATK_PATH} {SARS_REF} {input.bam} {output} {threads}"

def _gatk_multi_arg(flag, files):
    cmd_string = flag + flag.join(files)
    return cmd_string


## joint var calling. 
#rule genomics_db_import:
#    group: "vcf_output"
#    input:
#        gvcfs = expand("outputs/gvcfs/{intervals}/{sample}.g.vcf", sample=samples, intervals=intervals)
#    output:
#        "outputs/genomesdb/{intervals}_complete.txt",
#    log:
#        "logs/genomicsdb/{intervals}.log"
#    run:
#        shell("""{SCRIPTS_DIR}/make_genomicsdb.sh {GATK_PATH} {REFERENCE} {output} {intervals} {input.gvcfs}""")
#
#rule output_vcf:
#    group: "vcf_output"
#    input:
#        expand("outputs/genomesdb/{intervals}_complete.txt",intervals=intervals)
#    log:
#        "logs/vcfs/all_{intervals}.log"
#    output:
#        "outputs/vcfs/{intervals}/all.vcf.gz"
#    run:
#        shell("""{SCRIPTS_DIR}/extract_variants_from_genomes_db.sh {GATK_PATH} {input} {output} {log} {REFERENCE}""")
#
rule filter_vcf:
    """
        Filter VCF 
    """
    group: "vcf_output"
    input:
        #expand("outputs/bcftools_vcfs/{intervals}/{sample}.vcf",intervals=intervals,sample=samples)
        expand("outputs/vcfs/{intervals}/all.vcf.gz",intervals=intervals)

    output:
        out_snps_only="outputs/vcfs_filt/{intervals}/all.filt.snps.vcf.gz", 
        out_all_vars="outputs/vcfs_filt/{intervals}/all.filt.all_vars.vcf.gz"        
    log:
        "outputs/vcfs_filt/{intervals}/all.filt.log"
    shell:
        "{SCRIPTS_DIR}/filter_vcf.sh {input} {output.out_snps_only} {output.out_all_vars} {log} {MIN_QUAL} {MIN_DEPTH}"

rule subset_vcf:
#    group: "sample_vcf:
    group: "vcf_output"
    input:
        "outputs/vcfs_filt/{intervals}/all.filt.all_vars.vcf.gz"
    output:
        "outputs/vcfs_filt/{intervals}/{sample}.vcf"
    log:
        "outputs/vcf_filt/{intervals}/{sample}.vcf"
    params:
        rg_name = get_sample_name
    shell:
        """bcftools view -s {params.rg_name} {input} | bcftools view -i 'GT="1/1"' > {output} 2> {log}"""
#    vcf="outputs/vcfs_filt/{intervals}/{sample}.vcf",

rule joint_call_vcfs:
    group: "vcf_output"
    input:
        bams=expand("outputs/mapping_stats/bam_subset/merged/{sample}.bam",sample=mapped_uid.keys()),
        bam_index=expand("outputs/mapping_stats/bam_subset/merged/{sample}.bam.bai",sample=mapped_uid.keys())
    output:
        "outputs/quasi_species/merged/all_from_merged.vcf"
    log:
        "log/quasi_speciies/merged/all_from_merged.log"
    run: 
        shell("bcftools mpileup -a FMT/AD -a FMT/DP -f {SARS_REF} {input.bams} | bcftools call -Ov -mv | bgzip -c > test.vcf.gz && tabix -p vcf test.vcf.gz")
        shell("bcftools view -e 'QUAL < {MIN_QUAL}' test.vcf.gz | bgzip -c > test1.vcf.gz")
        shell("{SCRIPTS_DIR}/vcf_filter/allelic_balance.py --vcf test1.vcf.gz -o {output}")




rule create_quasi_species_vcf:
    group: "vcf_output"
    input:
        "outputs/final/merged/all.vcf.gz"
    output:
        "outputs/final/merged/quasi_species/all.vcf"
    shell:
        "{SCRIPTS_DIR}/vcf_filter/allelic_balance.py --vcf {input} -o {output}"

rule create_merged_quasi_species_vcf:
    group: "vcf_output"
    input:
        "outputs/bcftools_vcfs/merged/filt/sars2/all.vcf.gz"
    output:
        "outputs/quasi_species/merged/all.vcf"
    shell:
        "{SCRIPTS_DIR}/vcf_filter/allelic_balance.py --vcf {input} -o {output}"


rule convert_annovar_to_table:
    group: "annovar"
    input:
        multi_anno="outputs/final/annovar/all.NC_045512v2_multianno.txt",
        all_vcf="outputs/final/merged/all.vcf.gz"
    output:
        "outputs/final/merged/annovar/all_filtered_annovar_table.txt"
    run:
        sample_list_vcf = """ bcftools query -l {0} """.format(input.all_vcf)
        sample_list = subprocess.check_output(sample_list_vcf, shell=True)
        sample_list = (sample_list.decode().split("\n"))
        sample_list = sample_list[:(len(sample_list)-1)]
        multi_tab = pd.read_csv(input.multi_anno, sep="\t")
        columns_first = multi_tab.columns[:10]
        columns_to_drop = ["Otherinfo"+ str(i) for i in range(2,11)]
        multi_tab = multi_tab.drop(columns_to_drop,axis=1)
        columns_list = list(columns_first) + ["allele_frequency","information","gt"] + list(sample_list)
        multi_tab.columns = columns_list 
        multi_tab.to_csv(output[0],sep="\t",index=False)

rule run_annovar:
    group: "annovar"
    input:
        unfiltered_vcf="outputs/bcftools_vcfs/merged/all.vcf.gz",
        filtered_vcf="outputs/bcftools_vcfs/merged/filt/sars2/all.vcf.gz",
        final_vcf="outputs/final/merged/all.vcf.gz"
    output:
        filtered_output="outputs/annovar/all_filtered.NC_045512v2_multianno.txt",
        final_output="outputs/final/annovar/all.NC_045512v2_multianno.txt"
    run:
        directory_name = os.path.dirname(output.filtered_output)
        shell("{SOFTWARE_PATH}/annovar/table_annovar.pl --protocol avGene  --operation g  --buildver NC_045512v2 -vcfinput {input.unfiltered_vcf} --outfile {directory_name}/all {SOFTWARE_PATH}/annovar/sar2_db")
        shell("{SOFTWARE_PATH}/annovar/table_annovar.pl --protocol avGene  --operation g  --buildver NC_045512v2 -vcfinput {input.filtered_vcf} --outfile {directory_name}/all_filtered {SOFTWARE_PATH}/annovar/sar2_db ") 
        directory_name = os.path.dirname(output.final_output)
        shell("{SOFTWARE_PATH}/annovar/table_annovar.pl --protocol avGene  --operation g  --buildver NC_045512v2 -vcfinput {input.final_vcf} --outfile {directory_name}/all {SOFTWARE_PATH}/annovar/sar2_db ") 
#rule bcftools:
#    group: "vcf_output"
#    input:
#        bam="outputs/md/{intervals}/{sample}.sorted.md.bam",
#        bam_index="outputs/md/{intervals}/{sample}.sorted.md.bam.bai"
#    output:
#        "outputs/bcftools_vcfs/{intervals}/{sample}.vcf.gz"
#    log:
#        "log/bcftools_vcfs/{intervals}/{sample}.log"
#    params:
#    shell:
#        "bcftools mpileup -a FMT/DP -a FMT/AD -f {SARS_REF} {input.bam} | bcftools call -Ov -m | bgzip -c > {output} && tabix -p vcf {output}"
#
#rule bcftools_merge:
#    group: "vcf_output"
#    input:
#        bam="outputs/md/merged/{intervals}/{sample}.bam",
#        bam_index="outputs/md/merged/{intervals}/{sample}.bam.bai"
#    output:
#        "outputs/bcftools_vcfs/merged/{intervals}/{sample}.vcf.gz"
#    log:
#        "log/bcftools_vcfs/merged/{intervals}/{sample}.log"
#    shell:
#        "bcftools mpileup -a FMT/AD -a FMT/DP -f {SARS_REF} {input.bam} | bcftools call -Ov -m | bgzip -c > {output} && tabix -p vcf {output}"
rule bcftools:
    group: "vcf_output"
    input:
        bam="outputs/mapping_stats/bam_subset/{sample}.bam",
        bam_index="outputs/mapping_stats/bam_subset/{sample}.bam.bai"
    output:
        "outputs/bcftools_vcfs/{intervals}/{sample}.vcf.gz"
    log:
        "log/bcftools_vcfs/{intervals}/{sample}.log"
    params:
    shell:
        "bcftools mpileup -a FMT/DP -a FMT/AD -f {SARS_REF} {input.bam} | bcftools call -Ov -m | bgzip -c > {output} && tabix -p vcf {output}"

rule bcftools_merge:
    group: "vcf_output"
    input:
        bam="outputs/mapping_stats/bam_subset/merged/{sample}.bam",
        bam_index="outputs/mapping_stats/bam_subset/merged/{sample}.bam.bai"
    output:
        "outputs/bcftools_vcfs/merged/{intervals}/{sample}.vcf.gz"
    log:
        "log/bcftools_vcfs/merged/{intervals}/{sample}.log"
    shell:
        "bcftools mpileup -a FMT/AD -a FMT/DP -f {SARS_REF} {input.bam} | bcftools call -Ov -m | bgzip -c > {output} && tabix -p vcf {output}"
#rule bcftools_gvcf:
#    group: "vcf_output"
#    input:
#        bam="outputs/md/merged/{intervals}/{sample}.bam",
#        bam_index="outputs/md/merged/{intervals}/{sample}.bam.bai"
#    output:
#        # "outputs/bcftools_vcfs/merged/{intervals}/{sample}.vcf.gz"
#        "outputs/bcftools_vcfs/merged/{intervals}/{sample}.g.vcf"
#    log:
#        "log/bcftools_vcfs/merged/{intervals}/{sample}.log"
#    shell:
#        "bcftools mpileup -a FMT/AD -f {SARS_REF} {input.bam} | bcftools call -g -Ov -mv | bgzip -c > {output}"



rule combine_vcfs:
    group: "vcf_output"
    input:
        vcfs = expand("outputs/bcftools_vcfs/merged/{intervals}/{sample}.vcf.gz",intervals=intervals, sample=mapped_uid.keys())
    output:
        "outputs/bcftools_vcfs/merged/all.vcf.gz"
    shell:
        "bcftools merge {input.vcfs} | bgzip -c  > {output} && tabix -p vcf {output}"

rule combine_filtered_vcf:
    input:
        vcfs = expand("outputs/bcftools_vcfs/merged/filt/{intervals}/{sample}.snps_and_indels.vcf.gz",intervals=intervals, sample=mapped_uid.keys())
    output:
        "outputs/bcftools_vcfs/merged/filt/sars2/all.vcf.gz"
    shell:
        #"bcftools merge {input.vcfs} | bcftools view -v snps,indels | bedtools intersect -v -a stdin -b {SITES_TO_MASK} -header | bgzip -c > {output} &&  tabix -p vcf {output}"
        "bcftools merge {input.vcfs} | bcftools view -v snps,indels | bgzip -c > {output} &&  tabix -p vcf {output}"


rule annovar_vcf:
    input:
        "outputs/bcftools_vcfs/merged/filt/{intervals}/all.vcf.gz"
    output:
        "outputs/annovar/output.txt"
    shell:
        "annovar "
rule vcf_filter:
    group: "vcf_output"
    input:
        "outputs/bcftools_vcfs/{intervals}/{sample}.vcf.gz"
    output:
        out_snps_only="outputs/bcftools_vcfs/filt/{intervals}/{sample}.snps.vcf.gz",
        out_all_vars="outputs/bcftools_vcfs/filt/{intervals}/{sample}.snps_and_indels.vcf.gz"
    log:
        "log/bcftools_vcfs/merged/{intervals}/{sample}.log"
    shell:
        "{SCRIPTS_DIR}/filter_vcf.sh {input} {output.out_snps_only} {output.out_all_vars} {log} {MIN_QUAL} {MIN_DEPTH}"
       
rule vcf_merge:
    group: "vcf_output"
    input:
        "outputs/bcftools_vcfs/merged/{intervals}/{sample}.vcf.gz"
    output:
        out_snps_only="outputs/bcftools_vcfs/merged/filt/{intervals}/{sample}.snps.vcf.gz",
        out_all_vars="outputs/bcftools_vcfs/merged/filt/{intervals}/{sample}.snps_and_indels.vcf.gz"
    log:
        "log/bcftools_vcfs/merged/{intervals}/{sample}.log"
    shell:
        "{SCRIPTS_DIR}/filter_vcf.sh {input} {output.out_snps_only} {output.out_all_vars} {log} {MIN_QUAL} {MIN_DEPTH}"

rule aggregate_consensus_from_merged:
    """
        Aggregate the consensus genomes at a range of depths.
    """
    group: "ag"
    input:
        # Hardcoded depth filter
        input_fasta=expand("outputs/consensus/merged/{intervals}/{depth}/{sample}.fasta",sample=mapped_uid.keys(),intervals=intervals,depth=MIN_DEPTH),
        input_phylo=expand("outputs/consensus/merged/phylo/{intervals}/{depth}/{sample}.fasta",sample=mapped_uid.keys(),intervals=intervals,depth=MIN_DEPTH)
    output:
        fasta_one="outputs/consensus/merged/{intervals}/all_{depth}.fasta",
        fasta_two="outputs/consensus/merged/phylo/{intervals}/all_{depth}.fasta",
        all_fasta="outputs/consensus/merged/{intervals}/all_{depth}_none.fasta",
        all_fasta_phylo="outputs/consensus/merged/phylo/{intervals}/all_{depth}_none.fasta"
    run:
        coverage_file = ""
        with open(output.fasta_one,"w") as out_f:
            samples = []
            for f in input.input_fasta:
                coverage_file = f.split(".fasta")[0] + ".cov"
                with open(coverage_file) as in_cov:
                    coverage = float(in_cov.readline().strip())
                    if coverage >= min_cov:
                        with open(f) as in_f:
                            i = 0
                            for line in in_f:
                                out_f.write(line)
                                if i == 0:
                                    samples.append(line.strip().split(">")[1])
                                i = i + 1
        with open(output.all_fasta,"w") as out_f:
            for f in input.input_phylo:
                with open(f) as in_f:
                    i = 0
                    for line in in_f:
                        out_f.write(line)
                        i = i + 1
        with open(output.all_fasta_phylo,"w") as out_f:
            for f in input.input_phylo:
                with open(f) as in_f:
                    i = 0
                    for line in in_f:
                        out_f.write(line)
                        i = i + 1
        with open(output.fasta_two,"w") as out_f:
            for f in input.input_phylo:
                    with open(f) as in_f:
                        i = 0
                        for line in in_f:
                            if i == 0:
                                sample = (line.strip().split(">")[1])
                                if sample not in samples:
                                    break
                            out_f.write(line)
                            i = i + 1
        #TODO: ensure I only keep the columns I want. 
### TODO: Add new coverage pipeline
rule get_consensus_from_merged:
    shadow: "minimal"
    group: "vcf_output"
    input:
        interval = "outputs/bcftools_vcfs/merged/filt/{intervals}/{sample}.snps_and_indels.vcf.gz",
        quasi_in= "outputs/quasi_species/merged/all.vcf",
        bam="outputs/mapping_stats/bam_subset/merged/{sample}.bam",
        coverage="outputs/coverage/merged/sars2/{sample}.cov"
    output:
        fasta="outputs/consensus/merged/{intervals}/{depth}/{sample}.fasta",
        coverage="outputs/consensus/merged/{intervals}/{depth}/{sample}.cov",
        fasta_phylo="outputs/consensus/merged/phylo/{intervals}/{depth}/{sample}.fasta",
    params:
        sample=get_sample_name_for_merge
    run:
        shell("{SCRIPTS_DIR}/generate_consensus.py  --coverage-in {input.coverage} -v {input.interval} -d {wildcards.depth} -p 1.0 -o {output.fasta_phylo} -c /dev/null -r {SARS_REF}        -i {input.bam} -s {params.sample} -m {SITES_TO_MASK} --quasi-vcf {input.quasi_in}")
        shell("{SCRIPTS_DIR}/generate_consensus.py  --coverage-in {input.coverage} -v {input.interval} -d {wildcards.depth} -p 1.0 -o {output.fasta} -c {output.coverage} -r {SARS_REF}        -i {input.bam} -s {params.sample} --quasi-vcf {input.quasi_in}")


rule aggregate_consensus:
    """
        Aggregate the consensus genomes at a range of depths.
    """
    input:
        # Hardcoded depth filter
        expand("outputs/consensus/{intervals}/{{depth}}/{sample}.fasta",sample=samples,intervals=intervals)
    output:
        "outputs/consensus/{intervals}/all_{depth}.fasta",
    run:
        with open(output[0],"w") as out_f:
            samples = []
            for f in input:
                coverage = f.split(".fasta")[0] + ".cov"
                with open(coverage) as in_cov:
                    coverage = float(in_cov.readline().strip())
                    if coverage >= min_cov:
                        with open(f) as in_f:
                            i = 0
                            for line in in_f:
                                out_f.write(line)
                                if i == 0:
                                    samples.append(line.strip().split(">")[1])
                                i = i + 1
#checkpoint check_good_samples:



rule get_consensus:
    shadow: "minimal"
    group: "vcf_output"
    input:
        interval=  "outputs/bcftools_vcfs/filt/{intervals}/{sample}.snps_and_indels.vcf.gz",
        bam="outputs/mapping_stats/bam_subset/{sample}.bam",
        coverage="outputs/coverage/{intervals}/{sample}.cov"
    output:
        fasta="outputs/consensus/{intervals}/{depth}/{sample}.fasta",
        coverage="outputs/consensus/{intervals}/{depth}/{sample}.cov",
        fasta_phylo="outputs/consensus/phylo/{intervals}/{depth}/{sample}.fasta"
    params:
        sample=get_sample_name
    run:
        shell("{SCRIPTS_DIR}/generate_consensus.py  --coverage-in {input.coverage} -v {input.interval} -d {wildcards.depth} -p 1.0 -o {output.fasta_phylo} -c /dev/null -r {SARS_REF}        -i {input.bam} -s {params.sample} -m {SITES_TO_MASK} --quasi-vcf {input.interval}")
        shell("{SCRIPTS_DIR}/generate_consensus.py  --coverage-in {input.coverage} -v {input.interval} -d {wildcards.depth} -p 1.0 -o {output.fasta} -c {output.coverage} -r {SARS_REF}        -i {input.bam} -s {params.sample} --quasi-vcf {input.interval}")
        

rule map_consensus_to_genome:
    shadow: "minimal"
    group: "vcf_output"
    input:
        fasta="outputs/consensus/{intervals}/{depth}/{sample}.fasta"
    output:
        bam="outputs/consensus/bams/{intervals}/{sample}/{depth}.bam",
        bai="outputs/consensus/bams/{intervals}/{sample}/{depth}.bam.bai"
    params:
        rg=get_rg
        ### Could get the reference seq here in a function base on the the taxid extraction
    resources:
        mem_mb=16000,
        walltime=43200
    threads: 1
    shell:
        "{SCRIPTS_DIR}/quick_align.sh '{params.rg}' {SARS_REF} {input.fasta} {output.bam}"
