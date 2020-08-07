def get_pango_in_all(wildcards):

    return(f"outputs/consensus/merged/phylo/sars2/all_{MIN_DEPTH}_none.fasta")

rule merge_in_pango:
    input:
        merged_csv="outputs/final/merged/qc_report.tsv.inco",
        merged_csv_unfiltered="outputs/qc_report/merged/qc_report.tsv.inco",
        pangolin_out="outputs/pangolin/lineages_merged.csv",
    output:
        merged_csv="outputs/final/merged/qc_report.tsv",
        merged_csv_unfiltered="outputs/qc_report/merged/qc_report.tsv"
    run:
        pango_in = pd.read_csv(input.pangolin_out, sep="\t")
        out_df = pd.read_csv(input.merged_csv,sep="\t")
        out_df = out_df.merge(pango_in, left_on="read_group_name", right_on="taxon", how="outer")
        out_df2= pd.read_csv(input.merged_csv_unfiltered,sep="\t")
        out_df2 = out_df2.merge(pango_in, left_on="read_group_name", right_on="taxon", how="outer")
        out_df.to_csv(output.merged_csv, sep="\t",index=False)
        out_df2.to_csv(output.merged_csv_unfiltered, sep="\t",index=False)

rule run_pangolin:
    shadow: "shallow"
    conda: workflow.basedir + "/envs/pangolin.yaml"
    input:
        get_pango_in_all
    output:
        "outputs/pangolin/lineages.csv"
    threads:
        1
    shell:
        "{SCRIPTS_DIR}/pangolin.sh {input} {output} {threads} TRUE"

def get_pango_in_double_ups(wildcards):
    return(f"outputs/consensus/merged/phased/sars2/all_{MIN_DEPTH}_none.fasta")

rule run_pangolin_phased:
    shadow: "shallow"
    conda: workflow.basedir + "/envs/pangolin.yaml"
    input:
        get_pango_in_double_ups
    output:
        "outputs/pangolin/doubleups.csv"
    threads:
        1
    shell:
        "{SCRIPTS_DIR}/pangolin.sh {input} {output} {threads} TRUE "

rule assign_lineages_phased:
    conda: workflow.basedir +"/envs/nextstrain.yaml"
    input:
        input_fasta=get_pango_in_double_ups
    threads:
        2
    output:
        nextstrain_clades="outputs/clades/nextstrain/doubleups.csv",
        gisaid_clades="outputs/clades/gisaid/doubleups.csv"
    shell:
        "{SCRIPTS_DIR}/ncov/assign_lineages.sh {input.input_fasta} {output.nextstrain_clades} {output.gisaid_clades} {GISAID_CLADES} {NEXTSTRAIN_CLADES} {SARS2_REF_GENBANK} {threads}"

rule create_imputted_vcf:
    input:
        vcf_input="outputs/final/merged/quasi_species/all.vcf",
        qc_input="outputs/qc_output/merged/qc_report.tsv"
    output:
        impute_variants="outputs/impute/merged/all.vcf.gz"
    run:
        shell("{SCRIPTS_DIR}/impute/impute_vcf.R {input.vcf_input} {output.impute_variants} {meta_in}") 

rule get_consensus_from_merged_impute:
    shadow: "minimal"
    group: "vcf_output"
    input:
        interval = "outputs/vcf/filter/{sample}_{down}.vcf.gz", 
        quasi_in="outputs/impute/merged/all.vcf.gz",
        bam="outputs/downsamples/{sample}_{down}.bam",
        coverage="outputs/coverage/{sample}_{down}.cov"
    output:
        fasta="outputs/consensus/imputted/{sample}_{down}.fasta"
    params:
        sample=get_sample_depth
    run:
        shell("{SCRIPTS_DIR}/generate_consensus.py  --coverage-in {input.coverage} -v {input.interval} -d {MIN_DEPTH} -p 1.0 -o {output.fasta} -c /dev/null -r {SARS_REF}        -i {input.bam} -s {params.sample} --quasi-vcf {input.quasi_in} --imputted ")
#rule run_pangolin:
#    shadow: "shallow"
#    conda: workflow.basedir + "/envs/pangolin.yaml"
#    input:
#        "outputs/consensus/merged/phylo/sars2/all_3.fasta"
#    output:
#        "outputs/pangolin/lineages.csv"
#    threads:
#        20
#    shell:
#        "{SCRIPTS_DIR}/pangolin.sh {input} {output} {threads}"

rule align_genomes_and_generate_tree:
    input:
        "outputs/consensus/all.fasta"
    output:
        "outputs/trees/all.tree"
    shell:
        "{SCRIPTS_DIR}/generate_consensus_tree.py -i {input} -o {output} "


rule assign_lineages:
    conda: workflow.basedir +"/envs/nextstrain.yaml"
    input:
        input_fasta=get_pango_in_all
    output:
        nextstrain_clades="outputs/clades/nextstrain/lineages.csv",
        gisaid_clades="outputs/clades/gisaid/lineages.csv"
    threads:
        2
    shell:
        "{SCRIPTS_DIR}/ncov/assign_lineages.sh {input.input_fasta} {output.nextstrain_clades} {output.gisaid_clades} {GISAID_CLADES} {NEXTSTRAIN_CLADES} {SARS2_REF_GENBANK} {threads}"

##

##

rule phased_final:
    input:
        "outputs/final/merged/quasi_species/all.vcf",
    output:
        beagle5="outputs/quasi_species/merged/phased.vcf.gz",
        beagle5_tbi="outputs/quasi_species/merged/phased.vcf.gz.tbi",
        beagle4="outputs/quasi_species/merged/phased_beagle4.vcf.gz"
    run:
        shell("java -jar {SOFTWARE_PATH}/beagle/beagle.18May20.d20.jar gt={input} out=test")
        shell("cp test.vcf.gz {output.beagle5}")
        shell("tabix -p vcf {output.beagle5}")
        shell("java -jar {SOFTWARE_PATH}/beagle/beagle.27Jan18.7e1.jar gt={input} out=test")
        shell("cp test.vcf.gz {output.beagle4}")
        shell("tabix -p vcf {output.beagle4}")

rule create_phased_output:
    shadow: "minimal"
    input:
        interval = "outputs/bcftools_vcfs/merged/filt/{intervals}/{sample}.snps_and_indels.vcf.gz",
        quasi_in="outputs/quasi_species/merged/phased.vcf.gz",
        quasi_tbx="outputs/quasi_species/merged/phased.vcf.gz.tbi",
        bam="outputs/mapping_stats/bam_subset/merged/{sample}.bam",
        coverage="outputs/coverage/merged/sars2/{sample}.cov"
    output:
        fasta="outputs/consensus/merged/phased/{intervals}/{depth}/{sample}.fasta",
        fasta_two="outputs/consensus/merged/phased/{intervals}/{depth}/{sample}_2.fasta"
    params:
        sample=get_sample_name_for_merge
    run:
        shell("touch {output.fasta}")
        shell("touch {output.fasta_two}")
        shell("{SCRIPTS_DIR}/generate_consensus.py  --coverage-in {input.coverage} -v {input.interval} -d {wildcards.depth} -p 1.0 -o {output.fasta} -c /dev/null -r {SARS_REF}        -i {input.bam} -s {params.sample} -m {SITES_TO_MASK} --quasi-vcf {input.quasi_in} --phased")

rule aggregate_phased_output: 
    group: "ag"
    input:
        input_fasta=expand("outputs/consensus/merged/phased/{intervals}/{depth}/{sample}.fasta",sample=mapped_uid.keys(),intervals=intervals,depth=MIN_DEPTH)
    output:
        fasta_one="outputs/consensus/merged/phased/{intervals}/all_{depth}_none.fasta"
    run:
        with open(output.fasta_one,"w") as out_f:
            samples=[]
            for f in input.input_fasta:
                sample_name = os.path.basename(f).split(".fasta")[0]
                with open(f) as in_f:
                    for i, line in enumerate(in_f):
                        if i == 0:
                            out_f.write(">" + sample_name + ":1"+ "\n")
                        else:
                            out_f.write(line)
                f2 = f.split(".fasta")[0] + "_2.fasta"
                with open(f2) as in_f:
                    for i, line in enumerate(in_f):
                        if i == 0:
                            out_f.write(">" + sample_name + ":2"+ "\n")
                        else:
                            out_f.write(line)


#rule all:
#    input:
#        expand("outputs/consensus/merged/phased/{intervals}/all_{depth}_none.fasta", depth=depth_range, intervals=intervals),
#        "outputs/quasi_species/merged/all_from_merged.vcf",
#        "outputs/quasi_species/merged/phased.vcf.gz",
#        expand("outputs/mapping_stats/bam_subset/merged/{uid_lib}.bam", uid_lib=mapped_uid.keys()),
#        "outputs/amplicon_qc/umi/amplicon_df.txt",
#        expand("outputs/mapping_stats/md/merged/{uid_lib}.bam", uid_lib=mapped_uid.keys()),
#        #expand("outputs/kraken2/{sample}.kraken_report.txt",sample=samples),
#        expand("outputs/consensus/merged/{intervals}/all_{depth}.fasta", depth=depth_range, intervals=intervals),
#        expand("outputs/abundances/merged/{sample}.bracken", sample=mapped_uid.keys()), 
#        expand("outputs/coverage/merged/coverage_plots/{intervals}/{sample}.png", intervals=intervals, sample=mapped_uid.keys()), 
#        expand("outputs/mapping_stats/bam_subset/{sample}.bam", sample=samples),
#        "outputs/ncov/merged/all.fasta",
#        "outputs/bcftools_vcfs/merged/all.vcf.gz",
#        "outputs/quasi_species/merged/all.vcf",
##        expand("outputs/raw_fastqs/{sample}.R1.fastq.gz",sample=samples),
#       # TODO check if it works#
#      #"test.png",
#        expand("outputs/coverage/coverage_plots/{intervals}/{sample}.png", intervals=intervals, sample=samples), 
#        #expand("outputs/kraken_bams/{intervals}/{sample}.kraken.bam",intervals=intervals,sample=samples),
#        #expand("outputs/kraken_bams/{intervals}/{sample}.kraken.bam",sample=samples,intervals=intervals),
#        #expand("outputs/vcfs_filt/{intervals}/all.filt.vcf.gz",intervals=intervals),
#        #expand("outputs/consensus/{depth}/{intervals}/{sample}.fasta",depth=depth_range,sample=samples,intervals=intervals),
#        #expand("outputs/coverage/{depth}/{intervals}/{sample}.cov", depth=depth_range, sample=samples, intervals=intervals),
#        expand("outputs/coverage/{intervals}/{sample}.cov", sample=samples,intervals=intervals),
#        #expand("outputs/insert_sizes/{intervals}/{sample}.insert_sizes",sample=samples,intervals=intervals),
#        expand("outputs/consensus/{intervals}/all_{depth}.fasta", depth=depth_range, intervals=intervals),
#        #expand("outputs/abundances/{sample}.bracken",sample=samples),
#        #"{SCRIPTS_DIR}/quick_align.sh '{params.rg}' {SARS_REF} {input.read_one} {input.read_two} {output} {log}" 
#        expand("outputs/consensus/bams/{intervals}/{sample}/{depth}.bam",intervals=intervals,sample=samples, depth=depth_range),
#        expand("outputs/depleted_fastqs/human/{sample}.R1.fastq.gz",sample=samples),
#       # bai="outputs/consensus/bams/{intervals}/{sample}/{depth}.bam.bai"
#        #expand("outputs/igv-report/{intervals}/{sample}.html", intervals=intervals, sample=samples),
#        # TODO: Reimplement the mapping from 
#        #
#        # TODO: Combinne all sample information we can gather from the sequencing into a single file. 
#        "outputs/amplicon_qc/amplicon_df.txt",
#        "outputs/amplicon_qc/amplicon_merged_df.txt",
#        "outputs/final/abundances.bracken",
#        "outputs/pangolin/lineages.csv",
#        #"outputs/final/qc_report.tsv",
#        "outputs/qc_report/qc_report.tsv",
#        directory(expand("outputs/final/merged/coverage_plots/{intervals}/",intervals=intervals)),
#        #""outputs/final/filtered/metadata.csv",
#        "outputs/final/merged/qc_report.tsv",
#        expand("outputs/bcftools_vcfs/merged/filt/{intervals}/all.vcf.gz",intervals=intervals),
#        "outputs/final/merged/annovar/all_filtered_annovar_table.txt",
#        "outputs/final/merged/quasi_species/all.vcf"
#  #     "outputs/consensus/merged/sars2/all_5.tsv",
#
import shutil
import numpy as np

rule merge_pango_outputs:
    input:
        pangolin_out="outputs/pangolin/lineages.csv",
        pango_double_ups="outputs/pangolin/doubleups.csv",
        nextstrain_doubleups="outputs/clades/nextstrain/doubleups.csv",
        gisaid_doubleups="outputs/clades/gisaid/doubleups.csv",
        nextstrain_clade="outputs/clades/nextstrain/lineages.csv",
        gisaid_clades="outputs/clades/gisaid/lineages.csv"
    output:
        pangolin_out_merged="outputs/pangolin/lineages_merged.csv"
    run:
        shell("touch {output.pangolin_out_merged}")
        pango_in = pd.read_csv(input.pangolin_out, sep=",")
        pango_double_ups =pd.read_csv(input.pango_double_ups, sep=",")
        pango_double_ups["original_id"] = pango_double_ups["taxon"].str.replace(":[1-2]","").str.replace("_","-")
        nextstrain_in  = pd.read_csv(input.nextstrain_clade,sep="\t") 
        nextstrain_doubleups = pd.read_csv(input.nextstrain_doubleups,sep="\t") 
        nextstrain_doubleups["original_id"] = nextstrain_doubleups["name"].str.replace(":[1-2]","").str.replace("_","-")
        gisaid_in= pd.read_csv(input.gisaid_clades ,sep="\t") 
        gisaid_doubleups = pd.read_csv(input.gisaid_doubleups,sep="\t") 
        gisaid_doubleups["original_id"] = gisaid_doubleups["name"].str.replace(":[1-2]","").str.replace("_","-")
#        pango_in = pango_in.merge(pango_double_ups,left_on="taxon",right_on="original_id",how="outer",suffixes=("_u","_p"))
        
        gisaid_tuple = []
        nextstrain_tuple = []
        pango_tuple = []

        for sample in pango_in["taxon"]:
            pango_tmp = pango_double_ups[pango_double_ups["original_id"].isin([sample])]
            lineages = ",".join(pango_tmp["lineage"])
            probabilities = ",".join([str(o) for o in pango_tmp["probability"]])
            nextstrain_tmp = nextstrain_doubleups[nextstrain_doubleups["original_id"].isin([sample])]
            lineage_nextstrain = ",".join(nextstrain_tmp["clade"])
            gisaid_tmp= gisaid_doubleups[gisaid_doubleups["original_id"].isin([sample])]
            lineage_gisaid= ",".join(gisaid_tmp["clade"])
            pango_tuple.append([sample, lineages, probabilities, lineage_nextstrain, lineage_gisaid]) 
   
        pango_out_df = pd.DataFrame(pango_tuple, columns=["sample","pangolin_phased","prob_phased","nextstrain_phased","gisaid_phased"])
        pango_in = pango_in.merge(pango_out_df,left_on="taxon",right_on="sample")
 #       pango_in.to_csv(output.pangolin_out_merged,sep=",",header=True)
        # Merge the phased genome lineages into the main dataset.
        #



#        pango_in = pango_in.merge(pango_double_ups,left_on="taxon",right_on=["original_id"],suffixes=("_u","_p"))
        pango_in.to_csv(output.pangolin_out_merged,sep="\t",header=True, index=False)
        
        # Next strain in ##
        #

   #     gisaid_in= pd.read_csv(input.gisaid_clades ,sep="\t") 
  #      gisaid_doubleups = pd.read_csv(input.gisaid_doubleups,sep="\t") 
