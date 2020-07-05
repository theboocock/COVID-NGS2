rule picard_insert_size_metrics:
    """
        Run PicardTools CollectInsertSizeMetrics 

        https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics
    """
    message:
        "Calculating insert size metrics"
    group: "ag"
    input:
        bam="outputs/mapping_stats/md/{sample}.sorted.md.bam",
        bam_index="outputs/mapping_stats/md/{sample}.sorted.md.bam.bai"
        #bam="outputs/bam_subset/{sample}.sorted.md.bam",
        #bam_index="outputs/bam_subset/{sample}.sorted.md.bam.bai"
    group: "qc_sample"
    output:
        insert="outputs/insert_sizes/{sample}.insert_sizes",
        histogram=report("outputs/insert_sizes/{sample}.histogram.pdf",caption=workflow.basedir + "/report/insert_sizes.rst", category="insert_sizes")
    log:
        "logs/picard_insert_size/{sample}.log"
    params:
        cluster="-l h_rt=24:00:00 -l h_data=16G -cwd",
        tmpdir="-Djava.io.tmpdir=tmpdir"
    shell:
        "touch {output.insert} && touch {output.histogram} && {SCRIPTS_DIR}/insert_size_metrics.sh {params.tmpdir} {PICARD_PATH} {input.bam} {output.insert} {output.histogram}"

rule generate_coverage_plot:
    """
        #TODO: rename rule
        Use bedtools to get depth per bp.
    
    """
    shadow: "minimal"
    resources:
        mem_mb=16000,
        walltime=43200
    shadow: "minimal"
    group: "ag"
    input:
        bam="outputs/mapping_stats/bam_subset/{sample}.bam"
    output:
        "outputs/coverage/{intervals}/{sample}.cov"
    params:
        sample_type = get_sample_type
    run:
        if "neb" in params.sample_type:
            shell("samtools view -h -F 1024 -bS {input.bam} > test.bam")
            shell("samtools index test.bam")
            shell("bedtools genomecov -ibam test.bam -d > {output}")
        else:
            shell("bedtools genomecov -ibam {input.bam} -d > {output}")

rule make_coverage_plot:
    """
        Create coverage plots 
    """
    #report("outputs/coverage/coverage_plots/{intervals}/{sample}.png", caption="report/coverage_plots.rst", category="coverage_plots")
    group: "ag"
    conda: workflow.basedir + "/envs/r-covid.yaml"
    input:    
        "outputs/coverage/{intervals}/{sample}.cov"
    output:
        #report("outputs/coverage/coverage_plots/{intervals}/{sample}.png", caption=workflow.basedir+ "/report/coverage_plots.rst", category="coverage_plots")
        "outputs/coverage/coverage_plots/{intervals}/{sample}.png"
#        report("test.png", caption="report/coverage_plots.rst", category="coverage_plots")
    shell:
        "{SCRIPTS_DIR}/make_coverage_plots.R {input} {output} {wildcards.sample}_{wildcards.intervals}"

rule ag_qc_output:
    """
        Aggregate the per-sample QC metrics.
    """
    group: "ag"
    input:
        expand("outputs/qc_report/{intervals}/{sample}.qc",sample=samples, intervals=intervals)
    output:
        "outputs/qc_report/all.qc"
    run:
        with open(output[0],"w") as out_f:
            out_f.write("sample coverage_5x coverage_mean mapped_human mapped_rrna mapped_sars2 unmapped total_read_count median_insert_size mad_isert_size sars2_percent\n") 
            for f in input:
                sample = os.path.basename(f).split(".qc")[0]
                with open(f) as in_f:
                    for line in in_f:
                        line_split = line.split()
                        total_mapped =  (int(line_split[2]) + int(line_split[3]) + int(line_split[4])) 
                        if total_mapped == 0:
                            total = 0 
                        else:
                            total = 100 * int(line_split[4]) / int(float(line_split[6]))
                        line_split[6] = int(float(line_split[6]))
                        out_f.write(sample + " " + line.strip() + " " + str(total) + "\n")

import pandas as pd
rule merge_meta_data:
    """
        Merge the per sample QC metrics that were generated with the per sample metadata
    """
    group: "ag"
    input:
        in_qc = "outputs/qc_report/all.qc",
        pangolin_in = "outputs/pangolin/lineages.csv"
    output:
        "outputs/final/qc_report.tsv"
    run:
        in_qc = pd.read_csv(input[0], sep=" ")
        md_m= metadata.merge(in_qc,left_on="sample",right_on="sample")
        pango_in = pd.read_csv(input.pangolin_in,sep=",")
        md_m = md_m.merge(pango_in,left_on="strain",right_on="taxon",how="outer")
        md_m.to_csv(output[0],sep="\t", index=False)
        

def get_reference(wildcards):
    """
        Fuction to get the reference genome it should match the genome based on
        the interval. 
    """
    return(SARS_REF)

rule igv_report:
    """
        Create igv.js reports for each sample
    """
    group: "ag"
    input:
        fasta=get_reference,
        vcf="outputs/vcfs_filt/{intervals}/{sample}.vcf",
        # any number of additional optional tracks, see igv-reports manual
        tracks=[expand("outputs/consensus/bams/{{intervals}}/{{sample}}/{depth}.bam",depth=depth_range),"outputs/md/{intervals}/{sample}.sorted.md.bam"]
    output:
        ## Removed from the report because it was horrible. 
        ##report("outputs/igv-report/{intervals}/{sample}.html", caption=workflow.basedir+ "/report/igv-report.rst", category="igv-view")
        "outputs/igv-report/{intervals}/{sample}.html" 
    params:
        extra=""  # optional params, see igv-reports manual
    log:
        "logs/igv-report/{intervals}/{sample}/igv-report.log"
    wrapper:
        "0.55.1/bio/igv-reports"


rule calculate_mapping_statistics:
    """
        Extract human, sars2, rrna, and unmapped reads. 
    """
    group: "qc_sample"
    input:
        bam=("outputs/mapping_stats/md/{sample}.sorted.md.bam"),
        bai=("outputs/mapping_stats/md/{sample}.sorted.md.bam.bai"),
        bam_subset="outputs/mapping_stats/bam_subset/{sample}.bam"
    output:
        "outputs/mapping_stats/{sample}.mapstat"
    log:
        "logs/mapping_stats/{sample}.log"
    params:
        tmpdir="-Djava.io.tmpdir=tmpdir"
    shell:
        "{SCRIPTS_DIR}/calculate_mapping_statistics.py -b {input.bam} -c {CONFIGFILE}  -o {output} --bam-subset {input.bam_subset}"

rule qc_output:
    """
        combine QC metrics that were generated into a per sample QC metric file
    """
    group: "ag"
    input: 
        mapstat="outputs/mapping_stats/{sample}.mapstat",
        coverage=expand("outputs/consensus/{{intervals}}/{MIN_DEPTH}/{{sample}}.cov",MIN_DEPTH=MIN_DEPTH),
        insert_sizes="outputs/insert_sizes/{sample}.insert_sizes",
        coverage_full="outputs/coverage/{intervals}/{sample}.cov",
        read_count="outputs/read_counts/{sample}.txt"
    output:
        "outputs/qc_report/{intervals}/{sample}.qc"
    shell:
        "{SCRIPTS_DIR}/combine_metrics.py -m {input.mapstat} -c {input.coverage} -i {input.insert_sizes} -o {output} -a {input.coverage_full} --read-count {input.read_count} "


rule amplicon_qc_output:
    input:
        bam="outputs/mapping_stats/bam_subset/{sample}.bam",
        bai="outputs/mapping_stats/bam_subset/{sample}.bam.bai"
    output:
        temp("outputs/amplicon_qc/{sample}.check_run"),
    params:
        sample_type=get_sample_type, 
        oligos_used=get_oligos_used
    run:
        if params.sample_type == "amp":
            shell("{SCRIPTS_DIR}/amplicon/per_amplicon_depth.py -s amp -b {input.bam} -o {output} -q 1 -p {OLIGO_POOL_BED}")
        elif params.sample_type == "amp_complicated": 
            #### TODO: make this work for the new set of oligopools. 
            shell("{SCRIPTS_DIR}/amplicon/per_amplicon_depth.py -s amp -b {input.bam} -o {output} -q 1 -p {ALL_PRIMERS} -c --oligos-used {params.oligos_used}")
        else: 
            shell("echo NOTAMP > {output}")

rule combine_amplicon_qc_output:
    input:
        "outputs/amplicon_qc/amplicon_df.txt",
    output:
        "outputs/amplicon_qc/amplicon_merged_df.txt"
    run:
        pd_in = pd.read_csv(input[0],sep="\t")
        for key, items in mapped_uid.items():
            items_df = []
            for item in items:
                pd_tmp = pd_in[pd_in["sample"] == item]
                items_df.append(pd_tmp)
            out_df = items_df[0]
            out_df["sample"] = key
            if len(items_df) > 1:
                for df in items_df[1:]:
                    out_df["count"] = df["count"]
        out_df.to_csv(output[0], sep="\t",index=False)          

rule aggregate_amplicon_qc:
    input:
        expand("outputs/amplicon_qc/{sample}.check_run",sample=samples)
    output:
        amp_df="outputs/amplicon_qc/amplicon_df.txt",
        amp_complicated_df="outputs/amplicon_qc/amplicon_df_complicated.txt"
    run:
        with open(output.amp_df,"w") as out_f:
            with open(output.amp_complicated_df,"w") as out_f2:
                out_f.write("sample\tamplicon\tcount\tlibrary_type\n")
                out_f2.write("sample\tamplicon\tcount\toligo_set\toligos_used\n")
                for f in input:
                    sample_in = os.path.basename(f).split(".check_run")[0]
                    with open(f) as in_f:
                        line = in_f.readline()
                        if "NOTAMP" in line:
                            continue
                        if "oligos_used" in line:
                            for i, line in enumerate(in_f):
                                if i > 0:
                                    out_f2.write(sample_in + "\t" + line)
                        else:  
                            for i, line in enumerate(in_f):
                                if i > 0:
                                    out_f.write(sample_in + "\t" + line)
                        #out_f.write(sample_in +"\t" + line)
                # otherwise add the other files to the script
    

rule auspice_output:
    """
        Create the auspic output file at a minimum coverage. The defaults 
        are >5x for >90% of the sars2 genome 
    """
    input:
        fasta="outputs/final/merged/all.fasta",
        tsv="outputs/final/merged/qc_report.tsv"
    output:
        fasta="outputs/ncov/merged/all.fasta",
        tsv="outputs//ncov/merged/qc_report.tsv"
    shell:
        "{SCRIPTS_DIR}/combine_with_public_data.sh {input.fasta} {input.tsv} {GISAID_FASTA} {GISAID_METADATA} {output.fasta} {output.tsv}"

