def get_pango_in_all(wildcards):

    return(f"outputs/consensus/merged/phylo/sars2/all_{MIN_DEPTH}_none.fasta")

rule run_pangolin:
    shadow: "shallow"
    conda: workflow.basedir + "/envs/pangolin.yaml"
    input:
        get_pango_in_all
    output:
        "outputs/pangolin/lineages.csv"
    threads:
        20
    shell:
        "{SCRIPTS_DIR}/pangolin.sh {input} {output} {threads}"
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
    shell:
        "{SCRIPTS_DIR}/ncov/assign_lineages.sh {input.input_fasta} {output.nextstrain_clades} {output.gisaid_clades} {GISAID_CLADES} {NEXTSTRAIN_CLADES} {SARS2_REF_GENBANK}"

##

##

