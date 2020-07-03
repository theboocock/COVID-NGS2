rule run_pangolin:
    shadow: "shallow"
    conda: workflow.basedir + "/envs/pangolin.yaml"
    input:
        "outputs/consensus/merged/phylo/sars2/all_5.fasta"
    output:
        "outputs/pangolin/lineages.csv"
    threads:
        20
    shell:
        "{SCRIPTS_DIR}/pangolin.sh {input} {output} {threads}"

rule align_genomes_and_generate_tree:
    input:
        "outputs/consensus/all.fasta"
    output:
        "outputs/trees/all.tree"
    shell:
        "{SCRIPTS_DIR}/generate_consensus_tree.py -i {input} -o {output} "

