#### Auxilliary functions ######################################################
def return_bins_of_sample(wildcards):
    metawrap = pd.read_csv(sampletsv.at[wildcards.sample, 'metawrapreport'],
                           sep="\t")
    metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
    metawrap['sample_binID'] = [f"{wildcards.sample}_{i:03}" for i in metawrap['binid']]
    return [f"{config['tmpdir']}/fastas/{b.split('_')[0]}/{b}.fa" for b in metawrap['sample_binID']]
################################################################################

rule evaluation_samplelevel:
    input:
        expand("{resultdir}/{sample}/{sample}_GUNC_checkM.merged.tsv", resultdir=[config['resultdir']], sample=SAMPLES)
    output:
        touch(f"{config['tmpdir']}/evaluation_samplelevel.done")

rule decompress_filtered_contigs:
    input:
        lambda wildcards: f"{config['resultdir']}/{wildcards.bin.split('_')[0]}/bins/{wildcards.bin}.fasta.gz"
    output:
        temp("{tmpdir}/fastas/{sample}/{bin}.fa")
    message: "De-compress automatically filtered contigs: {wildcards.bin}"
    shell:
        "gunzip -c {input} > {output}"

rule gunc_run:
    input:
        fastas = lambda wildcards: return_bins_of_sample(wildcards),
        db = lambda wildcards: f"{config['resourcesdir']}/GUNC/db"
    output:
        "{tmpdir}/gunc/{sample}/GUNC.progenomes_2.1.maxCSS_level.tsv"
    message: "Run GUNC on bins: {wildcards.sample}"
    resources:
        mem = 32,
        cores = 8
    params:
        db_file = lambda wildcards: f"{config['resourcesdir']}/GUNC/db/gunc_db_progenomes2.1.dmnd",
        fadir = "{tmpdir}/fastas/{sample}",
        outdir = "{tmpdir}/gunc/{sample}",
        tmpdir = "{tmpdir}/gunc/{sample}/tmp"
    threads: 8
    wrapper:
        "file:workflow/wrappers/gunc_run"

rule checkM_lineage_wf:
    input:
        lambda wildcards: return_bins_of_sample(wildcards)
    output:
        "{tmpdir}/checkM/{sample}/storage/marker_gene_stats.tsv"
    message: "Run checkM using lineage-specific workflow on sample {wildcards.sample}"
    resources:
        mem = 80,
        cores = 8
    params:
        fadir = "{tmpdir}/fastas/{sample}",
        outputfd = "{tmpdir}/checkM/{sample}"
    log: "{tmpdir}/checkM/{sample}.checkM.log"
    threads: 8
    wrapper:
        "file:workflow/wrappers/checkm_lineage_wf"

rule checkM_qa:
    input:
        "{tmpdir}/checkM/{sample}/storage/marker_gene_stats.tsv"
    output:
        "{tmpdir}/checkM/{sample}.checkM.txt"
    message: "Generate extended checkM report for sample {wildcards.sample}"
    resources:
        mem = 20,
        cores = 1
    params:
        outputfd = "{tmpdir}/checkM/{sample}"
    threads: 1
    wrapper:
        "file:workflow/wrappers/checkm_qa"

rule gunc_merge:
    input:
        gunc = lambda wildcards: f"{config['tmpdir']}/gunc/{wildcards.sample}/GUNC.progenomes_2.1.maxCSS_level.tsv",
        checkm = lambda wildcards: f"{config['tmpdir']}/checkM/{wildcards.sample}.checkM.txt"
    output:
        "{resultdir}/{sample}/{sample}_GUNC_checkM.merged.tsv"
    message: "Merge GUNC and CheckM output: {wildcards.sample}"
    resources:
        mem = 4,
        cores = 1
    params:
        dir = "{resultdir}/{sample}"
    threads: 1
    wrapper:
        "file:workflow/wrappers/gunc_merge"
