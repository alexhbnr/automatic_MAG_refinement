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
        expand("{resultdir}/sample_stats/{sample}_GUNC_checkM_checkM2.merged.tsv", resultdir=[config['resultdir']], sample=SAMPLES),
    output:
        touch(f"{config['tmpdir']}/evaluation_samplelevel.done")

rule decompress_filtered_contigs:
    input:
        lambda wildcards: f"{config['resultdir']}/{wildcards.bin}/{wildcards.bin}.fasta.gz"
    output:
        temp("{tmpdir}/fastas/{sample}/{bin}.fa")
    message: "De-compress automatically filtered contigs: {wildcards.bin}"
    shell:
        "gunzip -c {input} > {output}"

#### GUNC ######################################################################

rule gunc_run:
    input:
        fastas = lambda wildcards: return_bins_of_sample(wildcards),
        db = lambda wildcards: f"{config['resourcesdir']}/GUNC/db"
    output:
        "{tmpdir}/gunc/{sample}/GUNC.progenomes_2.1.maxCSS_level.tsv"
    message: "Run GUNC on bins: {wildcards.sample}"
    container: "https://depot.galaxyproject.org/singularity/gunc:1.0.5--pyhdfd78af_0"
    resources:
        mem = 32,
        cores = 8
    params:
        db_file = lambda wildcards: f"{config['resourcesdir']}/GUNC/db/gunc_db_progenomes2.1.dmnd",
        fadir = "{tmpdir}/fastas/{sample}",
        outdir = "{tmpdir}/gunc/{sample}",
        tmpdir = "{tmpdir}/gunc/{sample}/tmp"
    threads: 8
    shell:
        """
        mkdir -p {params.tmpdir}
        gunc run \
            -r {params.db_file} \
            --input_dir {params.fadir} \
            -t {threads} \
            -o {params.outdir} \
            --temp_dir {params.tmpdir} \
            --detailed_output
        rmdir {params.tmpdir}
        """

################################################################################

#### checkM1 ###################################################################

rule checkM_lineage_wf:
    input:
        lambda wildcards: return_bins_of_sample(wildcards)
    output:
        "{tmpdir}/checkM/{sample}/storage/marker_gene_stats.tsv"
    message: "Run checkM using lineage-specific workflow on sample {wildcards.sample}"
    container: "https://depot.galaxyproject.org/singularity/checkm-genome:1.2.1--pyhdfd78af_0"
    resources:
        mem = 80,
        cores = 8
    params:
        fadir = "{tmpdir}/fastas/{sample}",
        outputfd = "{tmpdir}/checkM/{sample}"
    log: "{tmpdir}/checkM/{sample}.checkM.log"
    threads: 8
    shell:
        """
        rm -r {params.outputfd}    
        checkm lineage_wf \
            -x fa \
            -t {threads} \
            {params.fadir} \
            {params.outputfd} > {log}

        """

rule checkM_qa:
    input:
        "{tmpdir}/checkM/{sample}/storage/marker_gene_stats.tsv"
    output:
        "{tmpdir}/checkM/{sample}.checkM.txt"
    message: "Generate extended checkM report for sample {wildcards.sample}"
    container: "https://depot.galaxyproject.org/singularity/checkm-genome:1.2.1--pyhdfd78af_0"
    resources:
        mem = 20,
        cores = 1
    params:
        outputfd = "{tmpdir}/checkM/{sample}"
    threads: 1
    shell:
        """
        checkm qa -o 2 --tab_table -f {output[0]} \
            {params.outputfd}/lineage.ms \
            {params.outputfd}
        """

################################################################################

#### CheckM2 ###################################################################

rule checkm2:
    input:
        lambda wildcards: return_bins_of_sample(wildcards)
    output:
        "{tmpdir}/checkm2/{sample}/quality_report.tsv"
    message: "Evaluate the quality of the MAGs sample {wildcards.sample}"
    container: "https://depot.galaxyproject.org/singularity/checkm2:1.0.1--pyh7cba7a3_0"
    resources:
        mem = 36,
        cores = 8
    params:
        fadir = "{tmpdir}/fastas/{sample}",
        outdir = "{tmpdir}/checkm2/{sample}",
        db_path = f"{config['resourcesdir']}/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
    threads: 8
    shell:
        """
        checkm2 predict --threads {threads} \
            --input {params.fadir} \
            --output-directory {params.outdir} \
            --force \
            --database_path {params.db_path} \
            -x fa
        """

################################################################################

#### Summary ###################################################################

rule gunc_merge:
    input:
        gunc = lambda wildcards: f"{config['tmpdir']}/gunc/{wildcards.sample}/GUNC.progenomes_2.1.maxCSS_level.tsv",
        checkm = lambda wildcards: f"{config['tmpdir']}/checkM/{wildcards.sample}.checkM.txt"
    output:
        temp("{resultdir}/sample_stats/{sample}_GUNC_checkM.merged.tsv")
    message: "Merge GUNC and CheckM output: {wildcards.sample}"
    container: "https://depot.galaxyproject.org/singularity/gunc:1.0.5--pyhdfd78af_0"
    resources:
        mem = 4,
        cores = 1
    params:
        dir = "{resultdir}/sample_stats/{sample}"
    threads: 1
    shell:
        """
        mkdir -p {params.dir}
        gunc merge_checkm \
        -g {input.gunc} \
        -c {input.checkm} \
        -o {params.dir}
        mv {params.dir}/GUNC_checkM.merged.tsv {output}
        rmdir {params.dir}
        """

rule checkm2_merge:
    input:
        gunc = "{resultdir}/sample_stats/{sample}_GUNC_checkM.merged.tsv",
        checkm2 = lambda wildcards: f"{config['tmpdir']}/checkm2/{wildcards.sample}/quality_report.tsv"
    output:
        "{resultdir}/sample_stats/{sample}_GUNC_checkM_checkM2.merged.tsv"
    message: "Merge checkM2 with GUNC and checkM: {wildcards.sample}"
    resources:
        mem = 4,
        cores = 1
    run:
        pd.read_csv(input.gunc, sep="\t") \
            .merge(pd.read_csv(input.checkm2, sep="\t") \
                   .rename({'Name': 'genome',
                            'Completeness': 'checkM2.completeness',
                            'Contamination': 'checkM2.contamination',
                            'Completeness_Model_Used': 'checkM2.completeness_model_used',
                            'Coding_Density': 'checkM2.coding_density',
                            'Total_Coding_Sequences': 'checkM2.total_coding_sequences',
                            'Additional_Notes': 'checkM2.notes'}, axis=1) \
                   .drop(['Translation_Table_Used', 'Contig_N50', 'Average_Gene_Length',
                          'Genome_Size', 'GC_Content'], axis=1),
                   how="left", on="genome") \
            .to_csv(output[0], sep="\t", index=False, float_format="%.3f")

################################################################################
