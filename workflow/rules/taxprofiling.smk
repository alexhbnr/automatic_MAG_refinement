localrules: decompress_bins

def return_bins_for_taxprofiling(wildcards):
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
        metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
        bins.extend(metawrap['sample_binID'].tolist())
    return [f"{config['tmpdir']}/taxprofiling_fas/{b}.fa" for b in bins]

rule taxprofiling:
    input:
        f"{config['tmpdir']}/gtdbtk/gtdbtk.bac120.summary.tsv",
        f"{config['tmpdir']}/phylophlan/phylophlan_closestGenomes.tsv"
    output:
        "{resultdir}/tax_profiling.tsv"
    message: "Summarise the GTDBTK and PhyloPhlAn Metagenomics results"
    params:
        gtdbtk_dir = f"{config['tmpdir']}/gtdbtk"
    run:
        # PhyloPhlAn
        phylophlan = pd.read_csv(input[1], sep="\t", skiprows=1, header=None, usecols=[0,1],
                                 names=['sample_binID', 'res'])
        tax_results = phylophlan['res'].str.split(":", expand=True)
        tax_results.columns = ['SGBtype', 'SGBrank', 'SGBlineage', 'SGBavgdist']
        tax_results['SGBavgdist'] = tax_results['SGBavgdist'].astype(float)
        tax_results['sample_binID'] = phylophlan['sample_binID']

        # GTDBTK
        gtdbtk_bact = pd.read_csv(f"{params.gtdbtk_dir}/gtdbtk.bac120.summary.tsv", sep="\t")
        if os.path.isfile(f"{params.gtdbtk_dir}/gtdbtk.ar53.summary.tsv"):
            gtdbtk_ar = pd.read_csv(f"{params.gtdbtk_dir}/gtdbtk.ar53.summary.tsv", sep="\t")
            gtdbtk = pd.concat([gtdbtk_bact, gtdbtk_ar])
        else:
            gtdbtk = gtdbtk_bact

        tax_results.merge(gtdbtk.rename({'user_genome': 'sample_binID',
                                         'classification': 'GTDBlineage'},
                                         axis=1),
                          how="left", on="sample_binID") \
        .iloc[:, [4, 5, 2, 0, 1, 3] + list(range(6, 24))] \
        .sort_values(['sample_binID']) \
        .to_csv(output[0], sep="\t", index=False, float_format="%.3f")


#### Decompress bins ###########################################################

rule decompress_bins:
    input:
        lambda wildcards: f"{config['resultdir']}/{wildcards.bin}/{wildcards.bin}.fasta.gz"
    output:
        temp("{tmpdir}/taxprofiling_fas/{bin}.fa")
    message: "De-compress {wildcards.bin}"
    shell:
        "gunzip -c {input} > {output}"

################################################################################

#### GTDBTK ####################################################################

rule gtdbtk_classify:
    input:
        db = f"{config['resourcesdir']}/gtdbtk/gtdbtk_r207_v2/metadata/metadata.txt",
        fas = lambda wildcards: return_bins_for_taxprofiling(wildcards)
    output:
        "{tmpdir}/gtdbtk/gtdbtk.bac120.summary.tsv"
    message: "Run the GTDBTK's classify workflow"
    container: "docker://quay.io/biocontainers/gtdbtk:2.3.2--pyhdfd78af_0"
    resources:
        mem = 80,
        cores = 32
    params:
        fadir = f"{config['tmpdir']}/taxprofiling_fas",
        outdir = f"{config['tmpdir']}/gtdbtk",
        dbdir = f"{config['resourcesdir']}/gtdbtk/gtdbtk_r207_v2"
    threads: 32
    shell:
        """
        GTDBTK_DATA_PATH={params.dbdir} \
            gtdbtk classify_wf --cpu {threads} \
            --mash_db {params.dbdir}/mash_db \
            --extension fa --genome_dir {params.fadir} --out_dir {params.outdir}
        if [[ ! -s {output}  ]]; then
            echo "WARNING: no MAGs from the kingdom bacterium found!"
            touch {output}
        fi
        """

################################################################################

#### PhyloPhlAn ################################################################

rule phylophlan_metagenomics_closestgenomes:
    input:
        lambda wildcards: return_bins_for_taxprofiling(wildcards)
    output:
        "{tmpdir}/phylophlan/phylophlan_closestGenomes.tsv"
    message: "Run PhyloPhlAn metagenomics analysis on bins returning set of 10 closest genomes per bin"
    container: "https://depot.galaxyproject.org/singularity/phylophlan:3.0.3--pyhdfd78af_0"
    resources:
        mem = 100,
        cores = 24
    params:
        dbdir = lambda wildcards: f"{config['resourcesdir']}/phylophlan_databases",
        database = config['phylophlan_db_version'],
        bindir = f"{config['tmpdir']}/taxprofiling_fas",
        prefix = f"{config['tmpdir']}/phylophlan/phylophlan_closestGenomes"
    threads: 24
    shell:
        """
        phylophlan_metagenomic \
            -i {params.bindir} \
            --input_extension .fa \
            -o {params.prefix} \
            --how_many 3 \
            --nproc {threads} \
            --verbose \
            --database {params.database} \
            --database_folder {params.dbdir}
        """


################################################################################
