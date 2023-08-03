#### Auxilliary functions ######################################################

def return_bin_fas(wildcards):
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        bins.extend([f"{sampletsv.at[sample, 'metawrapreport'].replace('.stats', '')}/{b}.fa"
                     for b in metawrap['bin']])
    return bins

################################################################################

rule mmseqs2_gtdb_contigs:
    input:
        "{tmpdir}/mmseqs2/mmseqs2_splittsv.done"
    output:
        touch("{tmpdir}/mmseqs2_gtdb_contigs.done")
    run:
        for sample in SAMPLES:
            metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                                sep="\t")
            metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
            metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
            for b in metawrap['sample_binID']:
                if not os.path.isfile(f"{config['tmpdir']}/mmseqs2/{b}.mmseqs2_gtdb.annot.tsv"):
                    sys.exit(1)

rule concat_bins:
    output:
        temp("{tmpdir}/mmseqs2/concat/all_contigs.fasta")
    message: "Concatenate all contigs that were binned into a single FastA"
    resources:
        mem = 4,
        cores = 1
    params:
        fas = lambda wildcards: return_bin_fas(wildcards),
        type = "contigs"
    wrapper:
        "file:workflow/wrappers/concat_bins"

rule createdb_bins:
    input:
        "{tmpdir}/mmseqs2/concat/all_contigs.fasta"
    output:
        "{tmpdir}/mmseqs2/concat/all_contigs.contigs"
    message: "Create database of contigs"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 16,
        cores = 1
    shell:
        "mmseqs createdb {input} {output}"

rule screen:
    input:
        db = f"{config['resourcesdir']}/mmseqs2/gtdb/mmseqs2_gtdb_r207_db_mapping",
        contigs = "{tmpdir}/mmseqs2/concat/all_contigs.contigs"
    output:
        "{tmpdir}/mmseqs2/concat/all_contigs.mmseqs2_gtdb.index"
    message: "Assign taxonomy via the GTDB for contigs"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 750,
        cores = 36
    params:
        tmpdir = "{tmpdir}/mmseqs2_tmpdir",
        prefix = "{tmpdir}/mmseqs2/concat/all_contigs.mmseqs2_gtdb",
        gtdb_db = f"{config['resourcesdir']}/mmseqs2/gtdb/mmseqs2_gtdb_r207_db"
    threads: 24
    shell:
        """
        if [[ $(stat -c%s {input.contigs}) -ge 50 ]]; then
            mmseqs taxonomy \
                {input.contigs} \
                {params.gtdb_db} \
                {params.prefix} \
                {params.tmpdir} \
                -a \
                --tax-lineage 1 \
                --lca-ranks kingdom,phylum,class,order,family,genus,species \
                --majority 0.5 \
                --vote-mode 1 \
                --orf-filter 1 \
                --remove-tmp-files 1 \
                --threads {threads}
        else
            touch {output}
        fi
        """

rule create_tsv:
    input:
        contigs = "{tmpdir}/mmseqs2/concat/all_contigs.contigs",
        assignments = "{tmpdir}/mmseqs2/concat/all_contigs.mmseqs2_gtdb.index"
    output:
        pipe("{tmpdir}/mmseqs2/concat/all_contigs.mmseqs2_gtdb.tsv")
    message: "Convert MMSeqs2 GTDB results to TSV"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 8,
        cores = 1
    params:
        prefix = "{tmpdir}/mmseqs2/concat/all_contigs.mmseqs2_gtdb",
    shell:
        """
        if [[ $(stat -c%s {input.contigs}) -ge 50 ]]; then
            mmseqs createtsv {input.contigs} {params.prefix} {output}
        else
            touch {output}
        fi
        """

rule mmseqs2_annotatetsv:
    input:
        "{tmpdir}/mmseqs2/concat/all_contigs.mmseqs2_gtdb.tsv"
    output:
        "{tmpdir}/mmseqs2/concat/all_contigs.mmseqs2_gtdb.annot.tsv"
    message: "Add header to MMSeqs2 table"
    wrapper:
        "file:workflow/wrappers/mmseqs2_annotatetsv"

rule mmseqs2_splittsv:
    input:
        "{tmpdir}/mmseqs2/concat/all_contigs.mmseqs2_gtdb.annot.tsv"
    output:
        touch("{tmpdir}/mmseqs2/mmseqs2_splittsv.done")
    message: "Split MMSeqs2 result table into bins"
    params:
        sampletsv = config['sampletsv'],
        dir = "{tmpdir}/mmseqs2"
    wrapper:
        "file:workflow/wrappers/mmseqs2_splittsv"
