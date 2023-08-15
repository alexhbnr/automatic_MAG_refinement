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
    run:
        with open(output[0], "wt") as outfile:
            for fn in params.fas:
                input_empty = True
                if params.type == "contigs":
                    sample = os.path.basename(os.path.dirname(os.path.dirname(fn)))
                elif params.type == "genes":
                    sample = os.path.basename(fn).replace("-contigs_superkingdom_phylum.fa.gz", "")
                # Test whether input FastA file is empty
                if params.type == "genes":
                    with open(fn, "rb") as f:
                        fstart = f.read(3)
                        if fstart.startswith(b"\x1f\x8b\x08"):  # gzipped
                            if os.stat(fn).st_size > 50:
                                input_empty = False
                        else:
                            if os.stat(fn).st_size > 0:
                                input_empty = False
                # Write to file if not empty
                if not input_empty or params.type == "contigs":
                    for name, seq in pyfastx.Fasta(fn, build_index=False):
                        outfile.write(f">{sample}:{name}\n{seq}\n")

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
        cores = 24
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
    script:
        "../scripts/mmseqs2_annotatetsv.py"

rule mmseqs2_splittsv:
    input:
        "{tmpdir}/mmseqs2/concat/all_contigs.mmseqs2_gtdb.annot.tsv"
    output:
        touch("{tmpdir}/mmseqs2/mmseqs2_splittsv.done")
    message: "Split MMSeqs2 result table into bins"
    params:
        sampletsv = config['sampletsv'],
        dir = "{tmpdir}/mmseqs2"
    run:
        # Import MMSeqs2 results
        mmseqs2 = pd.read_csv(input[0], sep="\t")
        mmseqs2['sample'] = mmseqs2['contig'].str.split(":").str[0]
        mmseqs2['contig'] = mmseqs2['contig'].str.split(":").str[1]

        # Import sample TSV
        sampletsv = pd.read_csv(params.sampletsv, sep="\t", index_col=[0])
        sampletsv.index = sampletsv.index.astype(str)

        for sample in mmseqs2['sample'].unique():
            # Load overview list of bins
            sample_bins = pd.read_csv(sampletsv.at[sample, 'metawrapreport'], sep="\t")
            sample_bins['binid'] = sample_bins['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
            sample_bins['sample_binID'] = [f"{sample}_{i:03}" for i in sample_bins['binid']]
            sample_bins = sample_bins.set_index(['sample_binID'])
            # Load map of contigs per bin
            contigs = pd.read_csv(sampletsv.at[sample, 'metawrapreport'].replace(".stats", ".contigs"),
                                sep="\t", header=None, names=['contig', 'bin'], index_col=['bin'])
            for b in sample_bins.index:
                c = contigs.loc[sample_bins.at[b, 'bin']]['contig'].tolist()
                mmseqs2.loc[(mmseqs2['sample'] == sample) & (mmseqs2['contig'].isin(c))] \
                    .drop(['sample'], axis=1) \
                    .to_csv(f"{params.dir}/{b}.mmseqs2_gtdb.annot.tsv", sep="\t", index=False)
