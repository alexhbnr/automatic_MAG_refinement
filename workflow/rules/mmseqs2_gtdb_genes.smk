#### Auxilliary functions ######################################################

def return_genes_fas(wildcards):
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
        metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
        bins.extend(metawrap['sample_binID'].tolist())
    return [f"{config['tmpdir']}/mmseqs2_genes/{b}-contigs_superkingdom_phylum.fa.gz" for b in bins]

def path_to_bin_fa(wildcards):
    binid = f"bin.{int(wildcards.bin.split('_')[1])}"
    metawrap_path = sampletsv.at[wildcards.bin.split('_')[0], 'metawrapreport'] \
        .replace(".stats", "")
    return f"{metawrap_path}/{binid}.fa"

def return_filter_fasta_files(wildcards):
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
        metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
        bins.extend(metawrap['sample_binID'].tolist())
    return [f"{config['resultdir']}/{b}/{b}.fasta.gz" for b in bins]

################################################################################

checkpoint mmseqs2_gtdb_genes:
    input:
        return_filter_fasta_files
    output:
        touch(f"{config['tmpdir']}/filter_contigs.done")

rule pyrodigal:
    output:
        gff = "{tmpdir}/pyrodigal/{bin}.gff",
        ffn = "{tmpdir}/pyrodigal/{bin}.ffn",
    container: "docker://quay.io/biocontainers/pyrodigal:2.3.0--py311h031d066_0"
    resources:
        mem = 8,
        cores = 4
    params:
        fa = lambda wildcards: path_to_bin_fa(wildcards)
    threads: 4
    shell:
        """
        pyrodigal \
            -i {params.fa} \
            -o {output.gff} \
            -d {output.ffn} \
            -p meta \
            -j {threads}
        """

rule identify_kingdom_phylum_contigs:
    input:
        mmseqs2 = "{tmpdir}/mmseqs2_gtdb_contigs.done",
        gff = "{tmpdir}/pyrodigal/{bin}.gff",
        ffn = "{tmpdir}/pyrodigal/{bin}.ffn"
    output:
        fa = "{tmpdir}/mmseqs2_genes/{bin}-contigs_superkingdom_phylum.fa.gz",
        mclineage = "{tmpdir}/mmseqs2_genes/{bin}-most_common_lineage.txt"
    message: "Identify contigs that could only be assigned to taxonomic rank superkingdom or phylum: {wildcards.bin}"
    resources:
        mem = 4,
        cores = 1
    params:
        fa = lambda wildcards: path_to_bin_fa(wildcards),
        prokkadir = "{tmpdir}/pyrodigal",
        tsv = "{tmpdir}/mmseqs2/{bin}.mmseqs2_gtdb.annot.tsv",
    script:
        "../scripts/identify_kingdom_phylum_contigs.py"

rule concat_bins_kp_contigs:
    input:
        fas = lambda wildcards: return_genes_fas(wildcards)
    output:
        temp("{tmpdir}/mmseqs2_genes/kp_contigs.fasta")
    message: "Concatenate all contigs that were assigned to ranks kingdom and phylum into a single FastA"
    resources:
        mem = 4,
        cores = 1
    params:
        fas = lambda wildcards: return_genes_fas(wildcards)
    run:
        with open(output[0], "wt") as outfile:
            for fn in params.fas:
                input_empty = True
                sample = os.path.basename(fn).replace("-contigs_superkingdom_phylum.fa.gz", "")
                # Test whether input FastA file is empty
                with open(fn, "rb") as f:
                    fstart = f.read(3)
                    if fstart.startswith(b"\x1f\x8b\x08"):  # gzipped
                        if os.stat(fn).st_size > 50:
                            input_empty = False
                    else:
                        if os.stat(fn).st_size > 0:
                            input_empty = False
                # Write to file if not empty
                if not input_empty:
                    for name, seq in pyfastx.Fasta(fn, build_index=False):
                        outfile.write(f">{sample}:{name}\n{seq}\n")

rule createdb_bins_kp_contigs:
    input:
        "{tmpdir}/mmseqs2_genes/kp_contigs.fasta"
    output:
        "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_genes"
    message: "Create database of genes on contigs that were assigned to ranks kingdom and phylum"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 16,
        cores = 1
    shell:
        "mmseqs createdb {input} {output}"

rule screen_kp_contigs:
    input:
        db = f"{config['resourcesdir']}/mmseqs2/gtdb/mmseqs2_gtdb_r207_db_mapping",
        contigs = "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_genes"
    output:
        "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes.index"
    message: "Assign taxonomy via the GTDB for genes on contigs that were assigned to ranks kingdom and phylum"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 750,
        cores = 24
    params:
        tmpdir = "{tmpdir}/mmseqs2_tmpdir",
        prefix = "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes",
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

rule create_tsv_kp_contigs:
    input:
        contigs = "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_genes",
        assignments = "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes.index"
    output:
        temp("{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes.tsv")
    message: "Convert MMSeqs2 GTDB results to TSV"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 16,
        cores = 1
    params:
        prefix = "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes",
    shell:
        """
        if [[ $(stat -c%s {input.contigs}) -ge 50 ]]; then
            mmseqs createtsv {input.contigs} {params.prefix} {output}
        else
            touch {output}
        fi
        """

rule annotate_tsv_kp_contigs:
    input:
        "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes.tsv"
    output:
        "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes.annot.tsv"
    message: "Add header to MMSeqs2 table"
    script:
        "../scripts/mmseqs2_annotatetsv.py"

rule mmseqs2_splittsv_genes:
    input:
        "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes.annot.tsv"
    output:
        touch("{tmpdir}/mmseqs2_genes/mmseqs2_splittsv.done")
    message: "Split MMSeqs2 result table into bins"
    params:
        sampletsv = config['sampletsv'],
        dir = "{tmpdir}/mmseqs2_genes"
    run:
        # Import MMSeqs2 results
        mmseqs2 = pd.read_csv(input[0], sep="\t")
        mmseqs2['sample_binID'] = mmseqs2['contig'].str.split(":").str[0]
        mmseqs2['sample'] = mmseqs2['sample_binID'].str.split("_").str[0]
        mmseqs2['contig'] = mmseqs2['contig'].str.split(":").str[1]

        # Generate list of all expected bins
        sampletsv = pd.read_csv(params.sampletsv, sep="\t", index_col=[0])
        expected_bins = []
        for sample in sampletsv.index:
            # Load overview list of bins
            sample_bins = pd.read_csv(sampletsv.at[sample, 'metawrapreport'], sep="\t")
            sample_bins['binid'] = sample_bins['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
            sample_bins['sample_binID'] = [f"{sample}_{i:03}" for i in sample_bins['binid']]
            expected_bins.extend(sample_bins['sample_binID'].tolist())

        # Iterate overall bins: if genes were screened, export results, others create
        # an empty file
        for binid in expected_bins:
            if binid in mmseqs2['sample_binID'].tolist():
                mmseqs2.loc[mmseqs2['sample_binID'] == binid] \
                    .drop(['sample_binID'], axis=1) \
                    .to_csv(f"{params.dir}/{binid}.mmseqs2_gtdb.annot.tsv", sep="\t", index=False)
            else:
                Path(f"{params.dir}/{binid}.mmseqs2_gtdb.annot.tsv").touch()

rule filter_contigs:
    input:
        genes_mmseqs2 = lambda wildcards: f"{config['tmpdir']}/mmseqs2_genes/mmseqs2_splittsv.done",
        depth = lambda wildcards: f"{config['tmpdir']}/depth/{wildcards.bin}.breadth_depth.txt"
    output:
        fa = "{resultdir}/{bin}/{bin}.fasta.gz",
        contigs = "{resultdir}/{bin}/{bin}.contigs.txt"
    message: "Automatically filter contigs that have a high chance to not belong to major lineage: {wildcards.bin}"
    params:
        fa = lambda wildcards: path_to_bin_fa(wildcards),
        contigs_mmseqs2 = lambda wildcards: f"{config['tmpdir']}/mmseqs2/{wildcards.bin}.mmseqs2_gtdb.annot.tsv",
        genes_mmseqs2 = lambda wildcards: f"{config['tmpdir']}/mmseqs2_genes/{wildcards.bin}.mmseqs2_gtdb.annot.tsv",
        prokkadir = lambda wildcards: f"{config['tmpdir']}/pyrodigal"
    run:
        # Read coverage
        coverage = pd.read_csv(input.depth, sep="\t")
        # Read MMSeqs2 GTDB assignments on contig level
        mmseqs2 = pd.read_csv(params.contigs_mmseqs2, sep="\t")
        mmseqs2 = mmseqs2.loc[~mmseqs2['lineage'].isnull()]

        # Identify most common taxunit and infer path to root
        extract_lowest_rank = re.compile(r'.*;*([a-z]_)[A-Za-z0-9_\- ]+$')
        rank_order = ['d_', 'p_', 'c_', 'o_', 'f_', 'g_', 's_']
        ranks = [extract_lowest_rank.search(l).group(1)
                for l in mmseqs2['lineage'].tolist()
                if l != "-_root"]
        lowest_rank = max([rank_order.index(r) for r in ranks])
        if lowest_rank > 4:
            lowest_rank = 4
        lineage_counts = mmseqs2.loc[mmseqs2['lineage'].str.contains(rank_order[lowest_rank])]['lineage'].value_counts()
        most_common_lineage = lineage_counts.index[0]
        ranks = most_common_lineage.split(";")
        taxpaths = [";".join(ranks[:i]) for i in range(1, len(ranks) + 1)]

        # Read GFF file for gene analysis for contigs that were assigned to kingdom or phylum by LCA
        gff = pd.DataFrame([line.rstrip().split("\t")
                            for line in open(f"{params.prokkadir}/{wildcards.bin}.gff", "rt")
                            if not line.startswith("#")])
        gff.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
        # Read MMSeqs2 results on gene level
        if os.stat(params.genes_mmseqs2).st_size > 0:
            mmseqs2_genes = pd.read_csv(params.genes_mmseqs2, sep="\t")
            mmseqs2_genes = mmseqs2_genes.loc[~mmseqs2_genes['lineage'].isnull()]

        # Filter contigs
        ## 1. Retain all contigs falling on path
        contigs_on_path = mmseqs2.loc[((mmseqs2['lineage'].isin(taxpaths)) |
                                    (mmseqs2['lineage'].str.contains(taxpaths[-1])))]
        ## 2. Remove contigs with genes from different organisms - chimeric contigs
        if os.stat(params.genes_mmseqs2).st_size > 0:
            gff['contig'] = gff['attribute'].str.split(";").str[0].str.replace("ID=", "")
            mmseqs2_genes = mmseqs2_genes.merge(gff[['seqname', 'contig']], how="left", on="contig")
            mmseqs2_genes['consensus_lineage'] = mmseqs2_genes['lineage'].isin(taxpaths)
            kp_consensus = mmseqs2_genes.groupby(['seqname'])['consensus_lineage'].agg(lambda x: x.sum() / x.count()) \
                .reset_index()
            kp_consensus = kp_consensus.rename({'seqname': 'contig'}, axis=1)
            contigs_on_path = contigs_on_path.merge(kp_consensus[['contig', 'consensus_lineage']], how="left", on="contig")
            contigs_on_path['consensus_lineage'] = contigs_on_path['consensus_lineage'].fillna(value=1) == 1.0
            contigs_on_path = contigs_on_path.loc[contigs_on_path['consensus_lineage']]
        ## 3. Remove kingdom and phylum assigned contigs if coverage > 2 std of
        ## mean of median coverage across all class, order, family, genus, and
        ## species contigs
        contigs_on_path = contigs_on_path.merge(coverage[['Contig', 'Depth_(median)']],
                                                how="left", left_on="contig", right_on="Contig")
        contigs_on_path['cofgs'] = contigs_on_path['lineage'].str.count(";") > 1
        cofgs_cov_avg = contigs_on_path.loc[contigs_on_path['cofgs']]['Depth_(median)'].mean()
        cofgs_cov_std = contigs_on_path.loc[contigs_on_path['cofgs']]['Depth_(median)'].std()
        contigs_on_path['depth_filter'] = (contigs_on_path.loc[~contigs_on_path['cofgs']]['Depth_(median)'] - cofgs_cov_avg).abs() < (2 * cofgs_cov_std)
        contigs_on_path['depth_filter'] = contigs_on_path['depth_filter'].fillna(value=True)
        contigs_on_path = contigs_on_path.drop(['Contig'], axis=1)

        # Generate final list of contigs
        contigs = set(contigs_on_path.loc[contigs_on_path['depth_filter'],
                                        'contig'].tolist())

        # Write FastA file
        with bgzf.BgzfWriter(output.fa, "wb") as outfile:
            for name, seq in pyfastx.Fasta(params.fa, build_index=False):
                if name in contigs:
                    outfile.write(str(f">{name}\n{seq}\n").encode("utf-8"))

        contigs_on_path.to_csv(output.contigs, sep="\t", index=False)
