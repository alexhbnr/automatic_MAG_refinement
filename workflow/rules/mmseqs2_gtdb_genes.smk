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

rule mmseqs2_gtdb_genes:
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
        type = "genes",
        fas = lambda wildcards: return_genes_fas(wildcards)
    wrapper:
        "file:workflow/wrappers/concat_bins"

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
        pipe("{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes.tsv")
    message: "Convert MMSeqs2 GTDB results to TSV"
    container: "https://depot.galaxyproject.org/singularity/mmseqs2:14.7e284--pl5321hf1761c0_0"
    resources:
        mem = 8,
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
    wrapper:
        "file:workflow/wrappers/mmseqs2_annotatetsv"

rule mmseqs2_splittsv_genes:
    input:
        "{tmpdir}/mmseqs2_genes/kp_contigs.mmseqs2_gtdb_genes.annot.tsv"
    output:
        touch("{tmpdir}/mmseqs2_genes/mmseqs2_splittsv.done")
    message: "Split MMSeqs2 result table into bins"
    params:
        sampletsv = config['sampletsv'],
        dir = "{tmpdir}/mmseqs2_genes"
    wrapper:
        "file:workflow/wrappers/mmseqs2_splittsv_genes"

rule filter_contigs:
    input:
        genes_mmseqs2 = lambda wildcards: f"{config['tmpdir']}/mmseqs2_genes/mmseqs2_splittsv.done",
        depth = lambda wildcards: f"{config['tmpdir']}/depth/{wildcards.bin}.breadth_depth.txt"
    output:
        fa = "{resultdir}/{bin}/{bin}.fasta.gz",
        contigs = "{resultdir}/{bin}/{bin}.contigs.txt"
    message: "Automatically filter contigs that have a high chance to not belong to major lineage: {wildcards.bin} for sample {wildcards.sample}"
    params:
        fa = lambda wildcards: path_to_bin_fa(wildcards),
        contigs_mmseqs2 = lambda wildcards: f"{config['tmpdir']}/mmseqs2/{wildcards.bin}.mmseqs2_gtdb.annot.tsv",
        genes_mmseqs2 = lambda wildcards: f"{config['tmpdir']}/mmseqs2_genes/{wildcards.bin}.mmseqs2_gtdb.annot.tsv",
        prokkadir = lambda wildcards: f"{config['tmpdir']}/pyrodigal"
    wrapper:
        "file:workflow/wrappers/filter_contigs"
