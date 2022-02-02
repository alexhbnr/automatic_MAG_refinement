import pandas as pd

wildcard_constraints:
    sample = "[A-Za-z0-9\.]+",
    bin = "[A-Za-z0-9\.]+_[0-9]+"

#### Auxilliary functions ######################################################

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
    return [f"{config['resultdir']}/{b.split('_')[0]}/bins/{b}.fasta.gz" for b in bins]

################################################################################

rule mmseqs2_gtdb_genes:
    input:
        return_filter_fasta_files
    output:
        touch(f"{config['tmpdir']}/filter_contigs.done")

rule prokka_mag:
    output:
        "{tmpdir}/prokka/{bin}.gff"
    message: "Run Prokka on contigs of bin {wildcards.bin}"
    resources:
        mem = 16,
        cores = 8
    params:
        fa = lambda wildcards: path_to_bin_fa(wildcards),
        tmpdir = "{tmpdir}/prokka_bin_{bin}",
        outdir = "{tmpdir}/prokka",
    threads: 8
    wrapper:
        "file:workflow/wrappers/prokka"

rule identify_kingdom_phylum_contigs:
    input:
        tsv = "{tmpdir}/mmseqs2/{bin}.mmseqs2_gtdb.annot.tsv",
        prokka = "{tmpdir}/prokka/{bin}.gff"
    output:
        fa = "{tmpdir}/mmseqs2_genes/{bin}-contigs_superkingdom_phylum.fa.gz",
        mclineage = "{tmpdir}/mmseqs2_genes/{bin}-most_common_lineage.txt"
    message: "Identify contigs that could only be assigned to taxonomic rank superkingdom or phylum: {wildcards.bin}"
    resources:
        mem = 4,
        cores = 1
    params:
        fa = lambda wildcards: path_to_bin_fa(wildcards),
        prokkadir = "{tmpdir}/prokka"
    script:
        "../scripts/identify_kingdom_phylum_contigs.py"

rule createdb_bins_kp_contigs:
    input:
        "{tmpdir}/mmseqs2_genes/{bin}-contigs_superkingdom_phylum.fa.gz"
    output:
        "{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_genes"
    message: "Create database of genes of bin {wildcards.bin}"
    resources:
        mem = 16,
        cores = 1
    params:
        fa = "{tmpdir}/mmseqs2_genes/{bin}-contigs_superkingdom_phylum.fa.gz"
    wrapper:
        "file:workflow/wrappers/mmseqs2_createdb"

rule screen_kp_contigs:
    input:
        db = f"{config['resourcesdir']}/mmseqs2/gtdb/mmseqs2_gtdb_mapping",
        contigs = "{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_genes"
    output:
        "{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_gtdb_genes.index"
    message: "Assign taxonomy via the GTDB for genes of bin: {wildcards.bin}"
    resources:
        mem = 500,
        cores = 24
    params:
        tmpdir = "{tmpdir}/mmseqs2_tmpdir",
        prefix = "{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_gtdb_genes",
        gtdb_db = f"{config['resourcesdir']}/mmseqs2/gtdb/mmseqs2_gtdb"
    threads: 24
    wrapper:
        "file:workflow/wrappers/mmseqs2_taxonomy"

rule create_tsv_kp_contigs:
    input:
        contigs = "{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_genes",
        assignments = "{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_gtdb_genes.index"
    output:
        pipe("{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_gtdb_genes.tsv")
    message: "Convert MMSeqs2 GTDB results to TSV: {wildcards.bin}"
    resources:
        mem = 8,
        cores = 1
    params:
        prefix = "{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_gtdb_genes",
    wrapper:
        "file:workflow/wrappers/mmseqs2_createtsv"

rule annotate_tsv_kp_contigs:
    input:
        "{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_gtdb_genes.tsv"
    output:
        "{tmpdir}/mmseqs2_genes/{bin}.mmseqs2_gtdb_genes.annot.tsv"
    message: "Add header to MMSeqs2 table: {wildcards.bin}"
    wrapper:
        "file:workflow/wrappers/mmseqs2_annotatetsv"

rule filter_contigs:
    input:
        contigs_mmseqs2 = lambda wildcards: f"{config['tmpdir']}/mmseqs2/{wildcards.bin}.mmseqs2_gtdb.annot.tsv",
        genes_mmseqs2 = lambda wildcards: f"{config['tmpdir']}/mmseqs2_genes/{wildcards.bin}.mmseqs2_gtdb_genes.annot.tsv",
        depth = lambda wildcards: f"{config['tmpdir']}/depth/{wildcards.bin}.breadth_depth.txt"
    output:
        fa = "{resultdir}/{sample}/bins/{bin}.fasta.gz",
        contigs = "{resultdir}/{sample}/bins/{bin}.contigs.txt"
    message: "Automatically filter contigs that have a high chance to not belong to major lineage: {wildcards.bin} for sample {wildcards.sample}"
    params:
        fa = lambda wildcards: path_to_bin_fa(wildcards),
        prokkadir = lambda wildcards: f"{config['tmpdir']}/prokka"
    wrapper:
        "file:workflow/wrappers/filter_contigs"
