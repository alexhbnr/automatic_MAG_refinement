import pandas as pd

wildcard_constraints:
    bin = "[A-Za-z0-9\.]+_[0-9]+"

#### Auxilliary functions ######################################################

def path_to_bin_fa(wildcards):
    binid = f"bin.{int(wildcards.bin.split('_')[1])}"
    metawrap_path = sampletsv.at[wildcards.bin.split('_')[0], 'metawrapreport'] \
        .replace(".stats", "")
    return f"{metawrap_path}/{binid}.fa"

def return_mmseqs2_gtdb_tsv_files(wildcards):
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
        metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
        bins.extend(metawrap['sample_binID'].tolist())
    return [f"{wildcards.tmpdir}/mmseqs2/{b}.mmseqs2_gtdb.annot.tsv" for b in bins]

################################################################################

rule mmseqs2_gtdb_contigs:
    input:
        return_mmseqs2_gtdb_tsv_files
    output:
        touch("{tmpdir}/mmseqs2_gtdb_contigs.done")

rule createdb_bins:
    output:
        "{tmpdir}/mmseqs2/{bin}.contigs"
    message: "Create database of contigs of bin {wildcards.bin}"
    resources:
        mem = 16,
        cores = 1
    params:
        fa = lambda wildcards: path_to_bin_fa(wildcards)
    wrapper:
        "file:workflow/wrappers/mmseqs2_createdb"

rule screen:
    input:
        db = f"{config['resourcesdir']}/mmseqs2/gtdb/mmseqs2_gtdb_mapping",
        contigs = "{tmpdir}/mmseqs2/{bin}.contigs"
    output:
        "{tmpdir}/mmseqs2/{bin}.mmseqs2_gtdb.index"
    message: "Assign taxonomy via the GTDB for contigs of bin: {wildcards.bin}"
    resources:
        mem = 500,
        cores = 24
    params:
        tmpdir = "{tmpdir}/mmseqs2_tmpdir",
        prefix = "{tmpdir}/mmseqs2/{bin}.mmseqs2_gtdb",
        gtdb_db = f"{config['resourcesdir']}/mmseqs2/gtdb/mmseqs2_gtdb"
    threads: 24
    wrapper:
        "file:workflow/wrappers/mmseqs2_taxonomy"

rule create_tsv:
    input:
        contigs = "{tmpdir}/mmseqs2/{bin}.contigs",
        assignments = "{tmpdir}/mmseqs2/{bin}.mmseqs2_gtdb.index"
    output:
        pipe("{tmpdir}/mmseqs2/{bin}.mmseqs2_gtdb.tsv")
    message: "Convert MMSeqs2 GTDB results to TSV: {wildcards.bin}"
    resources:
        mem = 8,
        cores = 1
    params:
        prefix = "{tmpdir}/mmseqs2/{bin}.mmseqs2_gtdb",
    wrapper:
        "file:workflow/wrappers/mmseqs2_createtsv"

rule mmseqs2_annotatetsv:
    input:
        "{tmpdir}/mmseqs2/{bin}.mmseqs2_gtdb.tsv"
    output:
        "{tmpdir}/mmseqs2/{bin}.mmseqs2_gtdb.annot.tsv"
    message: "Add header to MMSeqs2 table: {wildcards.bin}"
    wrapper:
        "file:workflow/wrappers/mmseqs2_annotatetsv"
