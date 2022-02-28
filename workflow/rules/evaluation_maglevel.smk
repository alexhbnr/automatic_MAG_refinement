wildcard_constraints:
    sample = "[A-Za-z0-9\.]+",
    bin = "[A-Za-z0-9\.]+_[0-9]+"

localrules: decompress_bins

def return_bins_of_cmseq(wildcards):
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
        metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
        bins.extend(metawrap['sample_binID'].tolist())
    return [f"{config['resultdir']}/{b.split('_')[0]}/bins/{b}.{analysis}.txt" for b in bins for analysis in ['breadth_depth', 'polymut']]

def return_bins_of_phylophlan(wildcards):
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
        metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
        bins.extend(metawrap['sample_binID'].tolist())
    return [f"{config['tmpdir']}/phylophlan/{b}.fa" for b in bins]

rule evaluation_maglevel:
    input:
        phylophlan = lambda wildcards: f"{config['resultdir']}/phylophlan_closestGenomes.tsv",
        cmseq = return_bins_of_cmseq
    output:
        touch("{tmpdir}/evaluation_maglevel.done")

rule bedfile_filtered_contigs:
    input:
        lambda wildcards: f"{config['resultdir']}/{wildcards.bin.split('_')[0]}/bins/{wildcards.bin}.fasta.gz"
    output:
        temp("{tmpdir}/filtered_bins/{bin}.contiglist.txt")
    message: "Generate BED file for the filtered contigs belonging to bin {wildcards.bin} for subsetting the BAM file"
    resources:
        mem = 4,
        cores = 1
    wrapper:
        "file:workflow/wrappers/bedfile_filtered_contigs"

rule subset_bam_to_filtered_mag:
    input:
        "{tmpdir}/filtered_bins/{bin}.contiglist.txt"
    output:
        temp("{tmpdir}/filtered_bins/{bin}.mag.bam")
    message: "Subset BAM file for contigs of bin {wildcards.bin}"
    resources:
        mem = 2,
        cores = 1
    params:
        bam = lambda wildcards: sampletsv.at[wildcards.bin.split("_")[0], 'bamfn']
    wrapper:
        "file:workflow/wrappers/subset_bam_to_mag"

rule remove_extra_headerlines_filtered_contigs:
    input:
        contiglist = "{tmpdir}/filtered_bins/{bin}.contiglist.txt",
        bam = "{tmpdir}/filtered_bins/{bin}.mag.bam"
    output:
        bam = temp("{tmpdir}/filtered_bins/{bin}.bam"),
        bai = temp("{tmpdir}/filtered_bins/{bin}.bam.bai")
    message: "Remove contigs not belonging to bin to allow for faster iteration: {wildcards.bin}"
    resources:
        mem = 4,
        cores = 1
    params:
        fasta = lambda wildcards: sampletsv.at[wildcards.bin.split("_")[0], 'fastafn'],
        assembler = config['assembler']
    wrapper:
        "file:workflow/wrappers/remove_extra_headerlines"

rule breadth_depth_filtered_contigs:
    input: 
        bam = lambda wildcards: f"{config['tmpdir']}/filtered_bins/{wildcards.bin}.bam",
        bai = lambda wildcards: f"{config['tmpdir']}/filtered_bins/{wildcards.bin}.bam.bai"
    output:
        "{resultdir}/{sample}/bins/{bin}.breadth_depth.txt"
    message: "Investigate the breadth and depth of sample {wildcards.bin}"
    resources:
        mem = 4,
        cores = 1
    wrapper:
        "file:workflow/wrappers/breadth_depth"

rule prokka_filtered_mag:
    input:
        "{resultdir}/{sample}/bins/{bin}.fasta.gz"
    output:
        "{resultdir}/{sample}/bins/{bin}.gff"
    message: "Run Prokka on contigs of bin {wildcards.bin}"
    resources:
        mem = 16,
        cores = 8
    params:
        fa = "{resultdir}/{sample}/bins/{bin}.fasta.gz",
        tmpdir = lambda wildcards: f"{config['tmpdir']}/prokka_fltbin_{wildcards.bin}",
        outdir = "{resultdir}/{sample}/bins",
    threads: 8
    wrapper:
        "file:workflow/wrappers/prokka"

rule polymut_filtered_contigs:
    # Used by Pasolli et al. (2019); returns three-column table: number of
    # non-syn. mutations, number of syn. mutations, total number of positions
    input:
        bam = lambda wildcards: f"{config['tmpdir']}/filtered_bins/{wildcards.bin}.bam",
        bai = lambda wildcards: f"{config['tmpdir']}/filtered_bins/{wildcards.bin}.bam.bai",
        gff = "{resultdir}/{sample}/bins/{bin}.gff"
    output:
        "{resultdir}/{sample}/bins/{bin}.polymut.txt"
    message: "Estimate the polymorphic rate of sample {wildcards.bin}"
    resources:
        mem = 4,
        cores = 1
    wrapper:
        "file:workflow/wrappers/polymut"

rule decompress_bins:
    input:
        lambda wildcards: f"{config['resultdir']}/{wildcards.bin.split('_')[0]}/bins/{wildcards.bin}.fasta.gz"
    output:
        temp("{tmpdir}/phylophlan/{bin}.fa")
    message: "De-compress {wildcards.bin}"
    shell:
        "gunzip -c {input} > {output}"

rule phylophlan_metagenomics_closestgenomes:
    input:
        lambda wildcards: return_bins_of_phylophlan(wildcards)
    output:
        "{resultdir}/phylophlan_closestGenomes.tsv"
    message: "Run PhyloPhlAn metagenomics analysis on bins returning set of 10 closest genomes per bin"
    resources:
        mem = 100,
        cores = 36
    params:
        dbdir = lambda wildcards: f"{config['resourcesdir']}/phylophlan_databases",
        database = "SGB.Jul20",
        bindir = f"{config['tmpdir']}/phylophlan",
        prefix = f"{config['tmpdir']}/phylophlan/phylophlan_closestGenomes"
    threads: 36
    wrapper:
        "file:workflow/wrappers/phylophlan_metagenomic"
