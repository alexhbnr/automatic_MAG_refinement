import pandas as pd

wildcard_constraints:
    bin = "[A-Za-z0-9\.]+_[0-9]+"

def return_bins_of_samples(wildcards):
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
        metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
        bins.extend(metawrap['sample_binID'].tolist())
    return [f"{wildcards.tmpdir}/depth/{b}.breadth_depth.txt" for b in bins]

rule depth_calculation:
    input:
        return_bins_of_samples
    output:
        touch("{tmpdir}/depth_calculation.done")

rule bedfile_contigs:
    output:
        temp("{tmpdir}/depth/{bin}.contiglist.txt")
    message: "Generate BED file for contigs belonging to bin {wildcards.bin} for subsetting the BAM file"
    resources:
        mem = 4,
        cores = 1
    params:
        contigs = lambda wildcards: sampletsv.at[wildcards.bin.split("_")[0], 'metawrapreport'].replace("stats", "contigs"),
        fasta = lambda wildcards: sampletsv.at[wildcards.bin.split("_")[0], 'fastafn']
    wrapper:
        "file:workflow/wrappers/bedfile_contigs_of_bin"

rule subset_bam_to_mag:
    input:
        "{tmpdir}/depth/{bin}.contiglist.txt"
    output:
        temp("{tmpdir}/depth/{bin}.mag.bam")
    message: "Subset BAM file for contigs of bin {wildcards.bin}"
    resources:
        mem = 2,
        cores = 1
    params:
        bam = lambda wildcards: sampletsv.at[wildcards.bin.split("_")[0], 'bamfn']
    wrapper:
        "file:workflow/wrappers/subset_bam_to_mag"

rule remove_extra_headerlines:
    input:
        contiglist = "{tmpdir}/depth/{bin}.contiglist.txt",
        bam = "{tmpdir}/depth/{bin}.mag.bam"
    output:
        bam = temp("{tmpdir}/depth/{bin}.bam"),
        bai = temp("{tmpdir}/depth/{bin}.bam.bai")
    message: "Remove contigs not belonging to bin to allow for faster iteration: {wildcards.bin}"
    resources:
        mem = 4,
        cores = 1
    params:
        fasta = lambda wildcards: sampletsv.at[wildcards.bin.split("_")[0], 'fastafn']
    wrapper:
        "file:workflow/wrappers/remove_extra_headerlines"

rule breadth_depth:
    input: 
        bam = "{tmpdir}/depth/{bin}.bam",
        bai = "{tmpdir}/depth/{bin}.bam.bai"
    output:
        "{tmpdir}/depth/{bin}.breadth_depth.txt"
    message: "Investigate the breadth and depth of sample {wildcards.bin}"
    resources:
        mem = 4,
        cores = 1
    wrapper:
        "file:workflow/wrappers/breadth_depth"

