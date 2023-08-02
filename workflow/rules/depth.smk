#### Auxilliary functions ######################################################

def return_bins_of_samples(wildcards):
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
        metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
        bins.extend(metawrap['sample_binID'].tolist())
    return [f"{wildcards.tmpdir}/depth/{b}.breadth_depth.txt" for b in bins]

################################################################################

rule depth:
    input:
        f"{config['tmpdir']}/depth_calculation.done"

rule depth_calculation:
    input:
        return_bins_of_samples
    output:
        touch("{tmpdir}/depth_calculation.done")

rule bedfile_contigs:
    output:
        temp("{tmpdir}/depth/{bin}.contiglist.txt")
    message: "Generate BED file for contigs belonging to bin {wildcards.bin} for subsetting the BAM file"
    group: "prep_bamfile_depth"
    resources:
        mem = 2,
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
    group: "prep_bamfile_depth"
    container: "https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0"
    resources:
        mem = 2,
        cores = 0
    params:
        bam = lambda wildcards: sampletsv.at[wildcards.bin.split("_")[0], 'bamfn']
    shell:
        """
        cat {input} | xargs samtools view -o {output} -hu {params.bam}
        """

rule remove_extra_headerlines:
    input:
        contiglist = "{tmpdir}/depth/{bin}.contiglist.txt",
        bam = "{tmpdir}/depth/{bin}.mag.bam"
    output:
        bam = temp("{tmpdir}/depth/{bin}.bam"),
        bai = temp("{tmpdir}/depth/{bin}.bam.bai")
    message: "Remove contigs not belonging to bin to allow for faster iteration: {wildcards.bin}"
    group: "prep_bamfile_depth"
    resources:
        mem = 2,
        cores = 0
    params:
        assembler = config['assembler']
    wrapper:
        "file:workflow/wrappers/remove_extra_headerlines"

rule breadth_depth:
    input: 
        bam = "{tmpdir}/depth/{bin}.bam",
        bai = "{tmpdir}/depth/{bin}.bam.bai"
    output:
        "{tmpdir}/depth/{bin}.breadth_depth.txt"
    message: "Investigate the breadth and depth of sample {wildcards.bin}"
    container: "https://depot.galaxyproject.org/singularity/cmseq:1.0--pyh5ca1d4c_0"
    resources:
        mem = 4,
        cores = 1
    shell:
        """
        breadth_depth.py -f \
            --minqual 30 \
            --mincov 1 \
            {input.bam} > {output}
        """

