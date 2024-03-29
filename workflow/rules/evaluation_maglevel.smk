def return_bins_of_cmseq(wildcards):
    with open(checkpoints.mmseqs2_gtdb_genes.get(**wildcards).output[0], "rt") as dummyfile:
        pass
    bins = []
    for sample in SAMPLES:
        metawrap = pd.read_csv(sampletsv.at[sample, 'metawrapreport'],
                            sep="\t")
        metawrap['binid'] = metawrap['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
        metawrap['sample_binID'] = [f"{sample}_{i:03}" for i in metawrap['binid']]
        bins.extend([b for b in metawrap['sample_binID'].tolist()
                     if os.stat(f"{config['resultdir']}/{b}/{b}.fasta.gz").st_size > 50])
    return [f"{config['resultdir']}/{b}/{b}.{analysis}.txt" for b in bins for analysis in ['breadth_depth', 'polymut']]

def sampling_fraction(wildcards):
    depth = pd.read_csv(checkpoints.breadth_depth_filtered_contigs.get(resultdir=config['resultdir'],
                                                                       bin=wildcards.bin).output[0], sep="\t")
    avg_depth = depth['Depth_(avg)'].mean()
    if avg_depth >= config['polymut_depth']:
        return f"-s 1.{str((50 / avg_depth).round(8))[2:]}"
    else:
        return ""

rule evaluation_maglevel:
    input:
        cmseq = return_bins_of_cmseq
    output:
        touch("{tmpdir}/evaluation_maglevel.done")

#### Annotate with Bakta #######################################################

rule bakta:
    input:
        db = f"{config['resourcesdir']}/bakta/downloaded",
        fasta = "{resultdir}/{bin}/{bin}.fasta.gz"
    output:
        "{resultdir}/{bin}/{bin}.gff3.gz"
    message: "Annotate contigs using BAKTA: {wildcards.bin}"
    container: "docker://quay.io/biocontainers/bakta:1.8.1--pyhdfd78af_0"
    resources:
        mem = 128,
        cores = 8
    params:
        prefix = lambda wildcards: wildcards.bin,
        outdir = "{resultdir}/{bin}",
        dbdir = f"{config['resourcesdir']}/bakta/db",
        extra = "--keep-contig-headers --meta --skip-crispr"
    threads: 8
    shell:
        """
        bakta -p {params.prefix} \
            --db {params.dbdir} \
            --force \
            --output {params.outdir} \
            {params.extra} \
            --threads {threads} \
            {input.fasta}
        for sfx in embl faa ffn fna gbff gff3 hypotheticals.faa hypotheticals.tsv json log tsv; do
            gzip -f {params.outdir}/{params.prefix}.${{sfx}}
        done
        """

################################################################################

#### Breadth and depth #########################################################

checkpoint breadth_depth_filtered_contigs:
    input: 
        depth = lambda wildcards: f"{config['tmpdir']}/depth/{wildcards.bin}.breadth_depth.txt",
        fasta = "{resultdir}/{bin}/{bin}.fasta.gz"
    output:
        "{resultdir}/{bin}/{bin}.breadth_depth.txt"
    message: "Investigate the breadth and depth of sample {wildcards.bin}"
    resources:
        mem = 4,
        cores = 1
    run:
        # Retained contigs
        if os.stat(input.fasta).st_size > 50:
            contigs = [name for name, _ in pyfastx.Fastx(input.fasta)]
        else:
            contigs = []
        # Subset depth table
        pd.read_csv(input.depth, sep="\t") \
            .query('Contig.isin(@contigs)') \
            .to_csv(output[0], sep="\t", index=False, float_format="%.1f")

################################################################################

rule bedfile_filtered_contigs:
    input:
        fasta = lambda wildcards: f"{config['resultdir']}/{wildcards.bin}/{wildcards.bin}.fasta.gz"
    output:
        temp("{tmpdir}/polymut/{bin}.contiglist.txt")
    message: "Generate BED file for contigs belonging to bin {wildcards.bin} for subsetting the BAM file"
    group: "prep_bamfile_polymut"
    resources:
        mem = 2,
        cores = 1
    run:
        contig_lengths = {name: len(seq)
                          for name, seq in pyfastx.Fasta(input[0], build_index=False)}

        # Write list of contigs for extraction
        with open(output[0], "wt") as outfile:
            for contig in contig_lengths:
                outfile.write(f"{contig}:1-{contig_lengths[contig]}\n")

rule subsample_bam_of_mag:
    input:
        contigs = "{tmpdir}/polymut/{bin}.contiglist.txt",
        depth = lambda wildcards: f"{config['resultdir']}/{wildcards.bin}/{wildcards.bin}.breadth_depth.txt"
    output:
        temp("{tmpdir}/polymut/{bin}.mag.bam")
    message: "Subset BAM file for contigs of bin {wildcards.bin}"
    group: "prep_bamfile_polymut"
    container: "https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0"
    resources:
        mem = 2,
        cores = 0
    params:
        bam = lambda wildcards: sampletsv.at[wildcards.bin.split("_")[0], 'bamfn'],
        subsampling_fraction = lambda wildcards: sampling_fraction(wildcards)
    shell:
        """
        cat {input.contigs} | xargs samtools view -o {output} {params.subsampling_fraction} -hu {params.bam}
        """

rule remove_extra_headerlines_polymut:
    input:
        contiglist = "{tmpdir}/polymut/{bin}.contiglist.txt",
        bam = "{tmpdir}/polymut/{bin}.mag.bam"
    output:
        bam = temp("{tmpdir}/polymut/{bin}.bam"),
        bai = temp("{tmpdir}/polymut/{bin}.bam.bai")
    message: "Remove contigs not belonging to bin to allow for faster iteration: {wildcards.bin}"
    group: "prep_bamfile_polymut"
    resources:
        mem = 2,
        cores = 0
    params:
        assembler = config['assembler']
    wrapper:
        "file:workflow/wrappers/remove_extra_headerlines"

rule polymut:
    # Used by Pasolli et al. (2019); returns three-column table: number of
    # non-syn. mutations, number of syn. mutations, total number of positions
    input:
        bam = lambda wildcards: f"{config['tmpdir']}/polymut/{wildcards.bin}.bam",
        bai = lambda wildcards: f"{config['tmpdir']}/polymut/{wildcards.bin}.bam.bai",
        gff = "{resultdir}/{bin}/{bin}.gff3.gz"
    output:
        "{resultdir}/{bin}/{bin}.polymut.txt"
    message: "Estimate the polymorphic rate of sample {wildcards.bin}"
    container: "https://depot.galaxyproject.org/singularity/cmseq:1.0--pyh5ca1d4c_0"
    group: "polymut"
    resources:
        mem = 4,
        cores = 1
    shell:
        """
        gunzip -c {input.gff} | \
        polymut.py -f \
            --gff_file /dev/stdin \
            --minqual 30 \
            --mincov 5 \
            --dominant_frq_thrsh 0.8 \
            {input.bam} > {output}
        """

################################################################################
