import gzip

import bgzip
import pandas as pd
import pyfastx

# Read GFF file
gff = pd.DataFrame([line.rstrip().split("\t")
                    for line in open(f"{snakemake.params.prokkadir}/{snakemake.wildcards.bin}.gff", "rt")
                    if line.startswith("gnl|Prokka")])
gff.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
gff = gff.loc[gff['source'].str.startswith("Prodigal")]

# Generate contig map
original_contig_names = [name for name, _ in pyfastx.Fasta(snakemake.params.fa, build_index=False)]
prokka_contig_names = [name for name, _ in pyfastx.Fasta(f"{snakemake.params.prokkadir}/{snakemake.wildcards.bin}.fna",
                                                         build_index=False)]
contig_name_map = {o: p for o, p in zip(original_contig_names, prokka_contig_names)}

# Read MMSeqs2 results
mmseqs2 = pd.read_csv(snakemake.input[0], sep="\t")
mmseqs2 = mmseqs2.loc[~mmseqs2['lineage'].isnull()]

# Identify most common taxunit and infer path to root
lineage_counts = mmseqs2.loc[mmseqs2['lineage'].str.contains("f_")]['lineage'].value_counts()
most_common_lineage = lineage_counts.index[0]
with open(snakemake.output.mclineage, "wt") as outfile:
    outfile.write(most_common_lineage + "\n")
ranks = most_common_lineage.split(";")
taxpaths = [";".join(ranks[:i]) for i in range(1, 3)]
contigs = set([contig_name_map[c]
               for c in mmseqs2.loc[mmseqs2['lineage'].isin(taxpaths)]['contig']])
genes = set(gff.loc[gff['seqname'].isin(contigs)]['attribute'].str.split(";").str[0].str.replace("ID=", "").tolist())

with open(snakemake.output.fa, "wb") as rawstream:
    with bgzip.BGZipWriter(rawstream) as outfile:
        for name, seq in pyfastx.Fasta(f"{snakemake.params.prokkadir}/{snakemake.wildcards.bin}.ffn",
                                       build_index=False):
            if name in genes:
                outfile.write(str(f">{name}\n{seq}\n").encode("utf-8"))
