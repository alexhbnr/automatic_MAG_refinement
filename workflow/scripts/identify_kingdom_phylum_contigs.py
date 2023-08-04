import gzip
import re

from Bio import bgzf
import pandas as pd
import pyfastx

# Read GFF file
gff = pd.DataFrame([line.rstrip().split("\t")
                    for line in open(f"{snakemake.params.prokkadir}/{snakemake.wildcards.bin}.gff", "rt")
                    if not line.startswith("#")])
gff.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']

# Read MMSeqs2 results
mmseqs2 = pd.read_csv(snakemake.params.tsv, sep="\t")
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
with open(snakemake.output.mclineage, "wt") as outfile:
    outfile.write(most_common_lineage + "\n")
ranks = most_common_lineage.split(";")
taxpaths = [";".join(ranks[:i]) for i in range(1, 3)]
contigs = set(mmseqs2.loc[mmseqs2['lineage'].isin(taxpaths)]['contig'].tolist())
genes = set(gff.loc[gff['seqname'].isin(contigs)]['attribute'].str.split(";").str[0].str.replace("ID=", "").tolist())

with bgzf.BgzfWriter(snakemake.output.fa, "wb") as outfile:
    for name, seq in pyfastx.Fasta(f"{snakemake.params.prokkadir}/{snakemake.wildcards.bin}.ffn",
                                    build_index=False):
        if name in genes:
            outfile.write(str(f">{name}\n{seq}\n").encode("utf-8"))
