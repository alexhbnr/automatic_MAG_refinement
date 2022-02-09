__author__ = "Alexander Hübner"
__license__ = "MIT"

import os

import bgzip
import pandas as pd
import pyfastx

# Read coverage
coverage = pd.read_csv(snakemake.input.depth, sep="\t")
# Read MMSeqs2 GTDB assignments on contig level
mmseqs2 = pd.read_csv(snakemake.input.contigs_mmseqs2, sep="\t")
mmseqs2 = mmseqs2.loc[~mmseqs2['lineage'].isnull()]

# Identify most common taxunit and infer path to root
lineage_counts = mmseqs2.loc[mmseqs2['lineage'].str.contains("f_")]['lineage'].value_counts()
most_common_lineage = lineage_counts.index[0]
ranks = most_common_lineage.split(";")
taxpaths = [";".join(ranks[:i]) for i in range(1, len(ranks) + 1)]

# Read GFF file for gene analysis for contigs that were assigned to kingdom or phylum by LCA
gff = pd.DataFrame([line.rstrip().split("\t")
                    for line in open(f"{snakemake.params.prokkadir}/{snakemake.wildcards.bin}.gff", "rt")
                    if line.startswith("gnl|Prokka")])
gff.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
gff = gff.loc[gff['source'].str.startswith("Prodigal")]
# Generate contig map
original_contig_names = [name for name, _ in pyfastx.Fasta(snakemake.params.fa, build_index=False)]
prokka_contig_names = [name for name, _ in pyfastx.Fasta(f"{snakemake.params.prokkadir}/{snakemake.wildcards.bin}.fna",
                                                         build_index=False)]
contig_name_map = {p: o for o, p in zip(original_contig_names, prokka_contig_names)}
# Read MMSeqs2 results on gene level
if os.stat(snakemake.params.genes_mmseqs2).st_size > 0:
    mmseqs2_genes = pd.read_csv(snakemake.params.genes_mmseqs2, sep="\t")
    mmseqs2_genes = mmseqs2_genes.loc[~mmseqs2_genes['lineage'].isnull()]

# Filter contigs
## 1. Retain all contigs falling on path
contigs_on_path = mmseqs2.loc[((mmseqs2['lineage'].isin(taxpaths)) |
                               (mmseqs2['lineage'].str.contains(taxpaths[-1])))]
## 2. Remove contigs with genes from different organisms - chimeric contigs
if os.stat(snakemake.params.genes_mmseqs2).st_size > 0:
    gff['contig'] = gff['attribute'].str.split(";").str[0].str.replace("ID=", "")
    mmseqs2_genes = mmseqs2_genes.merge(gff[['seqname', 'contig']], how="left", on="contig")
    mmseqs2_genes['consensus_lineage'] = mmseqs2_genes['lineage'].isin(taxpaths)
    kp_consensus = mmseqs2_genes.groupby(['seqname'])['consensus_lineage'].agg(lambda x: x.sum() / x.count()) \
        .reset_index()
    kp_consensus['contig'] = [contig_name_map[s] for s in kp_consensus['seqname']]
    contigs_on_path = contigs_on_path.merge(kp_consensus[['contig', 'consensus_lineage']], how="left", on="contig")
    contigs_on_path['consensus_lineage'] = contigs_on_path['consensus_lineage'].fillna(value=1) == 1.0
    contigs_on_path = contigs_on_path.loc[contigs_on_path['consensus_lineage']]
## 3. Remove kingdom and phylum assigned contigs if coverage > 2 std of
## mean of median coverage across all class, order, family, genus, and
## species contigs
contigs_on_path = contigs_on_path.merge(coverage[['Contig', 'Depth median']],
                                        how="left", left_on="contig", right_on="Contig")
contigs_on_path['cofgs'] = contigs_on_path['lineage'].str.count(";") > 1
cofgs_cov_avg = contigs_on_path.loc[contigs_on_path['cofgs']]['Depth median'].mean()
cofgs_cov_std = contigs_on_path.loc[contigs_on_path['cofgs']]['Depth median'].std()
contigs_on_path['depth_filter'] = (contigs_on_path.loc[~contigs_on_path['cofgs']]['Depth median'] - cofgs_cov_avg).abs() < (2 * cofgs_cov_std)
contigs_on_path['depth_filter'] = contigs_on_path['depth_filter'].fillna(value=True)
contigs_on_path = contigs_on_path.drop(['Contig'], axis=1)

# Generate final list of contigs
contigs = set(contigs_on_path.loc[contigs_on_path['depth_filter'],
                                  'contig'].tolist())

# Write FastA file
with open(snakemake.output.fa, "wb") as rawstream:
    with bgzip.BGZipWriter(rawstream) as outfile:
        for name, seq in pyfastx.Fasta(snakemake.params.fa, build_index=False):
            if name in contigs:
                outfile.write(str(f">{name}\n{seq}\n").encode("utf-8"))

contigs_on_path.to_csv(snakemake.output.contigs, sep="\t", index=False)

