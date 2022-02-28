__author__ = "Alexander Hübner"
__license__ = "MIT"

import os

import pandas as pd

# Import MMSeqs2 results
mmseqs2 = pd.read_csv(snakemake.input[0], sep="\t")
mmseqs2['sample'] = mmseqs2['contig'].str.split(":").str[0]
mmseqs2['contig'] = mmseqs2['contig'].str.split(":").str[1]

# Import sample TSV
sampletsv = pd.read_csv(snakemake.params.sampletsv, sep="\t", index_col=[0])
sampletsv.index = sampletsv.index.astype(str)

for sample in mmseqs2['sample'].unique():
    # Load overview list of bins
    sample_bins = pd.read_csv(sampletsv.at[sample, 'metawrapreport'], sep="\t")
    sample_bins['binid'] = sample_bins['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
    sample_bins['sample_binID'] = [f"{sample}_{i:03}" for i in sample_bins['binid']]
    sample_bins = sample_bins.set_index(['sample_binID'])
    # Load map of contigs per bin
    contigs = pd.read_csv(sampletsv.at[sample, 'metawrapreport'].replace(".stats", ".contigs"),
                          sep="\t", header=None, names=['contig', 'bin'], index_col=['bin'])
    for b in sample_bins.index:
        c = contigs.loc[sample_bins.at[b, 'bin']]['contig'].tolist()
        mmseqs2.loc[(mmseqs2['sample'] == sample) & (mmseqs2['contig'].isin(c))] \
            .drop(['sample'], axis=1) \
            .to_csv(f"{snakemake.params.dir}/{b}.mmseqs2_gtdb.annot.tsv", sep="\t", index=False)
