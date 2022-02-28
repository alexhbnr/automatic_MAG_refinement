__author__ = "Alexander Hübner"
__license__ = "MIT"

import os
from pathlib import Path

import pandas as pd

# Import MMSeqs2 results
mmseqs2 = pd.read_csv(snakemake.input[0], sep="\t")
mmseqs2['sample_binID'] = mmseqs2['contig'].str.split(":").str[0]
mmseqs2['sample'] = mmseqs2['sample_binID'].str.split("_").str[0]
mmseqs2['contig'] = mmseqs2['contig'].str.split(":").str[1]

# Generate list of all expected bins
sampletsv = pd.read_csv(snakemake.params.sampletsv, sep="\t", index_col=[0])
expected_bins = []
for sample in sampletsv.index:
    # Load overview list of bins
    sample_bins = pd.read_csv(sampletsv.at[sample, 'metawrapreport'], sep="\t")
    sample_bins['binid'] = sample_bins['bin'].str.extract(r'bin.([0-9]+)$').astype(int)
    sample_bins['sample_binID'] = [f"{sample}_{i:03}" for i in sample_bins['binid']]
    expected_bins.extend(sample_bins['sample_binID'].tolist())

# Iterate overall bins: if genes were screened, export results, others create
# an empty file
for binid in expected_bins:
    if binid in mmseqs2['sample_binID'].tolist():
        mmseqs2.loc[mmseqs2['sample_binID'] == binid] \
            .drop(['sample_binID'], axis=1) \
            .to_csv(f"{snakemake.params.dir}/{binid}.mmseqs2_gtdb.annot.tsv", sep="\t", index=False)
    else:
        Path(f"{snakemake.params.dir}/{binid}.mmseqs2_gtdb.annot.tsv").touch()
