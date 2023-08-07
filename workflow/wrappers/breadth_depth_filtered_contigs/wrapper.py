__author__ = "Alexander Hübner"
__license__ = "MIT"

import pandas as pd
import pyfastx

# Retained contigs
contigs = [name for name, _ in pyfastx.Fastx(snakemake.input.fasta)]

# Read depth
pd.read_csv(snakemake.input.depth, sep="\t") \
    .query('Contig.isin(@contigs)') \
    .to_csv(snakemake.output[0], sep="\t", index=False, float_format="%.1f")
