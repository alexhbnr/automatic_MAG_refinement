__author__ = "Alexander Hübner"
__license__ = "MIT"

import pandas as pd
import pyfastx

bincode = f"bin.{int(snakemake.wildcards.bin.split('_')[1])}"
contigs = pd.read_csv(snakemake.params.contigs, sep="\t",
                      header=None, names=['contig', 'bin'])
contigs = contigs.loc[contigs['bin'] == bincode]
contig_lengths = {name: len(seq)
                  for name, seq in pyfastx.Fasta(snakemake.params.fasta, build_index=False)
                  if name in contigs['contig'].tolist()}

# Write list of contigs for extraction
with open(snakemake.output[0], "wt") as outfile:
    for contig in contigs['contig']:
        outfile.write(f"{contig}:1-{contig_lengths[contig]}\n")
