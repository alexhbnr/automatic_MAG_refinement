__author__ = "Alexander Hübner"
__license__ = "MIT"

import pyfastx

contig_lengths = {name: len(seq)
                  for name, seq in pyfastx.Fasta(snakemake.input[0], build_index=False)}

# Write list of contigs for extraction
with open(snakemake.output[0], "wt") as outfile:
    for contig in contig_lengths:
        outfile.write(f"{contig}:1-{contig_lengths[contig]}\n")
