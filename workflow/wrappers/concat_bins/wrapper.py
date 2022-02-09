__author__ = "Alexander Hübner"
__license__ = "MIT"

import os

import pyfastx

with open(snakemake.output[0], "wt") as outfile:
    for fn in snakemake.params.fas:
        sample = os.path.basename(os.path.dirname(os.path.dirname(fn)))
        for name, seq in pyfastx.Fasta(fn, build_index=False):
            outfile.write(f">{sample}:{name}\n{seq}\n")
