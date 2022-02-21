__author__ = "Alexander Hübner"
__license__ = "MIT"

import os

import pyfastx

with open(snakemake.output[0], "wt") as outfile:
    for fn in snakemake.params.fas:
        input_empty = True
        if snakemake.params.type == "contigs":
            sample = os.path.basename(os.path.dirname(os.path.dirname(fn)))
        elif snakemake.params.type == "genes":
            sample = os.path.basename(fn).replace("-contigs_superkingdom_phylum.fa.gz", "")
        # Test whether input FastA file is empty
        if snakemake.params.type == "genes":
            with open(fn, "rb") as f:
                fstart = f.read(3)
                if fstart.startswith(b"\x1f\x8b\x08"):  # gzipped
                    if os.stat(fn).st_size > 50:
                        input_empty = False
                else:
                    if os.stat(fn).st_size > 0:
                        input_empty = False

        if not input_empty or snakemake.params.type == "contigs":
            for name, seq in pyfastx.Fasta(fn, build_index=False):
                outfile.write(f">{sample}:{name}\n{seq}\n")
