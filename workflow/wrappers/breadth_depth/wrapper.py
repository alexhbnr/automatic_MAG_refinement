__author__ = "Alexander Hübner"
__license__ = "MIT"

from snakemake.shell import shell

shell(
    "breadth_depth.py -f "
    "--minqual 30 "
    "--mincov 1 "
    "{snakemake.input.bam} > {snakemake.output[0]}"
)
