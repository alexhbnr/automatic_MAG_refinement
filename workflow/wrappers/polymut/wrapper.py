__author__ = "Alexander Hübner"
__license__ = "MIT"

from snakemake.shell import shell

shell(
    "polymut.py -f "
    "--gff_file {snakemake.input.gff} "
    "--minqual 30 "
    "--mincov 5 "
    "--dominant_frq_thrsh 0.8 "
    "{snakemake.input.bam} > {snakemake.output[0]}"
)
