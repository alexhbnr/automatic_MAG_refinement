__author__ = "Alexander Hübner"
__license__ = "MIT"

from snakemake.shell import shell

shell(
    "(cat {snakemake.input[0]} | xargs samtools view -hu {snakemake.params.bam} > {snakemake.output[0]})"
)
