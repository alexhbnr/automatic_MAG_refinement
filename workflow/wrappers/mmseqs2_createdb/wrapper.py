__author__ = "Alexander Hübner"
__license__ = "MIT"

from snakemake.shell import shell

shell(
    "(mmseqs createdb {snakemake.input[0]} {snakemake.output[0]})"
)
