__author__ = "Alexander Hübner"
__license__ = "MIT"

from snakemake.shell import shell

shell(
    "(checkm qa -o 2 --tab_table -f {snakemake.output[0]} "
    "    {snakemake.params.outputfd}/lineage.ms "
    "    {snakemake.params.outputfd})"
)

