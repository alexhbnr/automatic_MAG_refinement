__author__ = "Alexander Hübner"
__license__ = "MIT"

import shutil

from snakemake.shell import shell

shutil.rmtree(snakemake.params.outputfd)

shell(
    "(checkm lineage_wf "
    "    -x fa "
    "    -t {snakemake.threads} "
    "    {snakemake.params.fadir} "
    "    {snakemake.params.outputfd} > {snakemake.log})"
)

