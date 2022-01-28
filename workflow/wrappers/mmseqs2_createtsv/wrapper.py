__author__ = "Alexander Hübner"
__license__ = "MIT"

import os
from pathlib import Path

from snakemake.shell import shell

if os.stat(snakemake.input.contigs).st_size > 0:
    shell(
        "(mmseqs createtsv {snakemake.input.contigs} {snakemake.params.prefix} {snakemake.output[0]})"
    )
else:
    Path(snakemake.output[0]).touch()
