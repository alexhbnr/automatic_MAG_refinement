__author__ = "Alexander Hübner"
__license__ = "MIT"

import os
from pathlib import Path

from snakemake.shell import shell

if os.stat(snakemake.input.contigs).st_size > 0:
    shell(
        "(mmseqs taxonomy "
        " {snakemake.input.contigs} "
        " {snakemake.params.gtdb_db} "
        " {snakemake.params.prefix} "
        " {snakemake.params.tmpdir} "
        " -a "
        " --tax-lineage 1 "
        " --lca-ranks kingdom,phylum,class,order,family,genus,species "
        " --majority 0.5 "
        " --vote-mode 1 "
        " --orf-filter 1 "
        " --remove-tmp-files 1 "
        " --threads {snakemake.threads})"
    )
else:
    Path(snakemake.output[0]).touch()
