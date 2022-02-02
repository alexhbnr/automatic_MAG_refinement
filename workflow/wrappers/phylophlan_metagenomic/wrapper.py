__author__ = "Alexander Hübner"
__license__ = "MIT"

import shutil

from snakemake.shell import shell

shell(
        "phylophlan_metagenomic "
        "    -i {snakemake.params.bindir} "
        "    --input_extension .fa "
        "    -o {snakemake.params.prefix} "
        "    --how_many 50 "
        "    --nproc {snakemake.threads} "
        "    --verbose "
        "    --database {snakemake.params.database} "
        "    --database_folder {snakemake.params.dbdir}"
)
shutil.move(f"{snakemake.params.prefix}.tsv",
            snakemake.output[0])

