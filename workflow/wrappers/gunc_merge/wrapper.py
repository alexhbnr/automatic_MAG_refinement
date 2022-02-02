__author__ = "Alexander Hübner"
__license__ = "MIT"

import shutil

from snakemake.shell import shell

shell(
    "(gunc merge_checkm "
    "    -g {snakemake.input.gunc} "
    "    -c {snakemake.input.checkm} "
    "    -o {snakemake.params.dir})"
)

shutil.move(f"{snakemake.params.dir}/GUNC_checkM.merged.tsv",
            snakemake.output[0])
