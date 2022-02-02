__author__ = "Alexander Hübner"
__license__ = "MIT"

import os

from snakemake.shell import shell

os.makedirs(snakemake.params.tmpdir, exist_ok=True)

shell(
    "(gunc run "
    "    -r {snakemake.params.db_file} "
    "    --input_dir {snakemake.params.fadir} "
    "    -t {snakemake.threads} "
    "    -o {snakemake.params.outdir} "
    "    --temp_dir {snakemake.params.tmpdir} "
    "    --detailed_output)"
)

os.rmdir(snakemake.params.tmpdir)
