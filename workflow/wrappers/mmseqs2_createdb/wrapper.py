__author__ = "Alexander Hübner"
__license__ = "MIT"

import os
from pathlib import Path

from snakemake.shell import shell

# Test whether input FastA file is empty
input_empty = True
with open(snakemake.params.fa, "rb") as f:
    fstart = f.read(3)
    if fstart.startswith(b"\x1f\x8b\x08"):  # gzipped
        if os.stat(snakemake.params.fa).st_size > 50:
            input_empty = False
    else:
        if os.stat(snakemake.params.fa).st_size > 0:
            input_empty = False

if not input_empty:
    shell(
        "(mmseqs createdb {snakemake.params.fa} {snakemake.output[0]})"
    )
else:
    Path(snakemake.output[0]).touch()