__author__ = "Alexander Hübner"
__license__ = "MIT"

import os
import shutil

from snakemake.shell import shell

os.makedirs(snakemake.params.tmpdir, exist_ok=True)

with open(snakemake.params.fa, "rb") as f:
    fstart = f.read(3)
    if fstart.startswith(b"\x1f\x8b\x08"):  # gzipped
        shell("gunzip -c {snakemake.params.fa} > {snakemake.params.tmpdir}/{snakemake.wildcards.bin}.fa")
        fafilename = f"{snakemake.params.tmpdir}/{snakemake.wildcards.bin}.fa"
    else:
        fafilename = snakemake.params.fa

shell(
    "(prokka --outdir {snakemake.params.tmpdir}  "
    "  --prefix {snakemake.wildcards.bin}  "
    "  --force  "
    "  --compliant  "
    "  --metagenome  "
    "  --cpus {snakemake.threads}  "
    "  --debug  "
    "  {fafilename})"
)

for suffix in ['faa', 'ffn', 'fna', 'gbk', 'gff', 'tsv', 'txt']:
    shutil.copy(f"{snakemake.params.tmpdir}/{snakemake.wildcards.bin}.{suffix}",
                f"{snakemake.params.outdir}/{snakemake.wildcards.bin}.{suffix}")
shutil.rmtree(snakemake.params.tmpdir)
