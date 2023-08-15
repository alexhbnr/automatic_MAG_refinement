import os
from pathlib import Path

import pandas as pd

if os.stat(snakemake.input[0]).st_size > 0:
    pd.read_csv(snakemake.input[0], sep="\t", header=None,
                usecols=list(range(8)) + [9],
                names=['contig', 'NCBItaxID', 'NCBIrank', 'NCBItaxName',
                       'nFrags', 'retainedFrags', 'taxassignedFrags',
                       'fractionAgreement', 'lineage']) \
        .to_csv(snakemake.output[0], sep="\t", index=False)
else:
    Path(snakemake.output[0]).touch()
