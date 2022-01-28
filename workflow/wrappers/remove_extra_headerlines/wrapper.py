__author__ = "Alexander Hübner"
__license__ = "MIT"

import re

import pyfastx
import pysam

# Load BAM file
bamfile = pysam.AlignmentFile(snakemake.input.bam)
header = bamfile.header

# Generate subset of SQ tags from contig list
extract_sn_ln = re.compile(r'(k[0-9]+_[0-9]+):1-([0-9]+)')
with open(snakemake.input.contiglist, "rt") as contigfile:
    sq_info = [extract_sn_ln.search(line.rstrip()).groups()
               for line in contigfile]
header_sq = [{'SN': contig[0], 'LN': int(contig[1])}
             for contig in sq_info]

# Add contigs of mates aligning to different contig than in bin to header
contig_ids = [contig[0] for contig in sq_info]
contig_lengths = {name: len(seq)
                  for name, seq in pyfastx.Fasta(snakemake.params.fasta, build_index=False)}
for read in bamfile:
    if read.is_paired:
        if read.next_reference_name not in contig_ids and read.next_reference_name is not None:
            contig_ids.append(read.next_reference_name)
            header_sq.append({'SN': read.next_reference_name, 'LN': contig_lengths[read.next_reference_name]})
header['SQ'] = header_sq
bamfile.close()

# Write header and reads back to file and fix read referenceID
with pysam.AlignmentFile(snakemake.input.bam) as bamfile:
    with pysam.AlignmentFile(snakemake.output.bam, "wbu", header=header) as outfile:
        for read in bamfile:
            read.reference_id = contig_ids.index(read.reference_name)
            if read.is_paired and read.next_reference_name is not None:
                read.next_reference_id = contig_ids.index(read.next_reference_name)
            outfile.write(read)
pysam.index(snakemake.output.bam)
