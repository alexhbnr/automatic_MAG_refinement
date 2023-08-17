Automatic revision of MAGs to remove likely chimeric contigs
================

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#quick-start" id="toc-quick-start">Quick start</a>
- <a href="#why-this-pipeline-matters"
  id="toc-why-this-pipeline-matters">Why this pipeline matters</a>
- <a href="#references" id="toc-references">References</a>

## Introduction

This is basic pipeline to automatically identify and remove contigs from
metagenome-assembled genomes (MAGs) that have a high chance of being
chimeric and do not actually belong to the MAG. It is written in
[Snakemake](https://snakemake.readthedocs.io/) and makes use of
containers
([Singularity](https://docs.sylabs.io/guides/3.5/user-guide/index.html)
and [Docker](https://www.docker.com/)) and [conda
environments](https://docs.conda.io/en/latest/) to ensure its
reproducibility. It has been previously been used by Klapper et al.
(2023) and has received some further updates since then.

The pipeline implements three different steps:

1.  Taxonomic assignment of all contigs of a MAG using MMseqs2’s
    taxonomy workflow (Mirdita et al. 2021) against the GTDB and
    identifying of the most common lineage at the rank family or higher
2.  Taxonomic assignment of every gene of a contig that could only be
    assigned to the taxonomic rank kingdom or phylum using MMseqs2’s
    taxonomy workflow (Mirdita et al. 2021) against the GTDB to
    discriminate contigs carrying preserved genes from contigs carrying
    genes from multiple lineages
3.  Evaluation of the average depth of all contigs

Afterwards, the information of these three steps are combined and
contigs are removed because they are likely chimeric when:

1.  A contig was assigned to a lineage at the taxonomic level class or
    lower that does not overlap with the most common lineage of the MAG.
    E.g. the main lineage was
    `d_Bacteria;p_Firmicutes A;c_Clostridia;o_Lachnospirales;f_Lachnospiraceae`
    but a contig was assigned to
    `d_Bacteria;p_Firmicutes      A;c_Clostridia;o_Peptostreptococcales;f_Filifactoraceae;g_Peptoanaerobacter;s_Peptoanaerobacter      stomatis`.
2.  A contig could only be assigned to the taxomic rank kingdom and
    phylum and the per-gene analysis revealed that genes on this contig
    can be assigned to a lineage other than the main lineage.
3.  A contig could only be assigned to the taxomic rank kingdom and
    phylum and its coverage was deviating from the average coverage of
    all contigs assigned to the main lineage by more than two standard
    deviations.

Finally, the pipeline runs multiple quality evaluation steps on the
refined MAGs:

1.  Functional annotation using Bakta (Schwengers et al. n.d.)
2.  Evaluation of the quality of the MAG using - checkM (Parks et
    al. 2015) - checkM2 (Chklovski et al. 2023) - GUNC (Orakov et
    al. 2021) - the presence of alternative alleles with a minimal
    allele frequency of 20% within coding sequences that would lead to a
    non-synonymous mutation (Pasolli et al. 2019)
3.  Taxonomic assignment using - GTDBTK (Chaumeil et al. 2020) -
    PhyloPhlAn (Asnicar et al. 2020)

A detailed description of the methodology and some results can be found
in the [Supplementary
Material](https://www.science.org/doi/suppl/10.1126/science.adf5300/suppl_file/science.adf5300_sm.pdf)
Section 3 “Reference-free binning of the *de novo* assembled contigs” of
Klapper et al. (2023).

## Quick start

To be able to run the pipeline,
[Snakemake](https://snakemake.readthedocs.io/) with a minimal version of
7.0 is necessary. The easiest way to install the dependencies of this
program and to have reproducible results is to create a new
[conda](https://docs.conda.io/en/latest/) environment using the
environment file provided with it.

``` bash
wget https://raw.githubusercontent.com/alexhbnr/automatic_MAG_refinement/main/environment.yml
conda env create -f environment.yml
```

After activating the environment using

``` bash
conda activate automatic_MAG_refinement
```

the necessary Python environment has been created.

Next, the pipeline itself needs to be downloaded. This can be easily
down by cloning this repository to your computer via git

``` bash
git clone https://github.com/alexhbnr/automatic_MAG_refinement.git
```

or by downloading the zip file and extracting it:

``` bash
wget -O automatic_MAG_refinement.zip https://github.com/alexhbnr/automatic_MAG_refinement/archive/refs/heads/main.zip
unzip automatic_MAG_refinement.zip
```

Finally, the configuration file and sample table have to be provided.
Templates for these can be found in `config/config.yaml` for the
configuration file and in `test/samples.tsv` for the sample table.

To run a test case using the assembly results of sample FUM003, a
Neanderthal dental calculus sample, that was first published by Fellows
Yates et al. (2021) and later *de novo* assembled in Klapper et al.
(2023):

``` bash
wget -O test/FUM003-megahit.fasta.gz https://share.eva.mpg.de/index.php/s/nQ7Df5Z4T2EFQrA
wget -O test/FUM003.sorted.dedup.bam https://share.eva.mpg.de/index.php/s/fbQNLGs74AGit6J
wget -O test/FUM003.sorted.dedup.bam.bai https://share.eva.mpg.de/index.php/s/B2nMAWLZCw5kK6y
wget -O test/metawrap_50_10_bins.stats https://share.eva.mpg.de/index.php/s/dkqeA2fNMksdqsk
```

To start the pipeline, we run

``` bash
snakemake --use-conda --conda-prefix conda \
          --use-singularity --singularity-prefix singularity -j 8
```

This will automatically evaluate the entries in the configuration file
`config/config.yaml` and use the sample `FUM003` as input. The temporary
files are written into the folder `tmp` and the results in the folder
`results`.

By activating the options `--use-conda` and `--use-singularity`,
`snakemake` will download and install the programs necessary to run the
pipeline via conda or pull the container images via singularity. This
step will only happen once, at the first time or when you change the
folder for storing these conda environments using `--conda-prefix` or
`--singularity-prefix`, respectively.

## Why this pipeline matters

The *de novo* assembly of short-read metagenomic sequencing data has
become the *de-facto* standard when studying the microbiome of complex
samples. In contrast to genomes derived from cultured microbial
isolates, the task of inferring which contig belongs to which genome is
not trivial. Therefore, multiple approaches have been developed to use
the sequence composition and the sequencing depth along the contigs to
cluster similar contigs together (Alneberg et al. 2013; Wu, Simmons, and
Singer 2016; Kang et al. 2019). While using the sequence composition and
the sequencing depth allows to cluster all contigs in a short amount of
time, this information is not always sufficient to correctly separate
all contigs in the correct genome bins. To improve these clusters,
workflows have been developed that infer the presence of
lineage-specific, single-copy genes along the contigs and use these to
revise the assignment of contigs into clusters (Sieber et al. 2018;
Uritskiy, DiRuggiero, and Taylor 2018). However, the performance of the
different combination of tools is dependent on the underlying set of
sequencing data (Yue et al. 2020).

There were attempts to standardise the evaluation of the minimum
information that is to be provided for new metagenome-assembled genomes
(MAGs) (Bowers et al. 2017). CheckM (Parks et al. 2015) has established
itself as the *de-facto* standard for estimating the completeness and
the contamination of a MAG based on marker genes. However, Orakov et al.
(2021) could show that checkM’s approach only identifies the surplus of
contigs from other taxa but is not able to identify chimeric contigs and
therefore overestimates the quality. Instead, the authors developed
their own tool, GUNC, that can evaluate whether a MAG is likely
chimeric.

How to proceed with MAGs that were either assigned to high or medium
quality using checkM but returned a higher than expected GUNC score is
yet unclear. In the presence of a large number of additional genomes
from the same habitat some researchers tended to discard these MAGs as
chimeric (e.g. Saheb Kashaf et al. (2021)). However, for samples that
are limited in quantity and underlie strong ethical considerations, such
as ancient DNA samples, this is not an adequate solution. Suggestions
have been put forward to manually curate the contigs of chimeric MAGs
(Chen et al. 2020) in programs such like anvi’o (Eren et al. 2015) and
discard the problematic contigs. This manual process ranges from
time-consuming to infeasible, when a dataset consists out of many
samples with each a large number of MAGs.

Here, we developed an automatic workflow that is heavily influenced by
the suggestions by Chen et al. (2020) and automatises many steps that
can be manually done in anvi’o. In brief, the pipeline written in
Snakemake (Mölder et al. 2021) expects MAGs refined by MetaWRAP
(Uritskiy, DiRuggiero, and Taylor 2018) as input and identifies contigs
that are likely chimeric by inferring the majority lineage across all
contigs using MMSeqs2 (Steinegger and Söding 2017) against the GTDB
reference database (Parks et al. 2020) using the command
`mmseqs taxonomy` and discard contigs that diverge either by average
sequencing depth or lineage assignment. For the revised contigs, a
standard set of assembly information including an updated estimate for
the genome completeness and the contamination is determined and
reported.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Alneberg2013" class="csl-entry">

Alneberg, Johannes, Brynjar Smári Bjarnason, Ino de Bruijn, Melanie
Schirmer, Joshua Quick, Umer Z Ijaz, Nicholas J Loman, Anders F
Andersson, and Christopher Quince. 2013. “CONCOCT: Clustering Contigs on
Coverage and Composition.” *arXiv Preprint arXiv:1312.4038*.

</div>

<div id="ref-Asnicar2020" class="csl-entry">

Asnicar, Francesco, Andrew Maltez Thomas, Francesco Beghini, Claudia
Mengoni, Serena Manara, Paolo Manghi, Qiyun Zhu, et al. 2020. “Precise
Phylogenetic Analysis of Microbial Isolates and Genomes from Metagenomes
Using PhyloPhlAn 3.0.” *Nature Communications* 11 (1): 2500.
<https://doi.org/10.1038/s41467-020-16366-7>.

</div>

<div id="ref-Bowers2017" class="csl-entry">

Bowers, Robert M, Nikos C Kyrpides, Ramunas Stepanauskas, Miranda
Harmon-Smith, Devin Doud, TBK Reddy, Frederik Schulz, et al. 2017.
“Minimum Information about a Single Amplified Genome (MISAG) and a
Metagenome-Assembled Genome (MIMAG) of Bacteria and Archaea.” *Nature
Biotechnology* 35 (8): 725–31.

</div>

<div id="ref-Chaumeil2020" class="csl-entry">

Chaumeil, Pierre-Alain, Aaron J Mussig, Philip Hugenholtz, and Donovan H
Parks. 2020. “GTDB-Tk: A Toolkit to Classify Genomes with the Genome
Taxonomy Database.” *Bioinformatics* 36 (6): 1925–27.
<https://doi.org/10.1093/bioinformatics/btz848>.

</div>

<div id="ref-Chen2020" class="csl-entry">

Chen, Lin-Xing, Karthik Anantharaman, Alon Shaiber, A Murat Eren, and
Jillian F Banfield. 2020. “Accurate and Complete Genomes from
Metagenomes.” *Genome Research* 30 (3): 315–33.

</div>

<div id="ref-Chklovski2023" class="csl-entry">

Chklovski, Alex, Donovan H. Parks, Ben J. Woodcroft, and Gene W. Tyson.
2023. “CheckM2: A Rapid, Scalable and Accurate Tool for Assessing
Microbial Genome Quality Using Machine Learning.” *Nature Methods* 20
(8): 1203–12. <https://doi.org/10.1038/s41592-023-01940-w>.

</div>

<div id="ref-Eren2015" class="csl-entry">

Eren, A Murat, Özcan C Esen, Christopher Quince, Joseph H Vineis, Hilary
G Morrison, Mitchell L Sogin, and Tom O Delmont. 2015. “Anvi’o: An
Advanced Analysis and Visualization Platform for ‘Omics Data.” *PeerJ*
3: e1319.

</div>

<div id="ref-FellowsYates2021" class="csl-entry">

Fellows Yates, James A., Irina M. Velsko, Franziska Aron, Cosimo Posth,
Courtney A. Hofman, Rita M. Austin, Cody E. Parker, et al. 2021. “The
Evolution and Changing Ecology of the African Hominid Oral Microbiome.”
*Proceedings of the National Academy of Sciences of the United States of
America* 118 (20): e2021655118.
<https://doi.org/10.1073/pnas.2021655118>.

</div>

<div id="ref-Kang2019" class="csl-entry">

Kang, Dongwan D, Feng Li, Edward Kirton, Ashleigh Thomas, Rob Egan, Hong
An, and Zhong Wang. 2019. “MetaBAT 2: An Adaptive Binning Algorithm for
Robust and Efficient Genome Reconstruction from Metagenome Assemblies.”
*PeerJ* 7: e7359.

</div>

<div id="ref-Klapper2023" class="csl-entry">

Klapper, Martin, Alexander Hübner, Anan Ibrahim, Ina Wasmuth, Maxime
Borry, Veit G. Haensch, Shuaibing Zhang, et al. 2023. “Natural Products
from Reconstructed Bacterial Genomes of the Middle and Upper
Paleolithic.” *Science* 380 (6645): 619–24.
<https://doi.org/10.1126/science.adf5300>.

</div>

<div id="ref-Mirdita2021" class="csl-entry">

Mirdita, M, M Steinegger, F Breitwieser, J Söding, and E Levy Karin.
2021. “Fast and Sensitive Taxonomic Assignment to Metagenomic Contigs.”
*Bioinformatics* 37 (18): 3029–31.
<https://doi.org/10.1093/bioinformatics/btab184>.

</div>

<div id="ref-Molder2021" class="csl-entry">

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, Michael B Hall,
Christopher H Tomkins-Tinch, Vanessa Sochat, Jan Forster, et al. 2021.
“Sustainable Data Analysis with Snakemake.” *F1000Research* 10.

</div>

<div id="ref-Orakov2021" class="csl-entry">

Orakov, Askarbek, Anthony Fullam, Luis Pedro Coelho, Supriya Khedkar,
Damian Szklarczyk, Daniel R Mende, Thomas SB Schmidt, and Peer Bork.
2021. “GUNC: Detection of Chimerism and Contamination in Prokaryotic
Genomes.” *Genome Biology* 22 (1): 1–19.

</div>

<div id="ref-Parks2020" class="csl-entry">

Parks, Donovan H, Maria Chuvochina, Pierre-Alain Chaumeil, Christian
Rinke, Aaron J Mussig, and Philip Hugenholtz. 2020. “A Complete
Domain-to-Species Taxonomy for Bacteria and Archaea.” *Nature
Biotechnology* 38 (9): 1079–86.

</div>

<div id="ref-Parks2015" class="csl-entry">

Parks, Donovan H, Michael Imelfort, Connor T Skennerton, Philip
Hugenholtz, and Gene W Tyson. 2015. “CheckM: Assessing the Quality of
Microbial Genomes Recovered from Isolates, Single Cells, and
Metagenomes.” *Genome Research* 25 (7): 1043–55.

</div>

<div id="ref-Pasolli2019" class="csl-entry">

Pasolli, Edoardo, Francesco Asnicar, Serena Manara, Moreno Zolfo,
Nicolai Karcher, Federica Armanini, Francesco Beghini, et al. 2019.
“Extensive Unexplored Human Microbiome Diversity Revealed by Over
150,000 Genomes from Metagenomes Spanning Age, Geography, and
Lifestyle.” *Cell* 176 (3): 649–662.e20.
<https://doi.org/10.1016/j.cell.2019.01.001>.

</div>

<div id="ref-Saheb2021" class="csl-entry">

Saheb Kashaf, Sara, Diana M Proctor, Clay Deming, Paul Saary, Martin
Hölzer, Monica E Taylor, Heidi H Kong, Julia A Segre, Alexandre Almeida,
and Robert D Finn. 2021. “Integrating Cultivation and Metagenomics for a
Multi-Kingdom View of Skin Microbiome Diversity and Functions.” *Nature
Microbiology*, 1–11.

</div>

<div id="ref-Schwengers2021" class="csl-entry">

Schwengers, Oliver, Lukas Jelonek, Marius Alfred Dieckmann, Sebastian
Beyvers, Jochen Blom, and AlexanderYR 2021 Goesmann. n.d. “Bakta: Rapid
and Standardized Annotation of Bacterial Genomes via Alignment-Free
Sequence Identification.” *Microbial Genomics* 7 (11): 000685. Accessed
November 22, 2022. <https://doi.org/10.1099/mgen.0.000685>.

</div>

<div id="ref-Sieber2018" class="csl-entry">

Sieber, Christian MK, Alexander J Probst, Allison Sharrar, Brian C
Thomas, Matthias Hess, Susannah G Tringe, and Jillian F Banfield. 2018.
“Recovery of Genomes from Metagenomes via a Dereplication, Aggregation
and Scoring Strategy.” *Nature Microbiology* 3 (7): 836–43.

</div>

<div id="ref-Steinegger2017" class="csl-entry">

Steinegger, Martin, and Johannes Söding. 2017. “MMseqs2 Enables
Sensitive Protein Sequence Searching for the Analysis of Massive Data
Sets.” *Nature Biotechnology* 35 (11): 1026–28.

</div>

<div id="ref-Uritskiy2018" class="csl-entry">

Uritskiy, Gherman V, Jocelyne DiRuggiero, and James Taylor. 2018.
“MetaWRAP—a Flexible Pipeline for Genome-Resolved Metagenomic Data
Analysis.” *Microbiome* 6 (1): 1–13.

</div>

<div id="ref-Wu2016" class="csl-entry">

Wu, Yu-Wei, Blake A Simmons, and Steven W Singer. 2016. “MaxBin 2.0: An
Automated Binning Algorithm to Recover Genomes from Multiple Metagenomic
Datasets.” *Bioinformatics* 32 (4): 605–7.

</div>

<div id="ref-Yue2020" class="csl-entry">

Yue, Yi, Hao Huang, Zhao Qi, Hui-Min Dou, Xin-Yi Liu, Tian-Fei Han, Yue
Chen, Xiang-Jun Song, You-Hua Zhang, and Jian Tu. 2020. “Evaluating
Metagenomics Tools for Genome Binning with Real Metagenomic Datasets and
CAMI Datasets.” *BMC Bioinformatics* 21 (1): 1–15.

</div>

</div>
