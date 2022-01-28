Automatic revision of MAGs to remove likely chimeric contigs
================

  - [Background](#background)
  - [References](#references)

# Background

The *de-novo* assembly of short-read metagenomic sequencing data has
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

In the last five years, there have been made attempts to standardise the
minimum information that is to be provided for new metagenome-assembled
genomes (MAGs) (Bowers et al. 2017). Two of the most important criteria
here are the genome completeness and its contamination, which are hard
to estimate in the absence of a ground truth. CheckM (Parks et al. 2015)
has evolved to be the quasi standard tool to estimate these two
parameters by inferring single-copy marker genes on the contigs of each
MAG and evaluating these against an expected set of genes that is
determined by assigning the MAG to a microbial lineage. Based on the
suggestions by Bowers et al. (2017), a MAG is considered to be of
high-quality when the completeness estimate \>= 90% and the
contamination estimate \< 5%. A medium-quality MAG is required to have a
completeness \>= 50% and a contamination \< 10%.

However, in a recent publication, Orakov et al. (2021) could show that
while checkM is able to identify contigs that do not belong to the
consensus microbial lineage of the MAG, checkM fails to identify the
presence of chimeric contigs in MAGs and thus overestimates the quality
of MAG. To estimate the presence of chimeric contigs in a data set,
Orakov et al. (2021) uses a very similar approach to checkM. The authors
developed GUNC, which uses a database with genes and their known
taxonomic origin, and identifies their presence on the contigs. If the
authors observe a higher variety of taxonomic origins than expected from
known microbial species, they will flag the MAG as likely be chimeric.

How to proceed with MAGs that were either assigned to high or medium
quality using checkM but returned a higher than expected GUNC score is
yet unclear. In the presence of a large number of additional genomes
from the same habitat some researchers tended to discard these MAGs as
chimeric (e.g. Saheb Kashaf et al. (2021)). However, for samples that
are limited in quantity and underlie strong ethical considerations, such
as ancient DNA samples, this is not an adequate solution. Suggestions
have been put forward to manually curate the contigs of chimeric MAGs
(Chen et al. 2020) in programs such like anvi’o (Eren et al. 2015) and
discard the problematic contigs. This manual process ranges from
time-consuming to infeasible, when a dataset consists out of many
samples with each a large number of MAGs.

In the following, I present an automatic workflow that is heavily
influenced by the suggestions by Chen et al. (2020) and automatises many
steps that can be manually done in anvi’o. In brief, the pipeline
written in Snakemake (Mölder et al. 2021) expects MAGs refined by
MetaWRAP (Uritskiy, DiRuggiero, and Taylor 2018) as input and identifies
contigs that are likely chimeric by inferring the majority lineage
across all contigs using MMSeqs2
(<span class="citeproc-not-found" data-reference-id="Steinegger2018">**???**</span>)
against the GTDB reference database (Parks et al. 2020) using the
command `mmseqs taxonomy` and discard contigs that diverge either by
average sequencing depth or lineage assignment. For the revised contigs,
a standard set of assembly information including an updated estimate for
the genome completeness and the contamination is determined and
reported.

# References

<div id="refs" class="references">

<div id="ref-Alneberg2013">

Alneberg, Johannes, Brynjar Smári Bjarnason, Ino de Bruijn, Melanie
Schirmer, Joshua Quick, Umer Z Ijaz, Nicholas J Loman, Anders F
Andersson, and Christopher Quince. 2013. “CONCOCT: Clustering Contigs on
Coverage and Composition.” *arXiv Preprint arXiv:1312.4038*.

</div>

<div id="ref-Bowers2017">

Bowers, Robert M, Nikos C Kyrpides, Ramunas Stepanauskas, Miranda
Harmon-Smith, Devin Doud, TBK Reddy, Frederik Schulz, et al. 2017.
“Minimum Information About a Single Amplified Genome (Misag) and a
Metagenome-Assembled Genome (Mimag) of Bacteria and Archaea.” *Nature
Biotechnology* 35 (8): 725–31.

</div>

<div id="ref-Chen2020">

Chen, Lin-Xing, Karthik Anantharaman, Alon Shaiber, A Murat Eren, and
Jillian F Banfield. 2020. “Accurate and Complete Genomes from
Metagenomes.” *Genome Research* 30 (3): 315–33.

</div>

<div id="ref-Eren2015">

Eren, A Murat, Özcan C Esen, Christopher Quince, Joseph H Vineis, Hilary
G Morrison, Mitchell L Sogin, and Tom O Delmont. 2015. “Anvi’o: An
Advanced Analysis and Visualization Platform for ‘Omics Data.” *PeerJ*
3: e1319.

</div>

<div id="ref-Kang2019">

Kang, Dongwan D, Feng Li, Edward Kirton, Ashleigh Thomas, Rob Egan, Hong
An, and Zhong Wang. 2019. “MetaBAT 2: An Adaptive Binning Algorithm for
Robust and Efficient Genome Reconstruction from Metagenome Assemblies.”
*PeerJ* 7: e7359.

</div>

<div id="ref-Molder2021">

Mölder, Felix, Kim Philipp Jablonski, Brice Letcher, Michael B Hall,
Christopher H Tomkins-Tinch, Vanessa Sochat, Jan Forster, et al. 2021.
“Sustainable Data Analysis with Snakemake.” *F1000Research* 10.

</div>

<div id="ref-Orakov2021">

Orakov, Askarbek, Anthony Fullam, Luis Pedro Coelho, Supriya Khedkar,
Damian Szklarczyk, Daniel R Mende, Thomas SB Schmidt, and Peer Bork.
2021. “GUNC: Detection of Chimerism and Contamination in Prokaryotic
Genomes.” *Genome Biology* 22 (1): 1–19.

</div>

<div id="ref-Parks2020">

Parks, Donovan H, Maria Chuvochina, Pierre-Alain Chaumeil, Christian
Rinke, Aaron J Mussig, and Philip Hugenholtz. 2020. “A Complete
Domain-to-Species Taxonomy for Bacteria and Archaea.” *Nature
Biotechnology* 38 (9): 1079–86.

</div>

<div id="ref-Parks2015">

Parks, Donovan H, Michael Imelfort, Connor T Skennerton, Philip
Hugenholtz, and Gene W Tyson. 2015. “CheckM: Assessing the Quality of
Microbial Genomes Recovered from Isolates, Single Cells, and
Metagenomes.” *Genome Research* 25 (7): 1043–55.

</div>

<div id="ref-Saheb2021">

Saheb Kashaf, Sara, Diana M Proctor, Clay Deming, Paul Saary, Martin
Hölzer, Monica E Taylor, Heidi H Kong, Julia A Segre, Alexandre
Almeida, and Robert D Finn. 2021. “Integrating Cultivation and
Metagenomics for a Multi-Kingdom View of Skin Microbiome Diversity and
Functions.” *Nature Microbiology*, 1–11.

</div>

<div id="ref-Sieber2018">

Sieber, Christian MK, Alexander J Probst, Allison Sharrar, Brian C
Thomas, Matthias Hess, Susannah G Tringe, and Jillian F Banfield. 2018.
“Recovery of Genomes from Metagenomes via a Dereplication, Aggregation
and Scoring Strategy.” *Nature Microbiology* 3 (7): 836–43.

</div>

<div id="ref-Uritskiy2018">

Uritskiy, Gherman V, Jocelyne DiRuggiero, and James Taylor. 2018.
“MetaWRAP—a Flexible Pipeline for Genome-Resolved Metagenomic Data
Analysis.” *Microbiome* 6 (1): 1–13.

</div>

<div id="ref-Wu2016">

Wu, Yu-Wei, Blake A Simmons, and Steven W Singer. 2016. “MaxBin 2.0: An
Automated Binning Algorithm to Recover Genomes from Multiple Metagenomic
Datasets.” *Bioinformatics* 32 (4): 605–7.

</div>

<div id="ref-Yue2020">

Yue, Yi, Hao Huang, Zhao Qi, Hui-Min Dou, Xin-Yi Liu, Tian-Fei Han, Yue
Chen, Xiang-Jun Song, You-Hua Zhang, and Jian Tu. 2020. “Evaluating
Metagenomics Tools for Genome Binning with Real Metagenomic Datasets and
Cami Datasets.” *BMC Bioinformatics* 21 (1): 1–15.

</div>

</div>
