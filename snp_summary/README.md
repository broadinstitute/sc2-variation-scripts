# Introduction

The contents of this directory relate to summarizing SARS-CoV-2 SNP frequencies in the global population of COVID-19 cases sampled by genomic sequencing.

The sequence and codon numbering are relative to the coordinate space of `hCoV-19/Wuhan/WIV04/2019` ([NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)). The spike protein with D614G spans 21563—25384bp.

# Setup

 - Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
 - [Set up conda channels to source from bioconda](https://bioconda.github.io/user/install.html#set-up-channels)
 - Install analysis dependencies by creating a [conda](https://docs.conda.io/en/latest/miniconda.html) environment based on the environemt file in this repository, `sc2-analysis.yml`.
```conda env create --file sc2-analysis.yml```

# Input data

These scripts are intended to analyze the Nextstrain-curated data from [GISAID](https://www.gisaid.org/), which requires registration to download data. GISAID prohibits redistribution, so this repository cannot include input data. 

To obtain the input data, follow these steps:
 1. Log into GISAID’s [EpiCoV site](https://www.epicov.org/epi3/frontend)
 2. Click “Downloads” to bring up a modal window
 3. In this window click on “nextmeta” to download the file `metadata.tsv.bz2`. This should be decompressed and saved as `data/metadata.tsv`.
 4. From the same modal window, click on “nextfasta” to download the file `sequences.fasta.bz2`. This should be decompressed and saved as `data/sequences.fasta`.
 5. Download a fasta-format file for the SARS-CoV-2 reference sequence, [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2), and save as `data/NC_045512.2.fasta`.

The sequences were prepared via:
```../scripts/align_and_clean.sh data/sequences.fasta data/blank.fasta data/NC_045512.2.fasta data/cleaned.fasta```

# Summarizing SNPs

```./snp_summary.sh```