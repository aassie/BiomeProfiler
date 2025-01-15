# BiomeProfiler #

This is the temporary Readme file of the Biomeprofiler tool. This script is still in active development and more will come soon!

# Installation

## Pre requirements

To use Biomeprofiler you will need to have:

- Local installation of BLAST
- [R](https://cran.r-project.org/), with the following packages:
  - [tidyverse](https://www.tidyverse.org/)
  - [vegan](https://cran.r-project.org/web/packages/vegan/)
  - [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer)
  - [labdsv](https://cran.r-project.org/web/packages/labdsv/index.html)
  - [FSA](https://cran.r-project.org/web/packages/FSA/index.html)
  - Bioconductor(https://bioconductor.org/)
    - [Phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html)
    - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
  - [composition](https://cran.r-project.org/web/packages/compositions/index.html)
- Our consensus annotation file. It is a big file, you can download it from Zenodo [here](https://zenodo.org/records/10689575). Once you have it place the file in the `Data/` folder and rename it `WB.tsv`.

## Conda

There is currently no conda package available for Biomeprofiler. But we advise you to use conda to manually create a environment to run the pipeline

###  Create the environment

```
conda create --name biomeprofiler
conda activate biomeprofiler
conda config --add channels defaults && conda config --add channels bioconda && conda config --add channels conda-forge
```

### Install the tools

```
# basic dependencies
conda install bash r perl parallel gawk sed grep bc git coreutils wget

# install blast
conda install blast

# R-package dependencies (via conda repos)
conda install r-tidyverse r-readr r-reshape2 r-vegan r-RcolorBrewer r-labdsv r-FSA r-biocmanager bioconductor-phyloseq bioconductor-deseq2 r-composition
```
Once you have those package installed, clone the git repository

```
git clone https://aassiebcm@bitbucket.org/the-samuel-lab/biomeprofiler.git && cd biomeprofiler
```

# Run the Pipeline

## The files needed to run the pipeline

### Bacteria Collection

To run the pipeline you need to provide a file with the 16S rRNA sequences of the bacteria you used in your experiments. Once you have this fasta file, place it in the `Data/Collection` folder.

Currently there are a few example in the `Data/Collection` folder, mainly the CeMbio and Big68 communities.

### ASV table

This table is your standard output of a 16S rRNA pipeline, such as [Qiime2](https://qiime2.org/), Mothur(https://mothur.org/) or directly from [Deblur](https://github.com/biocore/deblur) or [DaDa2](https://benjjneb.github.io/dada2/)

At the moment the Sample ID column needs to be named `X.SampleID`. I am working on making it smarter.

### Representative ASV file

This is a fasta file with the representative sequences for the ASV name used in the ASV table above. It should be an output from your 16S rRNA pipeline.

### Metadata

At the moment the Sample ID column needs to be named `X.SampleID`. I am working on making it smart. Other than that it should be a tab delimited file with all the information you think is needed. To run the statistic module from the pipeline you will give the name of one of the column in your metadata file.

## Running the script

This part is still in active development, I hope to have all the core information updated here very soon!

```
~/Your/Path/biomeprofiler.sh -a ASV.table.tsv -m metadata.txt -s ASV.sequences.fasta -b BacteriaCollectionName -o Output_Folder -s1 TestGroup -v -st
```

`BacteriaCollectionName` correspond to the name of the Bacteria Collection file you create or use (see above)

## In details

Below is detailed explanation of what each option stands for

```
Options:
   -h|--help                Display this help

   -a|--ASVtable            Required. Biom/tsv/csv File

   -m|--ASVmeta             Required. Count table related Metadata

   -s|--ASVsequence         Required. related ASV sequences

   -b|--BacteriaList        Required. Collection Name or list of bacteria
                            name used in the experiment. One per line

   -o|--OutFolder           Output OutFolder
                            Default: ./OUT/

   -s1|--Selector1          Name of column of the metadata category that defines
                            your dataset. E.g: "Host Strain"

   -k|--kegg                Generate Kegg annotation table

   -p|--PathwayConvert      Generate Metacyc Metabolic Pathway Table.

   -pl|--PathwayTable       If Metabolic Pathway is already generated,
                            indicate the location of the file here
                            Incompatible with -p

   -kl|--keggTable          If Kegg annotation Table is already generated,
                            indicate the location of the file here
                            Incompatible with -k

   -bn|--BacteriaName       Bacteria Collection Name

   -st|--Stats              Run Deseq2 using the "selector1" as category to test

   -v|--Verbose             Verbose mode, the script will talk to you, a lot
```
