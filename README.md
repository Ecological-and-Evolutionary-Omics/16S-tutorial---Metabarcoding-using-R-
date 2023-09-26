# 16S tutorial - METABARCODING (Using R)
>[!WARNING]
>This tutorial will work with R code. Usage of R Studio is recommended for a better understading of the commands used

## Table of Contents

- [Download the tutorial data](#download the tutorial data)
- [Installing Packages](#installing Packages)
- [Quality Control](#quality Control)


## Download the tutorial data

:heavy_exclamation_mark: The next data can be downloaded manually if you are not using shell by searching the links shown in the next commands

```Bash
wget http://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip
wget https://zenodo.org/record/3731176/files/silva_nr_v138_train_set.fa.gz?download=1
wget https://zenodo.org/record/3731176/files/silva_species_assignment_v138.fa.gz?download=1

#unzip the compressed document MISEQSOPData.zip and remove the __MACOSX/ directory
#it is also recommended to move the SILVA files into the MiSeq_S0P directory
unzip MiSeqSOPData.zip
rm -r __MACOSX/
```
## Installing Packages
For installing all the packages needed for the 16S analysis first we need to install Bioconductor :inbox_tray:


```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.17")
```

Then install the next packages 
```R
BiocManager::install('dada2')
BiocManager::install('phyloseq')
BiocManager::install('DECIPHER')
install.packages('ggplot2')
install.packages('phangorn')
```

## Quality Control

Quality Control of the samples can be performed using R or using the Shell scripts explained in the [Quality Control turorial](https://github.com/Ecological-and-Evolutionary-Omics/Quality-and-Trimming)

```R
#if you are using Windows remember to change the "\" in the path for "/"
path <- 'MiSeq_SOP'
list.files(path)
```
We are going to automatize a huge part of the process by making a list of all the files containing _R1_001.fastq or _R2_001.fastq
```R
raw_forward <- sort(list.files(path, pattern="_R1_001.fastq.gz",
                               full.names=TRUE))

raw_reverse <- sort(list.files(path, pattern="_R2_001.fastq.gz",
                               full.names=TRUE))

# we also need the sample names
sample_names <- sapply(strsplit(basename(raw_forward), "_"),
                       `[`,  # extracts the first element of a subset
                       1)
# Another option will be
sample_names <- gsub("_R1_001.fastq.gz", "", raw_forward)
```
### Check Quality of the reads

```R
plotQualityProfile(raw_forward[1:2]) 
plotQualityProfile(raw_reverse[1:2])

#you can change the number of samples shown by changing the [1:2] parameter
#to any type of selection 
```
Now an image like this should be visible in your screen (!This is just an example picture don't panic :scream:, your image should look different :thumbsup:)
![IMAGE](C:/Users/1001926/Desktop/R_quality.png)
