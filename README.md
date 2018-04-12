# BiocNYC-ChIPSeqSpike

# Performing ChIP-Seq data scaling with spike-in control in R/Bioconductor

04/19/2018

## Introduction

In this workshop, we will look at performing ChIP-Seq data scaling with spike-in controls. We will use 'ChIPSeqSpike' that is a new package at an early stage of development. You will learn how to perform the spike-in scaling properly and to visualize the different transformations in R. You will also be able to compare the effect of an inhibitor treatment on the Histone 3 Lysine 79 di-methylation (H3K79me2).

This workshop will not cover the data pre-processing and already processed BAM and BIGWIG files will be provided.


### Covering today

* Background on ChIP-Seq spike-in (see [pdf](https://github.com/descostesn/BiocNYC-ChIPSeqSpike/blob/master/BioconductorMeetupAp2018.pdf))

* Construction of objects containing all files necessary for proper scaling

* Computation and evaluation of scaling factors 

* Data Visualization


### The data

In this workshop we will make use of one dataset of published data. This dataset is coming from the article:

Orlando, David A, Mei Wei Chen, Victoria E Brown, Snehakumari Solanki, Yoon J Choi, EricR Olson, Christian C Fritz, James E Bradner, and Matthew G Guenther. 2014. “**Quantitative ChIP-Seq Normalization Reveals Global Modulation of the Epigenome.**” Cell Reports 9 (3). Cell press:1163–70.


It represents a gold-standard example of the importance of using spike-in controls with ChIP-Seq experiments. It uses chromatin from Drosophila Melanogaster as exogenous spike-in control to correct experimental biases. Without a spike-in control and using only RPM normalization, proper differences in the H3K79me2 histone modification in human Jurkat cells upon EPZ5676 inhibitor treatment were not observed.
This dataset is made of bigwig and bam files of H3K79me2 ChIP-Seq data and corresponding input DNA controls. Bam files contain data aligned to the Human reference genome Hg19 or to the Drosophila reference genome dm3. The latest is used to compute external spike-in scaling factors. All above mentioned data are available at 0%, 50% and 100% EPZ5676 inhibitor treatment.

The whole dataset is accessible at [GSE60104](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60104). Specifically, the data used are H3K79me2 0% (GSM1465004), H3K79me2 50% (GSM1465006), H3K79me2 100% (GSM1465008), input DNA 0% (GSM1511465), input DNA 50% (GSM1511467) and input DNA 100% (GSM1511469).

You can download the material for this workshop below:

* H3K79me2_0: [hg19.bw](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/H3K79me2_0.bw), [dm3.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/H3K79me2_0_dm3.bam), [hg19.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/H3K79me2_0_hg19.bam)

* H3K79me2_50: [hg19.bw](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/H3K79me2_50.bw), [dm3.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/H3K79me2_50_dm3.bam), [hg19.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/H3K79me2_50_hg19.bam)

* H3K79me2_100: [hg19.bw](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/H3K79me2_100.bw), [dm3.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/H3K79me2_100_dm3.bam), [hg19.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/H3K79me2_100_hg19.bam)

* input_0: [hg19.bw](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/input_0.bw), [dm3.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/input_0_dm3.bam), [hg19.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/input_0_hg19.bam)

* input_50: [hg19.bw](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/input_50.bw), [dm3.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/input_50_dm3.bam), [hg19.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/input_50_hg19.bam)

* input_100: [hg19.bw](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/input_100.bw), [dm3.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/input_100_dm3.bam), [hg19.bam](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/input_100_hg19.bam)

* Info table: [info.csv](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/info.csv)

* Gene annotations: [refseq_hg19.gff](http://www.hpc.med.nyu.edu/~descon01/biocnycworkshop/refseq_hg19.gff)


Create a folder 'workshop_files' containing all the downloaded files.


### Required R packages and versions.

We will need the devel version of R v(3.5) and Bioconductor (v3.7) with the following packages. To install R-devel you can follow this [tutorial](http://singmann.org/installing-r-devel-on-linux/).

```{r packages, echo=TRUE,eval=FALSE,cache=FALSE}
install.packages("tools", dependencies=TRUE)
install.packages("stringr", dependencies=TRUE)
install.packages("ggplot2", dependencies=TRUE)
install.packages("LSD", dependencies=TRUE)
install.packages("corrplot", dependencies=TRUE)
install.packages("methods", dependencies=TRUE)
install.packages("stats", dependencies=TRUE)
install.packages("grDevices", dependencies=TRUE)
install.packages("graphics", dependencies=TRUE)
install.packages("utils", dependencies=TRUE)

source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer")
biocLite("Rsamtools")
biocLite("GenomicRanges")
biocLite("IRanges")
biocLite("seqplots")
biocLite("BiocGenerics")
biocLite("S4Vectors")```

After downloading ChIPSeqSpike [here](https://bioconductor.org/packages/3.7/bioc/html/ChIPSeqSpike.html), run the following command in a terminal:

```console
sudo bash R-devel CMD INSTALL ChIPSeqSpike_0.99.22.tar.gz
```


## S4 object generation, summary and controls

First define the different pathes and folders needed:

```{r packages, echo=TRUE,eval=FALSE,cache=FALSE}
library("ChIPSeqSpike")

info_file_csv <- "mypath/workshop_files/info.csv"
bam_path <- "mypath/workshop_files"
bigwig_path <- "mypath/workshop_files"
gff_vec <- "mypath/workshop_files/refseq_hg19.gff"
genome_name <- "hg19";
col_vec <- c("red", "blue", "green");
output_folder <- "mypath/workshop_files"

#Observe the structure of the info.csv file
info_file <- read.csv(info_file_csv)
head(info_file)
## expName    endogenousBam   exogenousBam
## 1 H3K79me2_0		H3K79me2_0_hg19-filtered.bam		H3K79me2_0_dm3-filtered.bam
## 2 H3K79me2_50	H3K79me2_50_hg19-filtered.bam		H3K79me2_50_dm3-filtered.bam
## 3 H3K79me2_100 	H3K79me2_100_hg19-filtered.bam 	H3K79me2_100_dm3-filtered.bam
## inputBam   bigWigEndogenous    bigWigInput
## 1 input_0_hg19-filtered.bam 		H3K79me2_0-filtered.bw 		input_0-filtered.bw
## 2 input_50_hg19-filtered.bam		H3K79me2_50-filtered.bw		input_50-filtered.bw
## 3 input_100_hg19-filtered.bam	H3K79me2_100-filtered.bw		input_100-filtered.bw
```
The different data necessary for proper spike-in scaling are provided in a csv or a tab separated txt file. The columns must contain proper names and are organized as follows: Experiment name (expName); bam file name of data aligned to the endogenous reference genome (endogenousBam); bam file name of data aligned to the exogenous reference genome (exogenousBam); the corresponding input DNA bam file aligned to the endogenous reference genome (inputBam); the fixed steps bigwig file name of data aligned to the endogenous reference genome (bigWigEndogenous) and the fixed steps bigwig file names of the corresponding input DNA experiment aligned to the endogenous reference genome (bigWigInput).

From the info file, two kinds of objects can be generated: either a ChIPSeqSpikeDataset or a ChIPSeqSpikeDatasetList depending upon the number of input DNA experiments. A ChIPSeqSpikeDatasetList object is a list of ChIPSeqSpikeDataset object that is created if several input DNA experiments are used. In this latter case, ChIP-Seq experiments are grouped by their corresponding input DNA. The function spikeDataset creates automatically the suitable object. The folder path to the bam and fixed steps bigwig files must be provided.
If one have access to high performance computing facilities, ChIPSeqSpike offers a boost mode. This mode stores binding scores for each experiment in a GRanges object. Refer to the package vignette for more details as this mode will not be used today.

```{r packages, echo=TRUE,eval=FALSE,cache=FALSE}
csds_test <- spikeDataset(info_file_csv, bam_path, bigwig_path)
is(csds_test)
## [1] "ChIPSeqSpikeDatasetList"
```

