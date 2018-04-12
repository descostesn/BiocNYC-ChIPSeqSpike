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
