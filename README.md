# Variant-Detection-and-Annotation-workshop
The repo is intended to accompany the NYUAD Core Bioinformatics hands-on workshop on Variant Detection and Annotation.

During this 2-day workshop, participants will learn about the steps involved in calling variants (SNPs and Indels) from short read high throughput sequencing data.

Participants are expected to have some basic knowledge of mutation analysis as well as some biological knowledge on the subject matter.

Although the material and the methods are designed to cater for the NYU Abu Dhabi High Performance Computing environment, it is possible to run this workshop on any other system provided that the neccessary software is installed, and the data is uploaded to that environment. It will also be neccessary to change the input and output directories/files to accomodate such an environment.

The workshop will take the format of case studies with the aim of discovering the underlying mutations in 3 separate patients that are exhibiting some disease phenotypes. Our starting point will be raw sequencing files (Illumina paired end short reads), and throughout the course of the workshop, we will learn how to process and analyze the data so that we end up with annotated VCF (Variant Calling Format) files.

More specifically, we will:

  1. Perform Quality Checking and Quality Trimming (QC/QT) on the raw data.
  2. Align the data to the reference human genome (version HG38).
  3. Carry out the necessary alignment post processing steps following established best practice guides.
  4. Call variants (SNPs and Indels) and filter the variants that have been called.
  5. Annotate the Variants and attempt to establish causative mutations.

By the end of this workshop, participants should be able to replicate these analysis steps on any DNA sequencing dataset (WGS/WES/Panels) originating from short read sequencing technologies.


Enjoy!



## Setting up the environment and copying the data
We will be using the NYUAD High Performance Computing (HPC) cluster for this workshop, however, you can certainly run all of the analysis on any stand alone machine (server, personal laptop/Desktop etc.) provided that you have pre-installed the necessay software packages.

Before starting the workshop, you would've been assigned a dataset, either p1, p2, or p3. Once you have completed your analysis on one of these datasets, you can certainy run the steps again for the other datasets.

### Connecting to the HPC using a MAC/Linux machine and copying the data.
1. Open the "Terminal" app and type `ssh NetID@jubail.abudhabi.nyu.edu`. Enter your NYU password when prompted and hit enter.
2. Once logged in, navigate to your personal "SCRATCH" directory `cd $SCRATCH`.
3. Create a directory and change into it `mkdir variant_detection && cd variant_detection`.
4. Copy the dataset that has been assigned to you (either p1, p2, or p3) e.g. for p1 `cp -r p1 .`.

### Connecting to the HPC using a Windows machine and copying the data.
1. Open the "Putty" app, and fill out the fields as follows **Host name**=jubail.abudhabi.nyu.edu, **Port**=22, and then click on "Open".
2. Enter your NetId, and your password when prompted.
3. Once logged in, navigate to your personal "SCRATCH" directory `cd $SCRATCH`.
4. Create a directory and change into it `mkdir variant_detection && cd variant_detection`.
5. Copy the dataset that has been assigned to you (either p1, p2, or p3) e.g. for p1 `cp -r p1 .`.


## The software stack
As you might've imagined, this sort of analysis involves multiple steps, and multiple tools. The tools and their links are provided below in case you want to run this analysis on your own setup. Since we are using the HPC, we don't have to install any of them since they are already there!

Just a quick note on installation. Whenever possible, we recommend using conda for installing and maintaining your software stack. Installing Bioinformatics software from source can be a painful experience sometimes, and conda takes care of most cases with relative ease.

- FastQC [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/]
- FastP [https://github.com/OpenGene/fastp]
- BWA [https://bio-bwa.sourceforge.net/]
- SAMTools [https://www.htslib.org/]
- PICARD tools [https://broadinstitute.github.io/picard/]
- GATK [https://gatk.broadinstitute.org/hc/en-us]
- Qualimap [http://qualimap.conesalab.org/]
- SNPEff and SnpSift [https://pcingola.github.io/SnpEff/]
- ANNOVAR [https://annovar.openbioinformatics.org/en/latest/]

## File Formats
We will be introducing a number of different file formats that are common in many types of omics analyses. Below is a list of these formats and links describing what they are, and although we will be covering all of these in the workshop and the accompanying slide deck, we thought it's a good idea for you to have a "Quick links" section to them.

- FASTA [https://en.wikipedia.org/wiki/FASTA_format]
- FASTQ [https://en.wikipedia.org/wiki/FASTQ_format]
- SAM [https://en.wikipedia.org/wiki/SAM_(file_format)] 
- BAM [https://en.wikipedia.org/wiki/Binary_Alignment_Map]
- VCF/GVCF [https://en.wikipedia.org/wiki/Variant_Call_Format]



## The data
The datasets that we will be using for this workshop are unfiltered FASTQ sequencing files. The are publicly available to download from the Short Read Archive (SRA [https://www.ncbi.nlm.nih.gov/sra]) using the following accessions ERR674822, SRR24555542, SRR5604284. Our p1, p2, and p3 datasets are reduced from their original size, meaning that not all of the sequencing reads are available. The reason for this is due to time constraints. Running the complete dataset takes much longer than the time we have during this workshop.

As mentioned earlier, each one of you has been preallocated a particular dataset. Your task is to complete the analysis of this data and find out which one of the "case files" relates to the dataset that you have been given.

### Case file A
This case file exhibits the following:
- History of illness in the family but both parents are healthy.
- Asthma like symptoms, and recurring lung infections.
- Chronic constipation.

### Case file B
This case file exhibits the following:
- No history of illness in the family.
- Poor ability to perform excercises as well as poor coordination.
- Signs of poor spinal development.
- Blood tests indicated that Serum Creatinine Kinase levels were 12x the normal values.

### Case file C
Unfortunately, all of the information of this case file has been lost! It is therefore up to you to find out the disease associated with your dataset is.

By analyzing your given dataset, you will be able to call mutations. Through careful examination and filtering of the data, you should be able to successfully complete these tasks.



