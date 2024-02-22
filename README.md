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
4. Copy the dataset that has been assigned to you (either p1, p2, or p3) e.g. for p1 `cp -r /scratch/gencore/nd48/workshop_variant_detection/data/p1 .`.

### Connecting to the HPC using a Windows machine and copying the data.
1. Open the "Putty" app, and fill out the fields as follows **Host name**=jubail.abudhabi.nyu.edu, **Port**=22, and then click on "Open".
2. Enter your NetId, and your password when prompted.
3. Once logged in, navigate to your personal "SCRATCH" directory `cd $SCRATCH`.
4. Create a directory and change into it `mkdir variant_detection && cd variant_detection`.
5. Copy the dataset that has been assigned to you (either p1, p2, or p3) e.g. for p1 `cp -r /scratch/gencore/nd48/workshop_variant_detection/data/p1 .`.


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
- IGV [https://www.igv.org/]

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


## Step 1: Data quality checking and quality trimming
The very first step in almost all genomics analyses involves checking the quality of the sequencing data and applying some initial quality filtering techniques (more on this in the slide deck).
Sequencing data is not bullet-proof. There are quite a few issues that can affect your data including, bad calling and low quality score calls, contamination, duplication (optical and PCR), low yield, off-targets (for capture kits) etc. 
In order to assess the quality of our data, as well filter our data to include only high quality sequencing reads, we will be using two popular tools, FastQC(for quality checking) and FastP(for quality filtering).

The process is as follows,
1. Run FastQC on the unfiltered data.
2. Run the data through FastP to apply quality filtering.
3. Run the data through FastQC again to determine if the filtering worked and how much data we lost.

In the p1, p2, or p3 folder that you copied earlier, you will find a script called **variant_detection.sh**. This is meant to either be submitted in the HPC SLURM queue, or for you to run directly on your own setup (don't worry, we will show you how).

You will need a text editor to edit this script as we go along. Some popular ones include Brackets [https://brackets.io/] and Atom [https://atom-editor.cc/], or if you prefer to use the command line, vi is a popular one [https://www.redhat.com/sysadmin/introduction-vi-editor].

For the purpose of this workshop, and so that we are all on the same page, we will be using the HPC interactive dashboard.
1. Open the following link (you need an NYU HPC account for this, otherwise just use your own text editor) [https://ood.hpc.abudhabi.nyu.edu/pun/sys/dashboard].
2. Log in to your NYU account when prompted.
3. Click on **Files** and select **/scratch/NetID/** from the drop down menu.
4. Click on the **variant_detection** folder, and then on your dataset folder (p1, p2, or p3).
5. To the right of the **variant_detection.sh**, click on the menu and select "edit" from the drop down list.

The first few lines of the script are SLURM specific lines, and they are intended to be interpreted by the job scheduling program on the HPC (They start with #SBATCH). These lines instruct the HPC on how many CPUs we need, how much Memory and for how long we are intending to run this script. If you are not running this on the HPC, then these would simply be ignored when you execute your program.

Following these lines, we have a number of module loading commands. These commands "load" the required software packages so that they will be available to us during our analysis. Again, if you are running this script on your own (not the HPC at NYUAD), then you will simply need to add Hashtags infront of them (e.g. #module load all instead of module load).

After these initial "setting up" lines, we come to the main part of our analysis, and that is running FastQC, currently they look like this,
```
#fastqc \
#    --extract p2_r1.fastq \
#    -o raw_qc/read1qc/ \
#    -t 28
```
So you need to delete the hashtags.

The command above will run fastqc on the read1s, we also run the same command for read2s. Remember, these are illumina paired end sequencing reads.

As you can see, we have passed additional arguments to these commands by using "flags", e.g. "-t 28" which instructs FastQC to ustilize 28 CPU threads (for faster processing). This number is not arbitrary and it needs to match your current setup, so if your machine has 8 CPUs, this will need to be adjusted accordingly. There are many flags for various software, and to find out what they mean, or to add additional flags, just type `fastqc --help`.

Now go ahead and remove the remaining hashtags from the read2 fastqc command, the fastp command, and the fastqc commands following fastp. So your script should look like this now,
```
# Run FastQC on the unfiltered read 1
fastqc \
    --extract p2_r1.fastq \
    -o raw_qc/read1qc/ \
    -t 28

# Run FastQC on the unfiltered read 2
fastqc \
    --extract p2_r2.fastq \
    -o raw_qc/read2qc/ \
    -t 28


# Perform quality trimming using FastP using the default parameters and using the adapter trim Fasta file
fastp \
    -i p2_r1.fastq \
    -I p2_r2.fastq \
    -o fastp_p2_r1.fastq \
    -O fastp_p2_r2.fastq \
    --adapter_fasta adapter_trim.fa \
    -h p2.fastp.html \
    -R P2.sample_QT \
    -w 16


# Create the directories for the quality filtered FASTQ reports (produced using FastQC)
mkdir -p qt_qc/read1qc
mkdir qt_qc/read2qc


# Run FastQC for the quality filtered read 1
fastqc \
    --extract fastp_p2_r1.fastq \
    -o qt_qc/read1qc/ \
    -t 28

# Run FastQC for the quality filtere read 2
fastqc \
    --extract \
    -o qt_qc/read2qc/ fastp_p2_r2.fastq \
    -t 28
```

Before running the script, let's take a few minutes to decipher what the flags that we have passed to the quality trimming tool FastP mean. To do that, run `fastp --help` on the command line (make sure that you have loaded or installed fastp first otherwise you will get a "command not found" error).
```
module purge
module load all
module load gencore/2
module load fastp/0.20.1
fastp --help
```

Great! So now that we understand what the parameters mean, and how these commands are linked together (the output of the fastp quality trimmed reads are then passed to fastqc for checking), let's go ahead and run them.

On the HPC `sbatch variant_detection.sh`

On a standalone machine (not the HPC e.g. laptop) 
```
chmod 755 variant_detection.sh
./variant_detection.sh
```
Once these steps complete, you should have two folder, **qt_qc** and **raw_qc**, as well as the quality trimmed fastq files for read 1 and 2 (e.g. fastp_p1.fastq and fastp_p2.fastq), and these will be the input for the next step of aligning the reads to the reference genome.


## Step 2: Aligning to the reference genome using BWA and post alignment processing with SAMTools
So now that we have our quality trimmed reads, we will go ahead and align our reads to the reference human genome version Hg38, which is included in the GATK reference bundle (for more details and availability look here [https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle]).

As with all the analysis steps, there are multiple tools to accomplish the same task, and the same can be said of alignment. Our approach is to use BWA mem, especially when it comes to aligning vs Hg38 and if you want to learn more as to why, have a look here []

# VINU stuff

