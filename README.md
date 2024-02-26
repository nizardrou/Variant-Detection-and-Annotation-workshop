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
As you might've imagined, this sort of analysis involves multiple steps, and multiple tools. The tools and their links are provided below in case you want to run this analysis on your own setup. We are using the HPC, so we don't have to install any of them since they are already there!

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
Unfortunately, all of the information of this case file has been lost! It is therefore up to you to find out what disease is associated with your dataset.

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
    --extract fastp_p2_r2.fastq \
    -o qt_qc/read2qc/ \
    -t 28
```

Before running the script, let's take a few minutes to decipher what the flags that we have passed to the quality trimming tool FastP mean. To do that, run `fastp --help` on the command line (make sure that you have loaded or installed fastp first otherwise you will get a "command not found" error).

Can you find out why we instructed FastP to only use 16 CPUs (and not 28) even though our machine has 28 CPUs available?
```
module purge
module load all
module load gencore/2
module load fastp/0.20.1
fastp --help
```

Great! So now that we understand what the parameters mean, and how these commands are linked together (the output of the fastp quality trimmed reads are then passed to fastqc for checking), let's go ahead and run them. Here are a few options.

1. Running the command by submitting to the SLURM queue.
If you are within the directory containing the script "variant_detection.sh", simply type `sbatch variant_detection.sh`

2. An interactive session on the HPC.
To start an interactive session type "srun -n 28 -m  90G -t 4:00:00 --pty /bin/bash"
This will allocate a compute node to you, with 28 CPUs, 90GB of memory, and it will be available for 4 hours. Once the session starts, navigate to the directory containing your script (see above), and you can simply copy and paste the commands interactively.

3. If you are using your own setup.
On a standalone machine (not the HPC e.g. laptop), and as long as you have installed the required software, you can just navigate to your directory containing the script, change the script permissions to make it executable, and run it, 
```
chmod 755 variant_detection.sh
./variant_detection.sh
```

Please note that you will have to edit the script a bit to make it work on your own setup. First you will have to delete (or comment out) all of the lines in the begining that start with "module", since these are HPC specific. You will also have to tweak the CPU and memory limits to make sure they match your hardware setup.

Once these steps complete, you should have two folder, **qt_qc** and **raw_qc**, as well as the quality trimmed fastq files for read 1 and 2 (e.g. fastp_p1.fastq and fastp_p2.fastq), and these will be the input for the next step of aligning the reads to the reference genome.


## Step 2: Aligning to the reference genome using BWA and post alignment processing with SAMTools
So now that we have our quality trimmed reads, we will go ahead and align our reads to the reference human genome version Hg38, which is included in the GATK reference bundle (for more details and availability look here [https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle]).

As with all the analysis steps, there are multiple tools to accomplish the same task, and the same can be said of alignment. Our approach is to use BWA mem, because it has been shown to be very accurate and efficient in mapping sequences of low divergence, as well as the ability to handle alternative mapping contigs in the human reference Hg38.

The output of BWA mem is a SAM alignment file, and there is little to no benefit in working with the SAM file, the BAM formats is better (less space on disk as well as the ability to index it for fast processing). The SAM to BAM conversion, BAM sorting, and BAM indexing will be performed using SAMTools.

### Alignment
Open the variant_detection.sh file again, comment out (add hashtags) to the FastQC/FastP lines from the previous section, and uncomment (remove the hashtags) the BWA and SAMTools lines so that they look like this,
```
bwa mem \
    -t 28 \
    -M /scratch/Reference_Genomes/Public/Vertebrate_mammalian/Homo_sapiens/GATK_reference_bundle_hg38/Homo_sapiens_assembly38.fasta \
    fastp_p1_r1.fastq \
    fastp_p1_r2.fastq \
    > p1.sam

# Convert the BWA MEM SAM format to BAM using SAMTools
samtools view \
    -b \
    -@ 28 \
    p1.sam \
    -o p1.bam


# Delete the SAM file once the BAM has been produced
rm p1.sam


# Sort the BAM alignment file by coordinates (default) using SAMTools
samtools sort \
    -@ 28 \
    p1.bam \
    > p1.sorted.bam

```
As you can see, we supply the full path to the reference genome (after the -M flag), if you are running this on your own machine, then change the path to point to the location where you downloaded your reference.
SAMTools is one of those important Bioinformatics software packages that we use on an almost daily basis. I would highly recommend that you spend sometime familiarizing yourself with the package.

Before we go ahead and run these lines, let's decipher quickly what is happening.
1. We are aligning the quality trimmed reads to the Hg38 reference genome using "BWA mem" and redirecting the output to a SAM file.
2. We then use "SAMTools view" to convert the SAM to a BAM.
3. Since we have a BAM, we no longer need the SAM, so we delete it with the "rm" command.
4. Finally, we use "SAMTools sort" to **coordinate sort** the alignment (Can you find out what are the different types of sorting a BAM file?).


## Step 3: Post processing the BAM file with PICARD/GATK
As we discussed in the introduction (see slide deck), The Broad Best Practise guides offer a "blueprint" as to how variant calling is ideally performed (at least in human samples). Before calling the variants from the alignments, it is important to get a handle, as well as correct, some of the errors that are caused by the sequencing technologies as well as the sequencing library preparation steps. We also need to add some additional information to our BAM file headers so that the file can be correctly processed by the variant caller, as well as allowing us the option to combine it with other samples (if we were performed joint genotyping, trio analysis, GWAS etc.).

So, just like the previous section (alignment), let's start by openning our script using your prefered text editor, commenting the previous section (BWA and SAMTools), and uncomment the picard AddReadGroups, SAMTools index, picard MarkDuplicates, gatk BaseRecalibrator, and gatk ApplyBQSR. It should look like this,
```
# Add Read groups to the sorted BAM alignments using PICARD
picard -Xmx80g AddOrReplaceReadGroups \
    I=p1.sorted.bam \
    O=p1.rg.sorted.bam \
    RGID=1 \
    RGLB=1 \
    RGPL=ILLUMINA \
    RGPU=1 \
    RGSM=p1


# Mark duplicates (PCR/Optical) in the BAM alignments using PICARD
picard -Xmx80g MarkDuplicates \
    M=p1.metrics.txt \
    I=p1.rg.sorted.bam \
    O=p1.md.sorted.bam


# Run base quality score recalibrations using GATK's BQSR
gatk --java-options "-Xmx80G" BaseRecalibrator \
    -R /scratch/Reference_Genomes/Public/Vertebrate_mammalian/Homo_sapiens/GATK_reference_bundle_hg38/Homo_sapiens_assembly38.fasta \
    -I p1.md.sorted.bam \
    --known-sites /scratch/Reference_Genomes/Public/Vertebrate_mammalian/Homo_sapiens/GATK_reference_bundle_hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
    -O p1.recal_data.txt


# Apply the recalibrations from the earlier step and produce a new BAM alignment file using GATK's ApplyBQSR
gatk --java-options "-Xmx80G" ApplyBQSR \
    -R /scratch/Reference_Genomes/Public/Vertebrate_mammalian/Homo_sapiens/GATK_reference_bundle_hg38/Homo_sapiens_assembly38.fasta \
    -I p1.md.sorted.bam \
    --bqsr-recal-file p1.recal_data.txt \
    -O p1.rc.sorted.bam


# Index the recalibrated BAM file using SAMTools
samtools index \
    -@ 28 \
    p1.rc.sorted.bam
```

As you can see, at every step, we are producing a new BAM file, so it can become challenging trying to remember what is what! To avoid confusion, figure out a naming system that makes sense to you. For us, we tend to include the intials of the command that produced a given BAM (in addition to the sample name). E.g. "p1.rc.sorted.bam" refers to a recalibrated (rc) and coordinate sorted BAM file.

Let's quickly figure out what happened in the steps above.

### Adding Read Groups
We added some information to our BAM file using the picard AddOrReplaceReadGroups command (Look at the SAM/BAM specifications page here for more information [https://samtools.github.io/hts-specs/]).
- RGID=1 (This tag identifies which read group each read belongs to, so each read group's ID must be unique).
- RGLB=1 (This tag identifies the library that the read groups belong to, this is useful when specifiying the your library is PCR-free for example).
- RGPL=ILLUMINA (This tag specifies the sequencing platform/technology and can be useful in tunning specific systematic error profiles).
- RGPU=1 The PU holds three types of information, the flowcell barcode, the sequencing lane, and the sample barcode.
- RGSM=p1 This tag refers to the actual sample name (p1 in this case) and will uniquely identify the BAM alignments, especially when the BAMs have been merged.

Let's breifly look at how the BAM file looks before adding the read groups and after (you can use the command below to quickly inspect any BAM for the presence of read groups),

Extract the read groups from the BAM file before adding the read groups (in p1.sorted.bam for example) 
`samtools view -H p1.sorted.bam | grep '^@RG'`
You will see that the commands returns nothing, and that is because we did not add any read groups at that point.
If we were to run the command aggain only this time supplying the bam with the read groups added, we will get,
```
samtools view -H p2.rg.sorted.bam | grep '^@RG'
@RG	ID:1	LB:1	PL:ILLUMINA	SM:p2	PU:1
```

### Marking Duplicates
Sequencing technologies sometimes introduce what are know as "Systematic errors", that is, errors that are inherent in the technology. Luckily, these technologies have been around for years now, and we understand these errors (and their profiles), and hence, we can account for them and mitigate their impact.
There are two main parts to this step of correcting and detecting these artefacts, detecting duplicate, and accounting for systematic errors.
Duplicates that we are interested in are either PCR duplicates or Optical duplicates. PCR duplicates arise during the sample preperation and amplication stage before sequencing a library. Optical duplicates arise when a single amplification cluster is incorrectly detected as multiple cluster by the optics in the sequencing instrument.
The PICARD command that we use to mitigate these sequencing duplicates is called MarkDuplicates, and below is a description of the tool from the Broad Institute, which explains how the command works [https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard],

>_The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method). The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024._

Please note that by default, the tool only marks the duplicate reads, but does not remove them. This is because the PICAR/GATK software ecosystem is capable of reading these bitwise flags and interpreting them. If your downstream analysis tool(s) do not interpret these flags, then it would be better to "remove" the duplicates during this stage (by passing the appropriate --REMOVE_DUPLICATES flag). In addition to producing a BAM with the duplicates marked, MarkDuplicates also outputs a metrics file that summerizes the detection results.

### Base Quality Score Recalibration
GATK's BaseRecalibrator will then take the duplicate marked BAM file and will then attempt to detect (and correct) the second class of artefacts, systematic errors, and these are errors that the sequencing instrument makes when it detects or "calls" a sequenced nucleotide, as well as the Q score (or accuracy) associated with that call. Just to be clear though, BaseRecalibrator **DOES NOT** change the nucleotide calls, it adjusted the quality score associated with a given call.

Since we don't like to repeat what has already been very nicely described by the developers themselves, below is the explaination from the Broad Institute (and the complete explaination of tool is here [https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR]), 

>_Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. For example we can identify that, for a given run, whenever we called two A nucleotides in a row, the next base we called had a 1% higher rate of error. So any base call that comes after AA in a read should have its quality score reduced by 1%. We do that over several different covariates (mainly sequence context and position in read, or cycle) in a way that is additive. So the same base may have its quality score increased for one reason and decreased for another. 
This allows us to get more accurate base qualities overall, which in turn improves the accuracy of our variant calls. To be clear, we can't correct the base calls themselves, i.e. we can't determine whether that low-quality A should actually have been a T -- but we can at least tell the variant caller more accurately how far it can trust that A. Note that in some cases we may find that some bases should have a higher quality score, which allows us to rescue observations that otherwise may have been given less consideration than they deserve. Anecdotally our impression is that sequencers are more often over-confident than under-confident, but we do occasionally see runs from sequencers that seemed to suffer from low self-esteem. 
This procedure can be applied to BAM files containing data from any sequencing platform that outputs base quality scores on the expected scale. We have run it ourselves on data from several generations of Illumina, SOLiD, 454, Complete Genomics, and Pacific Biosciences sequencers._

Once BaseRecalibrator has finished and it has assessed our BAM file, we are then ready to apply the recalibrated scores and produce yet another BAM file (this is the last BAM file I promise!). We do this using GATK's ApplyBQSR tool, which takes as input the reaclibrated data table produced by BaseRecalibrator.

The last step here is to then index this final recalibrated BAM file using "samtools index" to allow coordinate based searching and faster processing time.


## Step 4: Alignment specific Quality Checking
Now that we have our finalized BAM file, it is a good idea to run some alignment specific QC. We like the output of Qualimaps bamqc as it produces graphical file (HTMLs/PDFs) as well as text based files. It is also versatile and can accomodate various scenarios such WGS and capture kits efficiently.
Please not that since this is a reduced dataset, some values might seem .... weird! More specifically, the coverage plots (since we are missing data from the whole genome and our reduced datasets only cover parts of a chromosome).
Similar to previous sections, let's comment out the previous sections and just uncomment the qualimap part, it should look like this (for p1),
```
qualimap --java-mem-size=80G \
    bamqc \
    -bam p1.rc.sorted.bam \
    -c \
    -gd HUMAN \
    -nt 28 \
    -outdir BAMQC \
    -outfile p1.bamqc.html \
    -outformat HTML \
    -sd


# Zip the GVCF
bgzip -@ 28 p1.g.vcf


# Index the zipped GVCF
tabix p1.g.vcf.gz
```
When you have done so, go ahead and execute the script again (either directly or through SLURM's sbatch). Can you figure out what the flags mean? The output should be under the "BAMQC" folder in your current directory.

For comparison, we have added the p1 complete dataset qualimap report in this Github repo. Go ahead and download it and see how it compares with your report.


## Step 5: Calling Variants (Finally!)
Well, this is why we are all here! To call variants and identify mutations in our data, which should help us relate them to the appropriate case files.

Again, there are multiple tools that can call variants, in fact, you will find that the best approaches, especially when dealing with clinical samples, is to use multiple variant callers. It might seem a bit redundant, but think of the implications, especially in a clinical diagnostic setting. This is in order to avoid what is know as "computational bias", since there is not single tool that can excel in every possible scenario, and using multiple software for the same task can help you mitigate that and identify "true positive" accurately.
In our case we will be using GATK's HaplotypeCaller, so go ahead and comment the previous qualimap section, and uncomment the HaplotypeCaller section, and run the script. It should look like this,
```
gatk --java-options "-Xmx80G" HaplotypeCaller \
    -R /scratch/Reference_Genomes/Public/Vertebrate_mammalian/Homo_sapiens/GATK_reference_bundle_hg38/Homo_sapiens_assembly38.fasta \
    -ERC GVCF \
    --dbsnp /scratch/Reference_Genomes/Public/Vertebrate_mammalian/Homo_sapiens/GATK_reference_bundle_hg38/Homo_sapiens_assembly38.dbsnp138.vcf \
    -I p1.rc.sorted.bam \
    -L chr7:80800000-85000000 \
    -O p1.g.vcf
```
Let's understand the flags.
- -R mentions the reference location (change this to reflect your setup similar to BWA mem earlier if you are running this one your own machine).
- -ERC GVCF instructs the caller to output a "Genomic VCF" not just a VCF. A GVCF also contains information on the locations where our sample is homozygous to the reference (the same) and not just the differences. This is especially useful for joint analysis and cohort based analysis, but it does result in larger file sizes.
- -I is the input BAM.
- --dbsnp provides the path to the DBSNP file for Hg38 and the tool will annotate any know SNPs that match the records in the file/variant.
- -L instructs the caller to only focus on the specified region. We do this in our case because of our reduced dataset, but if your data is WGS, you don't have to include this flag. You should use it if you have and exome or panel with a specific capture kit, which you can also supply as a BED file of regions.
- -O the name of the output GVCF file.

Finally, we want to zip the GVCF and index it using BGZIP and TABIX respectively in order to reduce the file sizes (GVCF can be huge, so always do this). The benefit of zipping and indexing is that the GVCF don't have to be uncompressed in order to processed by downstream tools, such as annotation tools.


## Step 6: Annotating our Variant GVCF file
The last step in our script is to annotate the variants that we discovered. We will be using two packages, SNPEff and SnpSift for this purpose. Go ahead and comment the previous sections (HaploTypeCaller, bgzip, and tabix), and uncomment the snpEff and SnpSift commands, and then run your script. They should look like this,
```
# Run SnpEff to annotate the variants effect in the GVCF
snpEff -Xmx80g eff \
    -nodownload \
    -csvStats p1.stats.csv \
    -o vcf \
    -s p1.summary.html \
    -lof GRCh38.105 p1.g.vcf.gz \
    > p1.snpEff.g.vcf


# Annotate the GVCF using SnpSift to add clinical variation (clinvar) annotations
SnpSift annotate \
    -tabix /scratch/Reference_Genomes/Public/SnpEff/GRCh38.105/clinvar.vcf.gz \
    p1.snpEff.g.vcf \
    > p1.snpEff.clinvar.g.vcf


# Annotate the GVCF using SnpSift to add dbSNP IDs to the variants
SnpSift annotate \
    -tabix /scratch/Reference_Genomes/Public/SnpEff/GRCh38.105/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    -dbsnp p1.snpEff.clinvar.g.vcf \
    > p1.snpEff.clinvar.dbsnp.g.vcf
```

SnpEff will assess the impact of the variants in our data. It will look at the changes caused by the short variants (SNPs and Indel) and predict the effect at the amino acid level and how the codons change. Similar to other tools, it needs a database for the organism/genome version in question, which we have already download and supplied beforehand. Have a look at the snpEff software page to understand how to set it up for your own particular case. SnpEff will produce another GVCF with the annotations now included, and you can certainly bgzip and index this file (although we are not doing that here).

SnpSift can also filter and annotate our SnpEff GVCF even further. We can supply multiple annotation levels to compare against. However, in our case, we are just providing the Clinical Variation database (Clinvar) and dbsnp.

Once these have finished running, you are ready to proceed to the results interpretation steps, which will be summerized below.



# VCF Annotation Filtering and Other Resources


## Tools for VCF Annotation

[SnpEff](http://snpeff.sourceforge.net)

[wAnnovar](http://wannovar.wglab.org)

[VEP](https://asia.ensembl.org/Tools/VEP)

[Franklin](https://franklin.genoox.com/home)

[Moon](http://www.diploid.com/moon)


## Variant Filtering

[snpsift](http://pcingola.github.io/SnpEff/snpsift/filter/)

```
cat Variant.vcf | SnpSift filter "(DP >= 20)) | (QUAL >= 30 )" > filtered.vcf 

```

[vcftools](https://vcftools.github.io/examples.html)

```
vcftools --gzvcf $VCF_IN  --remove-indels --max-missing $MISS \

--minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH \

--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout  >$VCF_OUT 

```
[gatk VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)

Spliting the VCF to SNPS only
```
gatk -T SelectVariants -R ref -V raw_variants.vcf -selectType SNP -o snps.vcf

```

Variant Filtering 
```
gatk -T VariantFiltration -R ref -V snps.vcf --filterExpression ' QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 ||SOR > 4.0' --filterName "basic_snp_filter" -o filtered_snps.vcf

```

In gatk variant filtering the variants are not filtered out but is labelled accordingly unlike snpsift makes subset(filtered) of the actual vcf


[bcftools](https://samtools.github.io/bcftools/bcftools.html#filter)



### Interpretation of Variants 

[Standard and Guidelines for Variant Interpretation](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4544753/)







