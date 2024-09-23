# Research-Transcriptomics-Analysis-Workflows

# RNA-seq Workflow Setup.
## Perform installations for all the provided packages in your WSL or Ubuntu Linux Local system to ensure smooth running.

```bash
# Update system and packages
sudo apt update

# Install FASTQC for quality analysis of FASTQ File
sudo apt install fastqc  

# Install Trimmomatic for trimming low-quality reads from FASTQ files, which result from poor sequencing signal
sudo apt-get install -y trimmomatic  

# Install HISAT2 for aligning reads to the reference genome
sudo apt-get -y install hisat2  

# Install FeatureCounts (subread package) to compute gene counts based on ENSEMBL IDs
sudo apt install subread  

# Install unzip utility for extracting compressed files
sudo apt install unzip  

# Create a directory for storing files during this workflow
mkdir Reads\ to\ Counts 

# Note: Store your FASTQ file in this folder
