# ---------------------------------------------------------------------------------------------------------------------------------------

# References and Citations 

# FASTQC - Andrews, S. (2010). FastQC:  A Quality Control Tool for High Throughput Sequence Data [Online]. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Trimmomatic - Trimmomatic: a flexible trimmer for Illumina sequence data - https://doi.org/10.1093/bioinformatics/btu170

# HISAT2 - https://daehwankimlab.github.io/hisat2/

# Subread - Liao Y, Smyth GK and Shi W (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923-30.

# ----------------------------------------------------------------------------------------------------------------------------------------

# Script Developed and Recreated by Rajarshi Ray , 29-09-2024 , Tampere University 

# STEP 1 : Set up Linux WSL or local system Linux

sudo apt update # Update system and packages

sudo apt install fastqc  # Install FASTQC for quality analysis of FASTQ File
sudo apt-get install -y trimmomatic # Install Trimmomatic for trimming low quality reads from FASTQ file which result due to poor signalling
sudo apt-get -y install hisat2  # Install Hisat2 for read aligning on reference genome
sudo apt install subread  # Package containing Featurecounts module to compute genes based on ENSEMBL ID to observe which genes are found.
sudo apt install unzip  # For extracting zip files

mkdir Reads\ to\ Counts # Default repository for files storing as per this workflow 
# PS - Store your fastq file in this folder only

# STEP 2 : Download necessary biological files for analysis 

wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz  # Download the Genome indices for HISAT2 for alignment


wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz  # Download Gene Annotation for Human GRCh38

# STEP 3 : Use FASTQC to analyse the reads and base quality based on sequencer output

fastqc demo.fastq  # Stores output as html file for further analysis

# STEP 4 : Perform trimming of FASTQ reads using Trimmomatic tools to screen out low filter reads.
## RULES
### 1. SE represents Single End sequencing , use PE if the sequencing is done in Pair End basis

java -jar /usr/share/java/trimmomatic-0.39.jar SE -threads 4 Reads\ to\ Counts/demo.fastq Reads\ to\ Counts/demo_trimmed.fastq TRAILING:10 -phred33      

# Re run STEP 3 again for cross checking the fastq file after trimming to see the analysed quality.

# STEP 5 : Perform alignment of reads onto human reference genome using HISAT2.

tar -xzf grch38_genome.tar.gz   # Use this

## RULES
### 1. In the rna strandness R is Reverse stranded sequencing from SE, if the sequencing is done in SE Forward use F, if PE , Read1 comes from F and Read2 comes from R , use FR, in case of PE Read1 from R and Read2 from F, use FR.

hisat2 -q --rna-strandness R -x grch38/genome -U demo_trimmed.fastq | samtools sort -o demo_trimmed.bam # Align the FASTQ Reads on the reference genome. 

# STEP 6 : Perform gene annotation and counts of reads per gene based on gtf gene annotation files

gzip -d Homo_sapiens.GRCh38.106.gtf.gz  # Unzip the Gene Annotation file 

# P.S - If yu have multiple reads for which you have generated BAM/SAN files, in the below code, add the BAM files space separated for DESEQ Analysis

featureCounts -S 2 -a Homo_sapiens.GRCh38.106.gtf -o demo_counts.txt demo_trimmed.bam  # Perform the gene annotation and assess the transcript count per ENSEMBLE Entry genes.

# STEP 7 : Process the Featurecounts file - Use tab delimited setting to view the output file properly.
