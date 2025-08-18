# Download Miniconda
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run Miniconda
bash Miniconda3-latest-Linux-x86_64.sh

#List to see if installation is successful 
conda list 

#If not 
export PATH=~/miniconda3/bin:$PATH

#Import conda environment from yml file that you have unzipped and brought to the main directory
conda env create -n MOOC --file MOOC.yml

#Activate MOOC
conda activate MOOC

#### WEEK 1 Analysis

#Download sequence data
fastq-dump --split-files ERR5743893 # Reference = MN908947.fa

#Make a directory to save fastQC outputs
mkdir -p QC_Reports 

#Quality Check using FastQC
fastqc ERR5743893_1.fastq ERR5743893_2.fastq --outdir QC_Reports

#Move to QC_report directory
cd QC_Reports

#Generate and HTML file to show report
multiqc . 

#Go back to main directory
cd ..

#Create a directory to store results from BWA-MEM
mkdir Mapping

#Copy the reference genome & samples to mapping 
cp MN908947.fasta Mapping
cp ERR5743893_2.fastq Mapping/
cp ERR5743893_1.fastq Mapping/

#Move to the new directory
cd Mapping

#Index reference genome
bwa index MN908947.fasta 

#Map our sample to the reference
bwa mem MN908947.fasta ERR5743893_1.fastq ERR5743893_2.fastq > ERR5743893.sam

#Convert to BAM using samtools to save on space
samtools view -S -b ERR5743893.sam > ERR5743893.bam
#S - input is a SAM file  b - output should be a BAM file

#Store in a sorted order (sorted according to order of mapping)
samtools sort -o ERR5743893.sorted.bam ERR5743893.bam

#Index sorted BAM file for ease of access
samtools index ERR5743893.sorted.bam

#Then visualize via IGV using your fasta and fasta.fai as a reference
samtools faidx MN908947.fasta #run if the fasta.fai is not present
#https://www.youtube.com/watch?v=die924mosh0 tutorial for IGV

#Use FreeBayes to identify variants
freebayes -f MN908947.fasta ERR5743893.sorted.bam  > ERR5743893.vcf

#Compress VCF file to take up less space
bgzip ERR5743893.vcf
tabix ERR5743893.vcf.gz

#Use bcftools to display variants
bcftools query -f '%TYPE\n' ERR5743893.vcf.gz | sort | uniq -c


######## WEEK 2

#Set up conda channels (if you haven't already)
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#Create an environment called nextflow and install nextflow in it
conda create --name nextflow nextflow

#Activate nextflow environment
conda activate nextflow

#Activate MOOC
conda activate MOOC
 
#Create a for loop to run fastq-dump on each accessions in the txt containing samples
for i in $(cat samples.txt);do fastq-dump --split-files $i;done

#compress fastq with gzip
gzip *.fastq

#Download a python script to create a CSV containing sample names and locations of files
wget -L https://raw.githubusercontent.com/nf-core/viralrecon/master/bin/fastq_dir_to_samplesheet.py

#run python script
python3 fastq_dir_to_samplesheet.py data samplesheet.csv -r1 _1.fastq.gz -r2 _2.fastq.gz

#activate nextflow
conda activate nextflow

#run nextflow
nextflow run nf-core/viralrecon -profile conda \
  --max_memory '6.GB' --max_cpus 1 \
  --input samplesheet.csv \
  --outdir results/viralrecon \
  --protocol amplicon \
  --genome 'MN908947.3' \
  --primer_set artic \
  --primer_set_version 3 \
  --skip_kraken2 \
  --skip_assembly \
  --skip_pangolin \
  --skip_nextclade \
  --skip_mosdepth \
  --skip_asciigenome \
  --platform illumina 

#Show the outputs of the commands
cd results/viralrecon
ls

#clean up after nextflow
du -sh work #check how much space the work directory is taking (this folder is present to help resume work in case the pipeline is not complete)

#remove work
rm -rf work


