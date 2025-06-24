#Download and index Daphnia reference genome:

#STEP1: Download the Daphnia reference genome= ONLY DO THIS THE FIRST TIME=> from then on this file is present in directory
#You need a .fna file (FASTA format of the genome). For example, from NCBI:
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/631/705/GCF_020631705.1_ASM2063170v1.1/GCF_020631705.1_ASM2063170v1.1_genomic.fna.gz
gunzip GCF_020631705.1_ASM2063170v1.1_genomic.fna.gz

module load Bowtie2

#Step 2: Build the Bowtie2 index= this already takes into account both forward and reverse strand= 
#Bowtie2 creates 6 index files — these are binary-encoded versions of your genome that make alignment very fast and memory-efficient.
bowtie2-build GCF_020631705.1_ASM2063170v1.1_genomic.fna daphnia_index

#STEP 3: Map reads and extract non-host pairs= Reads that align well are likely Daphnia and should be removed.
#Reads that fail to align are likely microbiome and are retained.

#INSTALL BOWTIE2 AND SAMTOOLS WITH CONDA IN AN ENVIRONMENT (only do this once!!!)
conda create -n bowsam -c bioconda bowtie2 samtools


#PASTE ALL THIS in job composer:

#!/bin/bash
#SBATCH --job-name=remove_host
#SBATCH --output=logs/remove_host_%j.out
#SBATCH --error=logs/remove_host_%j.err

#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kristina.yefimak@kuleuven.be
#SBATCH -A lp_svbelleghem

# Load modules
# Load required toolchain first
source ~/.bashrc

conda activate bowsam
cd /scratch/leuven/363/vsc36367/metagenomics/


# Run your host-removal script

INDEX="daphnia_index"
mkdir -p results/nonhost 
for fwd in results/*_R1.clean.fq.gz; do
    
    base=$(basename "$fwd" _R1.clean.fq.gz)
    
    
    rev="results/${base}_R2.clean.fq.gz"
    
    sam="results/${base}.host.sam" bam="results/nonhost/${base}_nonhost.bam" 
    out1="results/nonhost/${base}_R1.nonhost.fq" out2="results/nonhost/${base}_R2.nonhost.fq" echo "Processing sample: 
    $base"
    
    bowtie2 -x $INDEX -1 "$fwd" -2 "$rev" \ --very-sensitive -S "$sam"
    samtools view -b -f 12 -F 256 "$sam" > "$bam"
    
    samtools fastq -1 "$out1" -2 "$out2" "$bam"
    
    
    rm "$sam"
done

#UNTIL HERE


#STOP PASTE HERE

#data is stored in RESULTs/nonhost
