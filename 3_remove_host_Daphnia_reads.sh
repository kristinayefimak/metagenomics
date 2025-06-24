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

#make a file
nano remove_host_reads.sh

#PASTE ALL THIS:
#!/bin/bash
# Remove Daphnia host reads from cleaned FASTQ files using Bowtie2 + SAMtools


# Path to Bowtie2 index (built from Daphnia .fna)
INDEX="daphnia_index"  # Adjust if located elsewhere

# Make sure output folder exists
mkdir -p results/nonhost

# Loop through all cleaned forward reads
for fwd in results/*_R1.clean.fq.gz; do
    # Extract base sample name (e.g., A1_AnOudin_C02)
    base=$(basename "$fwd" _R1.clean.fq.gz)
    
    # Define reverse read
    rev="results/${base}_R2.clean.fq.gz"

    # Define intermediate/output file names
    sam="results/${base}.host.sam"
    bam="results/nonhost/${base}_nonhost.bam"
    out1="results/nonhost/${base}_R1.nonhost.fq"
    out2="results/nonhost/${base}_R2.nonhost.fq"

    echo "Processing sample: $base"

    # Align to Daphnia reference genome
    bowtie2 -x $INDEX -1 "$fwd" -2 "$rev" \
      --very-sensitive -S "$sam"

#Step 4: Extract the non-host reads
#This gives you non-Daphnia read pairs → ready for microbiome profiling with kraken2, metaphlan, etc.
    # Extract non-host reads (read pairs where neither mate maps)
    samtools view -b -f 12 -F 256 "$sam" > "$bam"

    # Convert BAM to FASTQ
    samtools fastq -1 "$out1" -2 "$out2" "$bam"
    
    # Optional: clean up SAM file to save space
    rm "$sam"
done
#UNTILL HERE PASTE IN NANO
#make it executable
#copy this file to the job you are gonna submit e.g.:

cp /scratch/leuven/363/vsc36367/metagenomics/remove_host_reads.sh /data/leuven/363/vsc36367/.ondemand/data/projects/default/4/
ls /data/leuven/363/vsc36367/.ondemand/data/projects/default/4/


#RUN THIS JOB ON THE CLUSTER
#on vsc: submit job with this script:

#!/bin/bash
#SBATCH --job-name=remove_host
#SBATCH --output=logs/remove_host_%j.out
#SBATCH --error=logs/remove_host_%j.err

#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --cluster=wice
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kristina.yefimak@kuleuven.be
#SBATCH -A lp_svbelleghem

# Load modules
# Load required toolchain first

module load Bowtie2/2.4.4-GCC-10.3.0
module load SAMtools/1.18-GCC-12.3.0


# Run your host-removal script
bash remove_host_reads.sh




# Run your host-removal script
bash remove_host_reads.sh

#STOP PASTE HERE

#data is stored in RESULTs/nonhost
