#HUMANN
#1-Prescreening: Identifies known microbial genomes from your reads
#2-Nucleotide-level alignment: Maps reads to known pangenomes (ChocoPhlAn)
#3-Translated search (optional): Maps unmatched reads to protein families (UniRef)
#4-Quantification:Outputs gene families: *_genefamilies.tsv+ Reconstructs metabolic pathways: *_pathabundance.tsv, *_pathcoverage.tsv

nano run_humann_all.sh
#PASTE THIS IN NANO:
#!/bin/bash
# Batch HUMAnN functional profiling for all non-host paired-end reads

module load HUMAnN/3.0.0  # or use your specific version

# Set output directory
mkdir -p results/humann_out

# Loop through all forward read files
for fwd in results/nonhost/*_R1.nonhost.fq; do
    base=$(basename "$fwd" _R1.nonhost.fq)
    rev="results/nonhost/${base}_R2.nonhost.fq"
    outfile="results/humann_out/${base}"

    echo "Running HUMAnN on: $base"

    humann \
        --input "${fwd},${rev}" \
        --output results/humann_out \
        --output-basename "$base"
done

#STOP PASTING
#run this as a job:

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
