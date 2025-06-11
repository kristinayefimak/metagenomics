mkdir -p data results scripts envs

#transfer data on globus to /scratch/363/vsc36367/metagenomics/data/ so all data is in data folder

#prepare VSC

module load Anaconda3
conda create -n metagenome fastp bowtie2 samtools kraken2 -c bioconda -y
conda activate metagenome

