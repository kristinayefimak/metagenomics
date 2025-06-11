#Kraken2 is a k-mer–based classifier that: Breaks your reads into k-mers
#Compares them to a reference database of known sequences (e.g. bacteria, viruses, etc.)
#Assigns each read to the most specific taxonomic label based on where its k-mers match
#It’s fast and accurate, especially for microbiome reads.

#write code in file:
nano run_kraken2.sh

#!/bin/bash
# Batch classify all microbiome reads using Kraken2

module load kraken2

# Path to Kraken2 database
DB="/path/to/kraken2_db"  # Replace with actual path

# Create output folder if needed
mkdir -p results/kraken

# Loop over all non-host R1 reads
for fwd in results/nonhost/*_R1.nonhost.fq; do
    base=$(basename "$fwd" _R1.nonhost.fq)
    rev="results/nonhost/${base}_R2.nonhost.fq"

    report="results/kraken/${base}.report"
    output="results/kraken/${base}.kraken"

    echo "Running Kraken2 on sample: $base"

    kraken2 --db "$DB" \
            --paired "$fwd" "$rev" \
            --report "$report" \
            --output "$output"
done

SAVE FILE
bash run_kraken2.sh

#Make sure to replace "/path/to/kraken2_db" with your actual Kraken2 DB path (check with ls /path/to/kraken2_db/ to confirm hash.k2d etc. exist).
#Output files will go into results/kraken/ as:
#sample.report (summary table)
#sample.kraken (read-level assignments)
