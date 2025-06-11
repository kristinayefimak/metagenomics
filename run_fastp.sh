module load fastp

# Loop over all _1.gz files (forward reads)
for fwd in data/*_1.gz; do
    # Extract the base sample name (e.g., A1_AnOudin_C02)
    base=$(basename "$fwd" _1.gz)
    
    # Define reverse read path
    rev="data/${base}_2.gz"
    
    # Define output files
    out1="results/${base}_R1.clean.fq.gz"
    out2="results/${base}_R2.clean.fq.gz"
    html="results/${base}_fastp.html"
    json="results/${base}_fastp.json"
    
    echo "Running fastp on sample: $base"

    fastp \
      -i "$fwd" \
      -I "$rev" \
      -o "$out1" \
      -O "$out2" \
      -h "$html" \
      -j "$json"
done

