

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
