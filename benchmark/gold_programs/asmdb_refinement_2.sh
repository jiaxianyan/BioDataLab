fastq-dump SRR942022 --split-files -N 1 -X 100000 -O rna_seq

fastp -i rna_seq/SRR942022_1.fastq \
      -o SRR942022_clean.fastq \
      -h SRR942022.html -j SRR942022.json

hisat2-build ref/TAIR10_chr_all.fas TAIR10_index

hisat2 -x TAIR10_index \
       -U SRR942022_clean.fastq \
       -S SRR942022.sam

samtools sort -o SRR942022_sorted.bam SRR942022.sam