cd ~/biodatalab/benchmark/dataset/DIANA-miTED
# random select a SRA data point of DIANA-miTED database and down thhe origin file from NCBI SRA, here, we select ERR2731324
zcat ERR2731324.fastq.gz | head -n 400000 | gzip > sample.100k.fastq.gz
mkdir fastqc_raw
fastqc sample.100k.fastq.gz -o fastqc_raw
dnapi.py sample.100k.fastq.gz > adapter.txt
cutadapt \
  -a $(cat adapter.txt) \
  -q 10 \
  -m 18 \
  -o sample.100k.trimmed.fastq.gz \
  sample.100k.fastq.gz

bowtie-build chr22.fa chr22

bowtie chr22 \
  -q sample.100k.trimmed.fastq.gz \
  -S sample.100k_vs_chr22.sam \
  -v 1 \
  -k 1 \
  --best --strata

samtools view -c sample.100k_vs_chr22.sam
samtools view -c -F 4 sample.100k_vs_chr22.sam


samtools view -bS sample.100k_vs_chr22.sam > sample.bam
samtools sort sample.bam -o sample.sorted.bam
samtools view sample.sorted.bam > sample.sorted.sam

convert_bowtie_output.pl sample.sorted.sam > sample.arf

zcat sample.100k.trimmed.fastq.gz \
| awk 'NR%4==1{h=$1; sub(/^@/,">",h); print h} NR%4==2{print}' \
> sample.100k.trimmed.clean.fa

mapper.pl sample.100k.trimmed.clean.fa \
  -c -m \
  -s sample.100k.trimmed.collapsed.fa
  
mkdir -p quant_out

quantifier.pl \
  -p miRBase_v22/hairpin.fa \
  -m miRBase_v22/mature.fa \
  -r sample.100k.trimmed.collapsed.fa \
  -y demo_sample \
  -t hsa