fastq-dump SRR2131217 --split-files -N 1 -X 100000 -O toy_srr1
fastq-dump SRR2131222 --split-files -N 1 -X 100000 -O toy_srr2
fastq-dump SRR2131244 --split-files -N 1 -X 100000 -O toy_srr3

fastq-dump SRR2131217 --split-files -N 1 -X 100000 -O toy_srr4
fastq-dump SRR2131217 --split-files -N 1 -X 100000 -O toy_srr5

hisat2-build GCF_001039355.2_LinAna2.0_genomic.fna genome_index

hisat2 \
  -x genome_index \
  -1 /root/biodatalab/benchmark/dataset/MetazExp/toy_srr1/SRR2131217_1.fastq \
  -2 /root/biodatalab/benchmark/dataset/MetazExp/toy_srr1/SRR2131217_2.fastq \
  --dta \
  -p 4 \
  -S toy.sam

hisat2 \
  -x genome_index \
  -1 /root/biodatalab/benchmark/dataset/MetazExp/toy_srr2/SRR2131222_1.fastq \
  -2 /root/biodatalab/benchmark/dataset/MetazExp/toy_srr2/SRR2131222_2.fastq \
  --dta \
  -p 4 \
  -S toy.sam

hisat2 \
  -x genome_index \
  -1 /root/biodatalab/benchmark/dataset/MetazExp/toy_srr3/SRR2131244_1.fastq \
  -2 /root/biodatalab/benchmark/dataset/MetazExp/toy_srr3/SRR2131244_2.fastq \
  --dta \
  -p 4 \
  -S toy.sam