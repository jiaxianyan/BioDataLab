fastq-dump SRR31029759 --split-files -N 1 -X 1000000 -O toy_srr

REF_DIR=~/biodatalab/benchmark/dataset/CyanoOmicsDB/ncbi_dataset/data/GCF_000009725.1
READS_DIR=~/biodatalab/benchmark/dataset/CyanoOmicsDB/toy_srr

# 命令格式: bowtie2-build [参考基因组fasta] [索引文件前缀]
bowtie2-build $REF_DIR/GCF_000009725.1_ASM972v1_genomic.fna PCC6803_index

bowtie2 -x PCC6803_index \
  -1 $READS_DIR/SRR31029759_1.fastq \
  -2 $READS_DIR/SRR31029759_2.fastq \
  -S SRR31029759.sam \
  -p 8

samtools sort -n -O bam -o SRR31029759_nsorted.bam SRR31029759.sam

head $REF_DIR/genomic.gff

htseq-count -f bam -r name -s no -t gene -i locus_tag \
  SRR31029759_nsorted.bam \
  $REF_DIR/genomic.gff \
  > gene_counts.txt