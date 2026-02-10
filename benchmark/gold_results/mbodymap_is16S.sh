cd ~/biodatalab/benchmark/dataset/mBodyMap
mkdir -p metadata_info
cd ~/biodatalab/benchmark/gold_programs

RUN=SRR1521254 #$1
THREADS=4
# 查询从SRA中下载的序列的metainfo
python check_sra_strategy.py SRR1521254
# vdb-dump --info SRR758663

echo "SILVA blast test..."
#fasterq-dump $RUN --threads $THREADS --split-files
RUNQ="$RUN.fastq"
RUNA="$RUN.fasta"

seqtk seq -a $RUNQ > $RUNA
blastn -query $RUNA \
  -db ./mapseq_db/SILVA_138.2_SSURef_NR99_tax_silva.fasta \
  -outfmt 6 \
  -max_target_seqs 1 \
  -num_threads $THREADS \
  > hits.txt

HITS=$(wc -l < hits.txt)
echo "SILVA hit reads: $HITS / 1000"

if [ $HITS -gt 700 ]; then
  echo "Conclusion: This is very likely 16S amplicon data."
elif [ $HITS -lt 100 ]; then
  echo "Conclusion: This is likely shotgun metagenome or RNA-seq."
else
  echo "Conclusion: Mixed or unclear dataset."
fi