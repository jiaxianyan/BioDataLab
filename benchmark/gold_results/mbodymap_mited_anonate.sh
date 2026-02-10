cd ~/biodatalab/benchmark/dataset/mBodyMap
# 下载的数据分类单端和双端，16S和宏基因，对应四种处理方式
# download the origin file from NCBI SRA, Project from mbodymap.	PRJNA477678-META-SRR7586637-2。   PRJNA453406-16S-SRR7061526-2
# PRJNA313226-16S-SRR3189676-1	PRJNA246028-META-SRR1521254_2-1

#fasterq-dump  --split-files
mkdir -p fastqc_raw

R=SRR1521254_2.fastq #单条序列修改
R1=SRR7061526_1.fastq #双条序列修改
R2=SRR7061526_2.fastq
THREADS=4
DATA_TYPE="META"  #META or 16S
DATA_NUM="SE" #SE or PE

if [ "$DATA_NUM" == "PE" ]; then
    fastqc $R1 -o fastqc_raw
    fastqc $R2 -o fastqc_raw
    # 动态计算原始读长的 2/3 用于过滤
    # 获取第一条read的长度 (去除换行符)
    # READ_LEN=$(zcat $R1 | head -n 2 | tail -n 1 | tr -d '\n' | wc -c) #gz文件处理
    READ_LEN=$(head -n 2 $R1 | tail -n 1 | tr -d '\n' | wc -c)
    # 计算 2/3 长度 (整数运算)
    MIN_LEN=$(( READ_LEN * 2 / 3 ))

    echo "Detected Read Length: $READ_LEN bp"
    echo "Minimum Length Threshold (2/3): $MIN_LEN bp"

    CLEAN_R1="${R1%.fastq}_clean_R1.fastq.gz"
    CLEAN_R2="${R2%.fastq}_clean_R2.fastq.gz"
    UNPAIRED_R1="${R1%.fastq}_unpaired_R1.fastq.gz"
    UNPAIRED_R2="${R2%.fastq}_unpaired_R2.fastq.gz"

    # 运行 Trimmomatic
    # ILLUMINACLIP: 去除接头
    # SLIDINGWINDOW:4:20: 去除低质量碱基 (窗口4，平均质量20)
    # MINLEN: 排除短于原始长度2/3的序列
    trimmomatic PE -threads $THREADS \
        $R1 $R2 \
        $CLEAN_R1 $UNPAIRED_R1 \
        $CLEAN_R2 $UNPAIRED_R2 \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:20 \
        MINLEN:$MIN_LEN

elif [ "$DATA_NUM" == "SE" ]; then
    fastqc $R -o fastqc_raw
    # READ_LEN=$(zcat $R | head -n 2 | tail -n 1 | tr -d '\n' | wc -c)  #命令处理gz文件
    READ_LEN=$(head -n 2 $R | tail -n 1 | tr -d '\n' | wc -c)
    # 计算 2/3 长度 (整数运算)
    MIN_LEN=$(( READ_LEN * 2 / 3 ))

    echo "Detected Read Length: $READ_LEN bp"
    echo "Minimum Length Threshold (2/3): $MIN_LEN bp"
    # 输出文件名定义
    CLEAN_R="${R%.fastq}_clean_R.fastq.gz"

    # 运行 Trimmomatic
    # ILLUMINACLIP: 去除接头
    # SLIDINGWINDOW:4:20: 去除低质量碱基 (窗口4，平均质量20)
    # MINLEN: 排除短于原始长度2/3的序列
    trimmomatic SE \
      -threads ${THREADS} \
      ${R} \
      -phred33 \
      ${CLEAN_R} \
      ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
      SLIDINGWINDOW:4:20 \
      MINLEN:${MIN_LEN}
fi

if [[ "$DATA_TYPE" == "16S" && "$DATA_NUM" == "PE" ]]; then
    # 使用 Casper (v0.8.2) 默认参数合并对端读段; vsearch代替
    # casper $CLEAN_R1 $CLEAN_R2 -o casper_merged \
    #     --out-prefix=clean_merged
    FINAL_FASTQ="${R1%.fastq}_16S_clean_merged.fastq"
    vsearch \
      --fastq_mergepairs $CLEAN_R1 \
      --reverse $CLEAN_R2 \
      --fastqout $FINAL_FASTQ
    CLEAN_DATA_FASTA="${FINAL_FASTQ%.fastq}.fasta"
    
    echo ">>> 16S Paired-end reads merged into: $FINAL_FASTQ"

elif [[ "$DATA_TYPE" == "META" && "$DATA_NUM" == "PE" ]]; then
    FINAL_FASTQ=$CLEAN_R1
    FINAL_FASTQ2=$CLEAN_R2
    CLEAN_DATA_FASTA="${FINAL_FASTQ%.fastq.gz}.fasta"
else
    FINAL_FASTQ=$CLEAN_R
    CLEAN_DATA_FASTA="${FINAL_FASTQ%.fastq.gz}.fasta"
fi

# ================= 4. 格式转换 (FASTQ -> FASTA) =================
# "如有必要，使用Seqtk将FASTQ序列转换为FASTA格式"
echo ">>> Converting FastQ($FINAL_FASTQ) to Fasta using Seqtk..."

# 如果是 16S，FINAL_FASTQ 是 vsearch 合并后的文件
# 如果是 META，FINAL_FASTQ 是 Trimmomatic 输出的文件
if [ -f "$FINAL_FASTQ" ] || [ -f "${FINAL_FASTQ%.gz}" ]; then
    if [[ "$DATA_TYPE" == "META" && "$DATA_NUM" == "PE" ]]; then
        CLEAN_DATA_FASTA2="${FINAL_FASTQ2%.fastq.gz}.fasta"
        seqtk seq -a $FINAL_FASTQ > $CLEAN_DATA_FASTA
        seqtk seq -a $FINAL_FASTQ2 > $CLEAN_DATA_FASTA2
    else
        seqtk seq -a $FINAL_FASTQ > $CLEAN_DATA_FASTA
        echo "Conversion complete."
    fi
else
    # 处理 .gz 压缩文件的情况
    seqtk seq -a <(zcat $FINAL_FASTQ) > $CLEAN_DATA_FASTA
    echo "Conversion complete (from gzipped input)."
fi

echo ">>> All steps finished. Clean Data is ready: $CLEAN_DATA_FASTA"