cd ~/biodatalab/benchmark/dataset/mBodyMap
# 输入文件 (必须是 FASTA 格式)
INPUT_FASTA="SRR3189676_clean_R.fasta"
SAMPLE_ID="Sample_16S"
OUTPUT_DIR="taxonomy_result"
THREADS=4

# 数据类型: "16S" 或 "META" (宏基因组)
DATA_TYPE="16S"

# MAPseq 数据库 16S选择SILVA SSURef数据库
MAPSEQ_DB="./mapseq_db/SILVA_119_SSURef_Nr99_tax_silva.fasta"
# MetaPhlAn2 数据库索引通常由工具自动管理，但需确保环境变量已配置
mkdir -p $OUTPUT_DIR
# 初始化 QC 状态 (1=Pass, 0=Fail)
QC_STATUS=1
QC_REASON="Passed"
echo ">>> Processing Sample: $SAMPLE_ID | Type: $DATA_TYPE"

# 丰度计算
# 用于存储最终丰度表的文件名 (格式: Taxon \t Relative_Abundance_%)
ABUNDANCE_TABLE="${OUTPUT_DIR}/${SAMPLE_ID}_abundance.txt"

if [ "$DATA_TYPE" == "16S" ]; then
    echo ">>> Running MAPseq v1.2 for 16S..."

    MSEQ_OUT="${OUTPUT_DIR}/${SAMPLE_ID}.mseq"
    mapseq "$INPUT_FASTA" \
        -nthreads $THREADS \
        "$MAPSEQ_DB" > "$MSEQ_OUT"
    echo ">>> Calculating Relative Abundances (Filter: Genus score >= 0.4)..."
    
    # 注意: MAPseq 输出列可能随版本变化，以下基于常见格式：
    # 通常第14列附近是 score，具体列需根据 mapseq -show-header 确认
    # 这里假设输出包含表头，使用 awk 处理
    
    awk -v out_file="$ABUNDANCE_TABLE" '
    BEGIN { FS="\t"; OFS="\t" }
    /^#/ { next } # 跳过注释行
    
    {
        # 假设 MAPseq 输出格式包含分类层级和分数
        # $0 是整行。我们需要根据实际输出提取 Genus, Species 和 Score
        # 这是一个通用逻辑示例，假设第12列是Genus，第13列是Species，第14列是Genus Score
        
        # 实际操作中通常需要解析详细的 mapseq 输出结构
        # 这里模拟逻辑：如果 Genus Score >= 0.4，则纳入统计
        
        # 提取相关字段 (根据实际情况调整列号)
        genus = $12 
        species = $13
        score = $14
        
        if (score >= 0.4) {
            total_valid++
            if (genus != "") counts_genus[genus]++
            if (species != "") counts_species[species]++
        }
    }
    END {
        print "Level\tTaxon\tRelative_Abundance" > out_file
        # 输出 Genus 丰度
        for (g in counts_genus) {
            abundance = (counts_genus[g] / total_valid) * 100
            print "Genus", g, abundance >> out_file
        }
        # 输出 Species 丰度
        for (s in counts_species) {
            abundance = (counts_species[s] / total_valid) * 100
            print "Species", s, abundance >> out_file
        }
    }
    ' "$MSEQ_OUT"

elif [ "$DATA_TYPE" == "META" ]; then
    echo ">>> Running MetaPhlAn2 for Metagenomics..."
    
    # 运行 MetaPhlAn2
    # 输入类型指定为 fasta (--input_type fasta)
    # 默认参数计算相对丰度
    METAPHLAN_OUT="${OUTPUT_DIR}/${SAMPLE_ID}_metaphlan.txt"
    
    metaphlan2.py "$INPUT_FASTA" \
        --input_type fasta \
        --nproc $THREADS \
        -o "$METAPHLAN_OUT"
    
    # 格式化输出以匹配 QC 检查
    # MetaPhlAn2 输出已经包含相对丰度 (0-100)
    # 我们只需筛选出 Genus (g__) 和 Species (s__) 层级
    
    echo "Level\tTaxon\tRelative_Abundance" > "$ABUNDANCE_TABLE"
    
    grep "s__" "$METAPHLAN_OUT" | grep -v "t__" | awk -F '\t' '{print "Species", $1, $2}' >> "$ABUNDANCE_TABLE"
    grep "g__" "$METAPHLAN_OUT" | grep -v "s__" | awk -F '\t' '{print "Genus", $1, $2}' >> "$ABUNDANCE_TABLE"
    
fi

# 检查生成的丰度表中是否存在 >= 99.99 的值
# 忽略表头 (NR>1)
MAX_ABUNDANCE=$(awk 'NR>1 {if ($3 > max) max=$3} END {print max}' "$ABUNDANCE_TABLE")

# 使用 bc 进行浮点数比较
IS_DOMINANT=$(echo "$MAX_ABUNDANCE >= 99.99" | bc -l)

if [ "$IS_DOMINANT" -eq 1 ]; then
    QC_STATUS=0
    QC_REASON="Failed: Single taxon dominance detected (${MAX_ABUNDANCE}%)"
    echo "WARNING: $QC_REASON"
else
    echo "Pass: Max dominance is ${MAX_ABUNDANCE}%"
fi

echo "------------------------------------------------"
echo "Final QC Status for $SAMPLE_ID: $QC_STATUS"
echo "Reason: $QC_REASON"
echo "------------------------------------------------"

# 将 QC 结果写入文件
echo "SampleID: $SAMPLE_ID" > "${OUTPUT_DIR}/qc_report.txt"
echo "QC_Status: $QC_STATUS" >> "${OUTPUT_DIR}/qc_report.txt"
echo "Reason: $QC_REASON" >> "${OUTPUT_DIR}/qc_report.txt"