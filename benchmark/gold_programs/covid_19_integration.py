import pandas as pd
from io import StringIO
import json

def extract_matrix_from_series_matrix_text(text: str) -> pd.DataFrame:
    start_tag = "!series_matrix_table_begin"
    end_tag = "!series_matrix_table_end"

    start = text.find(start_tag)
    if start == -1:
        raise ValueError("找不到 !series_matrix_table_begin")

    # 从下一行开始读表格
    table_start = text.find("\n", start)
    if table_start == -1:
        raise ValueError("表格开始标记后没有换行")

    end = text.find(end_tag, table_start)
    if end == -1:
        # 你贴出来的片段可能没有包含 end tag，就读到文本末尾
        table_text = text[table_start:].strip()
    else:
        table_text = text[table_start:end].strip()

    # 读 TSV：第一列是 ID_REF
    df = pd.read_csv(StringIO(table_text), sep="\t", dtype={"ID_REF": str})
    df = df.set_index("ID_REF")

    # 转成数值
    df = df.apply(pd.to_numeric, errors="coerce")
    return df

def map_id_to_hgnc_symbol(df: pd.DataFrame, mapping: dict, how="keep_unmapped", agg="mean") -> pd.DataFrame:
    """
    how:
      - keep_unmapped: 没映射到的保留原ID
      - drop_unmapped: 没映射到的丢弃
    agg:
      - mean / median / max ... 用于多个探针映射到同一基因时合并
    """
    symbols = df.index.to_series().map(mapping)

    if how == "drop_unmapped":
        df2 = df.loc[symbols.notna()].copy()
        df2.index = symbols[symbols.notna()].values
    else:
        df2 = df.copy()
        # ================= 修改了这里 =================
        # 将 df2.index 转换为 Series 才能传给 fillna
        # 或者使用 symbols.where(symbols.notna(), df2.index)
        filled_symbols = symbols.fillna(df2.index.to_series())
        df2.index = filled_symbols.values
        # ============================================

    # 若多个ID映射到同一symbol，聚合合并
    if df2.index.duplicated().any():
        df2 = df2.groupby(df2.index).agg(agg)

    return df2
# ===== 使用示例 =====
with open('benchmark/dataset/COVID-19/id2hgnc.json', 'r') as f:
    mapping = json.load(f)

with open('benchmark/dataset/COVID-19/GSE153428/GSE153428_series_matrix.txt', 'r') as f:
    text = f.read()
expr = extract_matrix_from_series_matrix_text(text)
print(expr)
expr_hgnc = map_id_to_hgnc_symbol(expr, mapping, how="keep_unmapped", agg="mean")
expr_hgnc.to_csv("benchmark/dataset/COVID-19/GSE153428_expression_hgnc.csv")

print(expr_hgnc)