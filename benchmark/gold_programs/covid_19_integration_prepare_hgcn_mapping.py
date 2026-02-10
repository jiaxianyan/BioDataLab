import pandas as pd
import re
import json

def extract_hgnc(text):
    if pd.isna(text) or text == "---":
        return None
    
    # 1. 首先按 '///' 分割不同的记录条目
    entries = text.split(' /// ')
    
    symbols = set() # 使用集合去重
    
    for entry in entries:
        # 2. 按 '//' 分割条目内的详细描述
        parts = entry.split(' // ')
        
        # 逻辑 A: 匹配描述部分中的 (SYMBOL)
        # 通常在第 3 个部分（索引为 2）
        for p in parts:
            # 寻找类似 (OR4F5) 或 (SAMD11) 的模式
            # 匹配逻辑：找左括号，中间是不包含空格和括号的字母/数字，紧跟右括号
            match = re.search(r'\(([A-Z0-9\-]{2,})\)', p)
            if match:
                symbols.add(match.group(1))
            
            # 逻辑 B: 针对某些行直接标明了 [Source:HGNC Symbol;Acc:HGNC:28706]
            # 这种情况下通常前面就是基因名，如 "sterile alpha motif domain containing 11 [Source:..."
            if "Source:HGNC Symbol" in p:
                # 这种格式比较特殊，Symbol 可能在描述文本的最末尾
                # 如果逻辑 A 没抓到，可以尝试从这里提取
                pass

    return ",".join(symbols) if symbols else None

# --- 处理流程 ---

# 假设你的文件是 tab 分隔的 (TSV)

with open('benchmark/dataset/COVID-19/GPL23159-184565.txt', 'r') as f:
    lines = f.read().strip().split('\n')[11:]
    
IDs = [line.split('\t')[0] for line in lines if line.strip()]
SPOT_IDs = [line.split('\t')[-1] for line in lines if line.strip()]

print(IDs[0])
print(SPOT_IDs[0])

# # 模拟读取刚才的数据
data = {
    'ID': IDs,
    'SPOT_ID' :SPOT_IDs
}
# 注意：实际文件中 SPOT_ID 是最后一列，且包含了你提供的那长串内容
df = pd.DataFrame(data)

# 应用提取函数
df['HGNC_Symbol'] = df['SPOT_ID'].apply(extract_hgnc)

print(df[['ID', 'HGNC_Symbol']])

id2hgcn = {}
for id, hgnc in zip(df['ID'], df['HGNC_Symbol']):
    if hgnc is not None:
        id2hgcn[id] = hgnc
        
with open('benchmark/dataset/COVID-19/id2hgnc.json', 'w') as f:
    json.dump(id2hgcn, f, indent=4)
    
    