import os

def cancerscem_annotate(d):
    try:
        # 检查生成的文件是否存在
        return os.path.exists(f'{d}/cancerscem_annotate.json')
    
    except Exception as e:
        print(f"Error processing cancerscem annotation results: {e}")
        return False