import csv

def asmdb_annotate(d):
    try:
        pred_path = f"{d}/asmdb_annotate.csv"
        gold_path = "benchmark/gold_results/asmdb_annotate.csv"

        def load_csv_as_set(path):
            with open(path, newline='') as f:
                reader = csv.DictReader(f)
                fieldnames = reader.fieldnames
                rows = set()
                for row in reader:
                    # 保持大小写、值完全一致
                    rows.add(tuple((k, row[k]) for k in fieldnames))
            return set(fieldnames), rows

        pred_fields, pred_rows = load_csv_as_set(pred_path)
        gold_fields, gold_rows = load_csv_as_set(gold_path)

        # 列名必须完全一致（不区分顺序）
        if pred_fields != gold_fields:
            return False

        # 行内容完全一致（不区分顺序）
        return pred_rows == gold_rows

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
