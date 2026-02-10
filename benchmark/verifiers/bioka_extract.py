import csv

def bioka_extract(d):
    try:
        pred_path = f"{d}/bioka_extract.csv"
        gold_path = "benchmark/gold_results/bioka_extract.csv"

        def load_csv_as_set(path):
            with open(path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                rows = []
                for row in reader:
                    # 统一：key 和 value 都转成小写并去掉多余空格
                    normalized_row = {
                        k.strip().lower(): (v.strip().lower() if v is not None else "")
                        for k, v in row.items()
                    }
                    rows.append(frozenset(normalized_row.items()))
                return set(rows)

        pred_set = load_csv_as_set(pred_path)
        gold_set = load_csv_as_set(gold_path)

        return pred_set == gold_set

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
