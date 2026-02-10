import csv

def ddinter_integration_1(d):
    try:
        pred_path = f"{d}/ddinter_integration_1.csv"
        gold_path = "benchmark/gold_results/ddinter_integration_1.csv"

        def load_csv_as_records(path):
            with open(path, newline='', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                records = []
                for row in reader:
                    # 每一行转为排序后的 tuple
                    records.append(tuple(sorted(row.items())))
                return records

        pred_records = load_csv_as_records(pred_path)
        gold_records = load_csv_as_records(gold_path)

        # 行顺序不敏感
        return sorted(pred_records) == sorted(gold_records)

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
