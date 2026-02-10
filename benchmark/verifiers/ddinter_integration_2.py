import csv

def ddinter_integration_2(d):
    try:
        pred_path = f"{d}/ddinter_integration_2.csv"
        gold_path = "benchmark/gold_results/ddinter_integration_2.csv"

        def load_csv_as_set(path):
            with open(path, newline='', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                # Normalize each row into a frozenset of (key, value) pairs
                rows = set()
                for row in reader:
                    rows.add(frozenset(row.items()))
                return rows

        pred_rows = load_csv_as_set(pred_path)
        gold_rows = load_csv_as_set(gold_path)

        return pred_rows == gold_rows

    except Exception as e:
        print(f"Error processing CSV files: {e}")
        return False
